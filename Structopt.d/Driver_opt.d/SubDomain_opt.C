#include <Structopt.d/Element_opt.d/Element_opt.h>
#include <Math.d/SparseMatrix.h>
#include <boost/scoped_array.hpp>

template<class Scalar>
double GenSubDomain_opt<Scalar>::getStrainEnergy(const Scalar* sol, int* eleList, int listSize)
{
  double energy=0.0;
  for(int iele=0; iele<numele; ++iele) 
    {       
      if (packedEset[iele]->isPhantomElement()) { continue; }
      if (!checkList(packedEset[iele]->getGlNum(),eleList,listSize)) { continue; }

      int numEleDOFs      = allDOFs->num(iele);
     
      getElDisp<Scalar>()->zero();
      getElKelD<Scalar>()->zero();
      
      FullSquareMatrix karray(numEleDOFs,kelData.get());
     
      for(int k=0; k<numEleDOFs; ++k) 
	{
	  int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
	  if(cn >= 0)	    
	    { (*getElDisp<Scalar>())[k] = sol[cn]; }
	  else
	    { (*getElDisp<Scalar>())[k] = this->Bcx((*allDOFs)[iele][k]); }
	}
      
      karray = packedEset[iele]->stiffness(nodes,karray.data());
      
      for (int i=0;i<numEleDOFs;i++) 
	{
	  for (int j=0;j<numEleDOFs;j++) 
	    {
	      (*getElKelD<Scalar>())[i]+=karray[i][j]*(*getElDisp<Scalar>())[j];
	    }
	}
      energy += 0.5*ScalarTypes::norm((*getElKelD<Scalar>())*(*getElDisp<Scalar>()));
    }
  return energy;
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenSubDomain_opt<Scalar>::getStructureMass(int* eleList, int listSize)
{
  double totmas = 0.0; // total mass of structure
  for(int iele=0; iele < numele; ++iele) 
    {
      if (packedEset[iele]->isPhantomElement()) { continue; }
      if (!checkList(packedEset[iele]->getGlNum(),eleList,listSize)) { continue; }
      totmas += packedEset[iele]->getMass(nodes);
    }

  double *nodeCount = static_cast<double *>(dbg_alloca(sizeof(double)*numnodes));
  double *nodeMass  = static_cast<double *>(dbg_alloca(sizeof(double)*numnodes));

  std::fill(nodeCount, nodeCount + numnodes, 0.0);
  std::fill(nodeMass,  nodeMass + numnodes, 0.0);
  // ADD DISCRETE MASS
  DMassData *current = firstDiMass;
  while(current != 0) 
    {
      int n = current->node;
      nodeMass[n]  += current->diMass;
      nodeCount[n] += 1;
      current = current->next;
    }

  for(int n=0; n<numnodes; ++n) 
    {
      if(nodeCount[n] == 0.0) { continue; }
      if(nodeCount[n] > 6) { nodeCount[n] = 6; }
      totmas += nodeMass[n]/nodeCount[n];
    }

  return totmas;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenSubDomain_opt<Scalar>::getGradduStrainEnergy(const Scalar* sol, Scalar* adj, 
						     int* eleList, int listSize)
{
  for(int iele=0; iele<numele; ++iele)
    {
      if (packedEset[iele]->isPhantomElement()) { continue; }
      if (!checkList(packedEset[iele]->getGlNum(),eleList,listSize)) { continue; }
      
      int NodesPerElement = elemToNode->num(iele);
      int numEleDOFs      = allDOFs->num(iele);
      
      FullSquareMatrix karray(numEleDOFs,kelData.get());      
      
      for(int k=0; k<numEleDOFs; ++k) 
	{
	  int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
	  if(cn >= 0) 
	    { (*getElDisp<Scalar>())[k] = sol[cn]; }
	  else 
	    { (*getElDisp<Scalar>())[k] = this->Bcx((*allDOFs)[iele][k]); }
	}
      
      karray=packedEset[iele]->stiffness(nodes,karray.data());
      getElAdj<Scalar>()->zero();
      
      karray.multiply(*getElDisp<Scalar>(), *getElAdj<Scalar>(), 1.0);
      
      for(int k=0; k<numEleDOFs; ++k) 
	{
	  int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
	  if(cn >= 0) 
	    { adj[cn] += ScalarTypes::conj((*getElAdj<Scalar>())[k]); }
	}
    }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenSubDomain_opt<Scalar>::getGradPartStrainEnergy(const Scalar* sol, int* eleList, int listSize)
{
  double gradenergy = 0.0;
  int numActelm = actInf->size();  
  for(int actelm=0; actelm < numActelm; ++actelm) 
    {
      int glbIele = actInf->getEntry(actelm);
      int locIele = this->globalToLocalElement(glbIele);
      if(locIele == -1) { continue; }
      if (packedEset[locIele]->isPhantomElement()) { continue; }
      if (!checkList(glbIele,eleList,listSize)) { continue; }

      int NodesPerElement = elemToNode->num(locIele);
      int numEleDOFs      = allDOFs->num(locIele);
     
      FullSquareMatrix gradkarray(numEleDOFs,gkelData.get());

      // get element displacements
      getElDisp<Scalar>()->zero();
      for (int i=0;i<numEleDOFs;i++) 
	{
	  int cn = c_dsa->getRCN((*allDOFs)[locIele][i]);
	  if(cn >= 0) 
	    { (*getElDisp<Scalar>())[i] = sol[cn]; }
	  else 
	    { (*getElDisp<Scalar>())[i] = this->Bcx((*allDOFs)[locIele][i]); }
	}    
      // get contribution of dKel*getElDisp<Scalar>()       
      dynamic_cast<Element_opt*>(packedEset[locIele])->gradstiffness(nodes,*gradnodes,gradkarray);
      getElKelD<Scalar>()->zero();
      gradkarray.multiply(*getElDisp<Scalar>(),*getElKelD<Scalar>(), 1.0);				   
      gradenergy += 0.5*ScalarTypes::Real((*getElKelD<Scalar>())*(*getElDisp<Scalar>()));
    }
  return gradenergy;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenSubDomain_opt<Scalar>::buildDRHSforceDs(Scalar* adj)
{
  // compute derivative of external force
  for(int i=0; i < numNeuman; ++i) 
    {
      int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
      if(dof < 0) { continue; }
      adj[dof] += gradnbc[i].val;
    }
  // add gravity forces
  if (gravityFlag()) 
    { fprintf(stderr,"WATCH: Derivatives of Gravity Force neglected\n"); }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenSubDomain_opt<Scalar>::DStiffnesDSmDUAdj(const Scalar* sol, const Scalar* adj)
{
  // building D(elastic stiffness) / D (optimization variable) 
  // * Displacements 
  // * Adjoint Vector
  double alpha = sinfo.alphaDamp;
  double  beta = sinfo.betaDamp;
  bool isShifted = geoSource->isShifted() || solInfo().isAcousticHelm();
  bool isDamped = (alpha != 0.0) || (beta != 0.0);
  
  double omega=0, omega2=0;  
  if(isShifted) 
    {
      omega2 = geoSource->shiftVal();
      omega  = sqrt(omega2);
    }

  int size = sizeof(double)*maxNumDOFs*maxNumDOFs;
  double *karray = static_cast<double *>(dbg_alloca(size));
  double *marray = static_cast<double *>(dbg_alloca(size));
  double *izarray = static_cast<double *>(dbg_alloca(size));
  double mratio = geoSource->getMRatio();

  Scalar val = 0.0;
  int numActelm = actInf->size(); 
  for(int actelm=0; actelm < numActelm; ++actelm) 
    {
      int glbIele = actInf->getEntry(actelm);
      int locIele = this->globalToLocalElement(glbIele);
      if(locIele == -1) { continue; }
      int NodesPerElement = elemToNode->num(locIele);
      int numEleDOFs      = allDOFs->num(locIele);
      FullSquareMatrix gradkarray(numEleDOFs,karray);
      FullSquareMatrix gradmarray(numEleDOFs,marray);
      // ... DETERMINE ELEMENT DISPLACEMENT VECTOR 
     for(int k=0; k<numEleDOFs; ++k) 
       {       
	 int cn = c_dsa->getRCN((*allDOFs)[locIele][k]);
	 if(cn >= 0) 
	   { (*getElDisp<Scalar>())[k] = sol[cn]; }
	 else 
	   { (*getElDisp<Scalar>())[k] = this->Bcx((*allDOFs)[locIele][k]); }
       }
      Element_opt* el = dynamic_cast<Element_opt*>(packedEset[locIele]);
      assert(el != 0);
      getElRHS<Scalar>()->zero();     
      el->gradstiffness(nodes, *gradnodes, gradkarray);
      gradkarray.multiply(*getElDisp<Scalar>(),*getElRHS<Scalar>(),-1.0);
      if(isShifted)
	{
	  el->gradMassMatrix(nodes, *gradnodes, gradmarray,mratio);
	  gradmarray.multiply(*getElDisp<Scalar>(),*getElRHS<Scalar>(),omega2);
	  if(el->hasDamping())
	    {
	      FullSquareMatrix gradcarray(numEleDOFs,izarray);
	      el->gradDampMatrix(nodes, *gradnodes, gradcarray);
	      gradcarray.multiply(*getElDisp<Scalar>(),*getElRHS<Scalar>(),DComplex(0,-omega));
	    }
	}
      if(isDamped)
	{
	  gradkarray.multiply(*getElDisp<Scalar>(),*getElRHS<Scalar>(),DComplex(0,-omega*beta));
	  gradmarray.multiply(*getElDisp<Scalar>(),*getElRHS<Scalar>(),DComplex(0,-omega*alpha));
	}

      for(int k=0; k<numEleDOFs; ++k) 
	{
	  int cn = c_dsa->getRCN((*allDOFs)[locIele][k]);
	  if(cn >= 0) { val += (*getElRHS<Scalar>())[k] * adj[cn]; }
	}
    }  
  return val;
}


//------------------------------------------------------------------------------
template<class Scalar>
double GenSubDomain_opt<Scalar>::getGradStructureMass(int* eleList, int listSize)
{
  double gradtotmas = 0.0;
  int numActelm = actInf->size();  
  for(int actelm=0; actelm < numActelm; ++actelm) 
    {
      int glbIele = actInf->getEntry(actelm);
      int locIele = this->globalToLocalElement(glbIele);
      if(locIele == -1) { continue; }
      if (packedEset[locIele]->isPhantomElement()) { continue; }
      if (!checkList(glbIele,eleList,listSize) ) continue;
      double gradelementMass = dynamic_cast<Element_opt*>(packedEset[locIele])->getGradMass(nodes,*gradnodes);
      gradtotmas += gradelementMass;
    }
  //
  //  WATCH: derivatives w.r.t. variable nodal mass not implemented
  //  
  return gradtotmas;
}


//------------------------------------------------------------------------------
template<class Scalar>
void GenSubDomain_opt<Scalar>::makeGlobalToLocalElementMap()
{
  for(int iele=0; iele<numele; ++iele)
    {
      glToLocalElement[packedEset[iele]->getGlNum()] = iele;
    }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
int GenSubDomain_opt<Scalar>::globalToLocalElement(int glNum) const
{
  std::map<int, int>::const_iterator i = glToLocalElement.find(glNum);
  return (i==glToLocalElement.end())? -1: i->second;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenSubDomain_opt<Scalar>::getMidPointMass(double& totmas, double* midPoint)
{
  totmas = 0;
  midPoint[0] = midPoint[1] = midPoint[2] = 0.0;
  int numActelm = actInf->size();
  double xyz[3];
  for(int actelm=0; actelm < numActelm; ++actelm) 
    {
      int glbIele  = actInf->getEntry(actelm);
      int locIele  = this->globalToLocalElement(glbIele);
      if(locIele == -1) { continue; }
      totmas      += packedEset[locIele]->getMass(nodes);
      packedEset[locIele]->getMidPoint(nodes, xyz);
   
      midPoint[0] += xyz[0];
      midPoint[1] += xyz[1];
      midPoint[2] += xyz[2];
    }
  totmas /= this->actInf->size();
  midPoint[0] /= this->actInf->size();
  midPoint[1] /= this->actInf->size();
  midPoint[2] /= this->actInf->size();
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenSubDomain_opt<Scalar>::buildDStiffnesDSmDU(Scalar* rhs, const Scalar* sol, int flag)
{
  // building D(elastic stiffness) / D (optimization variable) * Displacements
  // flag = 0 means loop over only active elements
  // flag = 1 means loop over all elements
  double alpha = sinfo.alphaDamp;
  double  beta = sinfo.betaDamp;
  bool isShifted = geoSource->isShifted() || solInfo().isAcousticHelm();
  bool isDamped = (alpha != 0.0) || (beta != 0.0);
  
  double omega=0, omega2=0;  
  if(isShifted) 
    {
      omega2 = geoSource->shiftVal();
      omega  = sqrt(omega2);
    }

  int size = sizeof(double)*maxNumDOFs*maxNumDOFs;
  double *karray = static_cast<double *>(dbg_alloca(size));
  double *marray = static_cast<double *>(dbg_alloca(size));
  double *izarray = static_cast<double *>(dbg_alloca(size));
  double mratio = geoSource->getMRatio();
  
  int numActelm = (flag==0)? actInf->size() : numele;
  for(int actelm=0; actelm < numActelm; ++actelm) 
    {
      int locIele, glbIele;
      if(flag==0)
	{
	  glbIele = actInf->getEntry(actelm);
	  locIele = this->globalToLocalElement(glbIele);
	  if(locIele == -1) { continue; }
	}
      else
	{
	  locIele = actelm;
	  glbIele = packedEset[locIele]->getGlNum();
	}
      int numEleDOFs      = allDOFs->num(locIele);
      FullSquareMatrix gradkarray(numEleDOFs,karray);
      FullSquareMatrix gradmarray(numEleDOFs,marray);
      // ... DETERMINE ELEMENT DISPLACEMENT VECTOR 
     for(int k=0; k<numEleDOFs; ++k) 
       {
	 int cn = c_dsa->getRCN((*allDOFs)[locIele][k]);
	 if(cn >= 0) 
	   { (*getElDisp<Scalar>())[k] = sol[cn]; }
	 else 
	   { (*getElDisp<Scalar>())[k] = this->Bcx((*allDOFs)[locIele][k]); }
       }
      Element_opt* el = dynamic_cast<Element_opt*>(packedEset[locIele]);
      assert(el != 0);
      getElRHS<Scalar>()->zero();     
      el->gradstiffness(nodes, *gradnodes, gradkarray);
      gradkarray.multiply(*getElDisp<Scalar>(),*getElRHS<Scalar>(),-1.0);
      if(isShifted)
	{
	  el->gradMassMatrix(nodes, *gradnodes, gradmarray,mratio);
	  gradmarray.multiply(*getElDisp<Scalar>(),*getElRHS<Scalar>(),omega2);
	  if(el->hasDamping())
	    {
	      FullSquareMatrix gradcarray(numEleDOFs,izarray);
	      el->gradDampMatrix(nodes, *gradnodes, gradcarray);
	      gradcarray.multiply(*getElDisp<Scalar>(),*getElRHS<Scalar>(),DComplex(0,-omega));
	    }
	}
      if(isDamped)
	{
	  gradkarray.multiply(*getElDisp<Scalar>(),*getElRHS<Scalar>(),DComplex(0,-omega*beta));
	  gradmarray.multiply(*getElDisp<Scalar>(),*getElRHS<Scalar>(),DComplex(0,-omega*alpha));
	}


     for(int k=0; k<numEleDOFs; ++k) 
       {
	 int cn = c_dsa->getRCN((*allDOFs)[locIele][k]);
	 if(cn >= 0) 
	   { rhs[cn] += (*getElRHS<Scalar>())[k]; }
       }
    }
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
double GenSubDomain_opt<Scalar>::getGradStrainEnergy(const Scalar* sol, const Scalar* grad,
						     int* eleList, int listSize)
{
  double gradenergy = 0.0;
  for(int iele=0; iele < numele; ++iele) 
    {
      if (packedEset[iele]->isPhantomElement()) { continue; }
      if (!checkList(packedEset[iele]->getGlNum(),eleList,listSize)) { continue; }
      int NodesPerElement = elemToNode->num(iele);
      int numEleDOFs      = allDOFs->num(iele);
      
      FullSquareMatrix karray(numEleDOFs,kelData.get());
      FullSquareMatrix gradkarray(numEleDOFs,gkelData.get());
      
      // get element displacements     
      getElDisp<Scalar>()->zero();
      for(int i=0; i<numEleDOFs; ++i) 
	{
	  int cn = c_dsa->getRCN((*allDOFs)[iele][i]);
	  if(cn >= 0) 
	    { (*getElDisp<Scalar>())[i] = sol[cn]; }
	  else
	    { (*getElDisp<Scalar>())[i] = this->Bcx((*allDOFs)[iele][i]); }
	}     
      // get contribution of Kel*elGrad
      getElKelD<Scalar>()->zero();
      getElGrad<Scalar>()->zero();     
      for(int i=0; i<numEleDOFs; ++i) 
	{
	  int cn = c_dsa->getRCN((*allDOFs)[iele][i]);
	  if(cn >= 0)
	    { (*getElGrad<Scalar>())[i] = grad[cn]; }
	}     
      karray = packedEset[iele]->stiffness(nodes,karray.data());
      karray.multiply(*getElGrad<Scalar>(),*getElKelD<Scalar>(),2.0);
      // get contribution of dKel*getElDisp<Scalar>()
      int actelm = actInf->check(packedEset[iele]->getGlNum());
      if(actelm) 
	{
	  dynamic_cast<Element_opt*>(packedEset[iele])->gradstiffness(nodes,*gradnodes,gradkarray);
	  gradkarray.multiply(*getElDisp<Scalar>(),*getElKelD<Scalar>(),1.0);
	}
      gradenergy += 0.5*ScalarTypes::Real((*getElKelD<Scalar>())*(*getElDisp<Scalar>()));
    }
  return gradenergy;
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenSubDomain_opt<Scalar>::getNodalDisp(const Scalar* sol, int glNode, int dispTyp, int weight, bool dblfy)
{
  int node = this->globalToLocal(glNode);
  if(node < 0) { return 0; }
  assert(sol != 0);
  int loc, loc1;
  Scalar val = 0, dummy;
  double nd, sd=0.0;

  switch(dispTyp) 
    {    
    case 0: 
      loc  = c_dsa->locate(node, DofSet::Xdisp);
      loc1 =   dsa->locate(node, DofSet::Xdisp);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = this->Bcx(loc1); }
      else                  { val = 0.0; }      
      break;
    case 1: 
      loc  = c_dsa->locate(node, DofSet::Ydisp);
      loc1 =   dsa->locate(node, DofSet::Ydisp);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = this->Bcx(loc1); }
      else                  { val = 0.0; }
      break;
    case 2: 
      loc  = c_dsa->locate(node, DofSet::Zdisp);
      loc1 =   dsa->locate(node, DofSet::Zdisp);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = this->Bcx(loc1); }
      else                  { val = 0.0; }
      break;
    case 3: 
      loc  = c_dsa->locate( node, DofSet::Xrot);
      loc1 =   dsa->locate( node, DofSet::Xrot);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = this->Bcx(loc1); }
      else                  { val = 0.0; }
      break;
    case 4: 
      loc  = c_dsa->locate( node, DofSet::Yrot);
      loc1 =   dsa->locate( node, DofSet::Yrot);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = this->Bcx(loc1); }
      else                  { val = 0.0; }
      break;
    case 5: 
      loc  = c_dsa->locate( node, DofSet::Zrot);
      loc1 =   dsa->locate( node, DofSet::Zrot);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = this->Bcx(loc1); }
      else                  { val = 0.0; }
      break;
    case 6:
      loc  = c_dsa->locate( node, DofSet::Xdisp);
      loc1 =   dsa->locate( node, DofSet::Xdisp);
      if      ( loc >= 0  ) { dummy = sol[loc];  } 
      else if ( loc1 >= 0 ) { dummy = this->Bcx(loc1); }
      else                  { dummy = 0.0; }
      nd = ScalarTypes::norm(dummy);
      sd += nd*nd;
      loc  = c_dsa->locate( node, DofSet::Ydisp);
      loc1 =   dsa->locate( node, DofSet::Ydisp);
      if      ( loc >= 0  ) { dummy = sol[loc];  } 
      else if ( loc1 >= 0 ) { dummy = this->Bcx(loc1); }
      else                  { dummy = 0.0; }
      nd = ScalarTypes::norm(dummy);
      sd += nd*nd;
      loc  = c_dsa->locate( node, DofSet::Zdisp);
      loc1 =   dsa->locate( node, DofSet::Zdisp);
      if      ( loc >= 0  ) { dummy = sol[loc];  } 
      else if ( loc1 >= 0 ) { dummy = this->Bcx(loc1); }
      else                  { dummy = 0.0; }
      nd = ScalarTypes::norm(dummy);
      sd += nd*nd;
      val = sqrt(sd);
      break;
    default:
      printf("OptcritDisp::evaluate - Wrong Disptype !");
      assert(0);
    }
  return (dblfy? ScalarTypes::doublify(val) : val)/static_cast<Scalar>(weight);
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenSubDomain_opt<Scalar>::getGradduNodalDisp(const Scalar* sol, Scalar* adj,
						  int glNode, int dispTyp, int weight,
						  bool dblfy)
{
  // NOTE: we add to adj here!!! It has to be initialized to zero elsewhere
  int node = this->globalToLocal(glNode);
  if(node < 0) { return; }
  assert(sol != 0 && adj != 0);
  int loc, locx, locy, locz, loc1x, loc1y, loc1z;  
  Scalar dummyx, dummyy, dummyz, val = 0;
  double nd, sd=0.0;
  
  switch(dispTyp) 
    {    
    case 0: 
      loc  = c_dsa->locate(node, DofSet::Xdisp);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>())/static_cast<Scalar>(weight); }
      break;
    case 1: 
      loc  = c_dsa->locate(node, DofSet::Ydisp);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>())/static_cast<Scalar>(weight); }
      break;
    case 2: 
      loc  = c_dsa->locate(node, DofSet::Zdisp);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>())/static_cast<Scalar>(weight); }
      break;
    case 3: 
      loc  = c_dsa->locate(node, DofSet::Xrot);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>())/static_cast<Scalar>(weight); }
      break;
    case 4: 
      loc  = c_dsa->locate(node, DofSet::Yrot);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>())/static_cast<Scalar>(weight); }
      break;
    case 5: 
      loc  = c_dsa->locate(node, DofSet::Zrot);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>())/static_cast<Scalar>(weight); }
      break;
    case 6:
      locx  = c_dsa->locate( node, DofSet::Xdisp);
      loc1x =   dsa->locate( node, DofSet::Xdisp);
      if      ( locx >= 0  ) { dummyx = sol[locx];  } 
      else if ( loc1x >= 0 ) { dummyx = this->Bcx(loc1x); }
      else                   { dummyx = 0.0; }
      nd = ScalarTypes::norm(dummyx);
      sd += nd*nd;
      locy  = c_dsa->locate( node, DofSet::Ydisp);
      loc1y =   dsa->locate( node, DofSet::Ydisp);
      if      ( locy >= 0  ) { dummyy = sol[locy];  } 
      else if ( loc1y >= 0 ) { dummyy = this->Bcx(loc1y); }
      else                   { dummyy = 0.0; }
      nd = ScalarTypes::norm(dummyy);
      sd += nd*nd;
      locz  = c_dsa->locate( node, DofSet::Zdisp);
      loc1z =   dsa->locate( node, DofSet::Zdisp);
      if      ( locz >= 0  ) { dummyz = sol[locz];  } 
      else if ( loc1z>= 0 ) { dummyz = this->Bcx(loc1z); }
      else                  { dummyz = 0.0; }
      nd = ScalarTypes::norm(dummyz);
      sd += nd*nd;
      val = sqrt(sd);
      if(sd > 0)
	{
	  if(locx>=0) { adj[locx] += ScalarTypes::conj(dummyx)/val/static_cast<Scalar>(weight); }
	  if(locy>=0) { adj[locy] += ScalarTypes::conj(dummyy)/val/static_cast<Scalar>(weight); }
	  if(locz>=0) { adj[locz] += ScalarTypes::conj(dummyz)/val/static_cast<Scalar>(weight); }
	}
      break;
    default:
      printf("OptcritDisp::evaluate - Wrong Disptype !");
      assert(0);
    }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenSubDomain_opt<Scalar>::getGradNodalDisp(const Scalar* sol, const Scalar* grad,
						  int glNode, int dispTyp, int weight,
						  bool dblfy)
{
  int node = this->globalToLocal(glNode);
  if(node < 0) { return 0; }
  assert(sol != 0 && grad != 0);
  int loc, locx, locy, locz, loc1x, loc1y, loc1z;
  Scalar val = 0;
  Scalar dummyx, dummyy, dummyz, ddx, ddy, ddz;
  double nd, sd=0.0;

  switch(dispTyp) 
    {    
    case 0: 
      loc  = c_dsa->locate(node, DofSet::Xdisp);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc])/static_cast<Scalar>(weight); } 
      break;
    case 1: 
      loc  = c_dsa->locate(node, DofSet::Ydisp);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc])/static_cast<Scalar>(weight); } 
      break;
    case 2: 
      loc  = c_dsa->locate(node, DofSet::Zdisp);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc])/static_cast<Scalar>(weight); } 
      break;
    case 3: 
      loc  = c_dsa->locate(node, DofSet::Xrot);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc])/static_cast<Scalar>(weight); } 
      break;
    case 4: 
      loc  = c_dsa->locate(node, DofSet::Yrot);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc])/static_cast<Scalar>(weight); } 
      break;
    case 5: 
      loc  = c_dsa->locate(node, DofSet::Zrot);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc])/static_cast<Scalar>(weight); } 
      break;
    case 6:
      locx  = c_dsa->locate( node, DofSet::Xdisp);
      loc1x =   dsa->locate( node, DofSet::Xdisp);
      if      ( locx >= 0  ) { dummyx = sol[locx]; ddx = grad[locx]; } 
      else if ( loc1x >= 0 ) { dummyx = this->Bcx(loc1x); ddx = 0; }
      else                   { dummyx = 0.0; ddx = 0; }
      nd = ScalarTypes::norm(dummyx);
      sd += nd*nd;
      locy  = c_dsa->locate( node, DofSet::Ydisp);
      loc1y =   dsa->locate( node, DofSet::Ydisp);
      if      ( locy >= 0  ) { dummyy = sol[locy]; ddy = grad[locy]; } 
      else if ( loc1y >= 0 ) { dummyy = this->Bcx(loc1y); ddy = 0; }
      else                   { dummyy = 0.0; ddy = 0;}
      nd = ScalarTypes::norm(dummyy);
      sd += nd*nd;
      locz  = c_dsa->locate( node, DofSet::Zdisp);
      loc1z =   dsa->locate( node, DofSet::Zdisp);
      if      ( locz >= 0  ) { dummyz = sol[locz]; ddz = grad[locz]; } 
      else if ( loc1z>= 0 ) { dummyz = this->Bcx(loc1z); ddz = 0; }
      else                  { dummyz = 0.0; ddz = 0; }
      nd = ScalarTypes::norm(dummyz);
      sd += nd*nd;
      val = sd>0? ScalarTypes::Real(dummyx*ScalarTypes::conj(ddx) + dummyy*ScalarTypes::conj(ddy) + dummyz*ScalarTypes::conj(ddz))/
	sqrt(sd)/static_cast<Scalar>(weight) : 0;
      break;
    default:
      printf("OptcritDisp::evaluate - Wrong Disptype !");
      assert(0);
    }
  return val;
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenSubDomain_opt<Scalar>::getAverageDisp(const Scalar* sol,
						const int* nodes, const int* dtypes,
						int size, Connectivity& nodeToSub)
{
  Scalar result = 0;
  for(int i=0; i<size; ++i)
    {
      result += getNodalDisp(sol, nodes[i], dtypes[i], nodeToSub.num(nodes[i]), false);
    }
  return result/static_cast<Scalar>(size);
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenSubDomain_opt<Scalar>::getGradAverageDisp(const Scalar* sol, const Scalar* grad,
						    const int* nodes, const int* dtypes,
						    int size, Connectivity& nodeToSub)
{
  Scalar result = 0;
  for(int i=0; i<size; ++i)
    {
      result += getGradNodalDisp(sol, grad, nodes[i], dtypes[i], nodeToSub.num(nodes[i]), false);
    }
  return result/static_cast<Scalar>(size);
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenSubDomain_opt<Scalar>::getGradduAverageDisp(const Scalar* sol, Scalar* adj,
						    const int* nodes, const int* dtypes,
						    int size, Connectivity& nodeToSub)
{
  for(int i=0; i<size; ++i)
    {
      getGradduNodalDisp(sol, adj, nodes[i], dtypes[i], nodeToSub.num(nodes[i]), false);
    }
  for(int i=0; i<this->numUncon(); ++i) { adj[i] /= size; }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenSubDomain_opt<Scalar>::getDisplevel(const Scalar* sol,
					      int size, double powFac, 
					      const int* nodeid, const int* distyp, 
					      const double* refVal, const Scalar* difVal, 
					      Connectivity& nodeToSub) 
{
  double sumVal = 0.0;  
  for(int i=0; i<size; ++i)
    {
      int node   = nodeid[i];
      int typ    = distyp[i];
      if(this->globalToLocal(node)<0) { continue; }
      Scalar disVal = (getNodalDisp(sol, node, typ, 1, false) - difVal[i])/refVal[i];
      sumVal += pow(ScalarTypes::doublify(disVal), powFac)/nodeToSub.num(node);
    }
  return sumVal;
}			       

//------------------------------------------------------------------------------
template<class Scalar>
double GenSubDomain_opt<Scalar>::getGradDisplevel(const Scalar* sol, const Scalar* grad,
						  int size, double powFac, 
						  const int* nodeid, const int* distyp, 
						  const double* refVal, 
						  const Scalar* difVal, const Scalar* gradDifVal,
						  Connectivity& nodeToSub) 
{
  double sumVal = 0.0;  
  for(int i=0; i<size; ++i)
    {
      int node    = nodeid[i];
      int typ     = distyp[i];
      int locNode = this->globalToLocal(node);
      if(locNode < 0) { continue; }
      Scalar disVal  = (getNodalDisp(sol, node, typ, 1, false) - difVal[i])/refVal[i];
      Scalar dDisVal = (getGradNodalDisp(sol, grad, node, typ, 1, false) - gradDifVal[i])/refVal[i];
      sumVal += powFac*pow(ScalarTypes::doublify(disVal), powFac-1)/nodeToSub.num(node) * 
	ScalarTypes::d_doublify(disVal, dDisVal);
    }
  return sumVal;
}
