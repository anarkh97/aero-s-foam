#include <Structopt.d/Element_opt.d/Element_opt.h>

//------------------------------------------------------------------------------
template<class Scalar>
Scalar Domain_opt::getNodalDisp(GenVector<Scalar>& sol, const Scalar* bcx, 
				int node, int dispTyp, bool dblfy)
{
  int loc, loc1;
  Scalar val = 0, dummy;
  double nd, sd=0.0;

  switch(dispTyp) 
    {    
    case 0: 
      loc  = c_dsa->locate(node, DofSet::Xdisp);
      loc1 =   dsa->locate(node, DofSet::Xdisp);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }      
      break;
    case 1: 
      loc  = c_dsa->locate(node, DofSet::Ydisp);
      loc1 =   dsa->locate(node, DofSet::Ydisp);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }
      break;
    case 2: 
      loc  = c_dsa->locate(node, DofSet::Zdisp);
      loc1 =   dsa->locate(node, DofSet::Zdisp);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }
      break;
    case 3: 
      loc  = c_dsa->locate( node, DofSet::Xrot);
      loc1 =   dsa->locate( node, DofSet::Xrot);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }
      break;
    case 4: 
      loc  = c_dsa->locate( node, DofSet::Yrot);
      loc1 =   dsa->locate( node, DofSet::Yrot);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }
      break;
    case 5: 
      loc  = c_dsa->locate( node, DofSet::Zrot);
      loc1 =   dsa->locate( node, DofSet::Zrot);
      if      ( loc >= 0  ) { val = sol[loc];  } 
      else if ( loc1 >= 0 ) { val = bcx[loc1]; }
      else                  { val = 0.0; }
      break;
    case 6:
      loc  = c_dsa->locate( node, DofSet::Xdisp);
      loc1 =   dsa->locate( node, DofSet::Xdisp);
      if      ( loc >= 0  ) { dummy = sol[loc];  } 
      else if ( loc1 >= 0 ) { dummy = bcx[loc1]; }
      else                  { dummy = 0.0; }
      nd = ScalarTypes::norm(dummy);
      sd += nd*nd;
      loc  = c_dsa->locate( node, DofSet::Ydisp);
      loc1 =   dsa->locate( node, DofSet::Ydisp);
      if      ( loc >= 0  ) { dummy = sol[loc];  } 
      else if ( loc1 >= 0 ) { dummy = bcx[loc1]; }
      else                  { dummy = 0.0; }
      nd = ScalarTypes::norm(dummy);
      sd += nd*nd;
      loc  = c_dsa->locate( node, DofSet::Zdisp);
      loc1 =   dsa->locate( node, DofSet::Zdisp);
      if      ( loc >= 0  ) { dummy = sol[loc];  } 
      else if ( loc1 >= 0 ) { dummy = bcx[loc1]; }
      else                  { dummy = 0.0; }
      nd = ScalarTypes::norm(dummy);
      sd += nd*nd;
      val = sqrt(sd);
      break;
    default:
      printf("OptcritDisp::evaluate - Wrong Disptype !");
      assert(0);
    }
  return dblfy? ScalarTypes::doublify(val) : val;
}


//------------------------------------------------------------------------------
template<class Scalar>
Scalar Domain_opt::getGradNodalDisp(GenVector<Scalar>& sol, GenVector<Scalar>& grad,
				    const Scalar* bcx,
				    int node, int dispTyp, bool dblfy)
{
  int loc, locx, locy, locz, loc1x, loc1y, loc1z;
  Scalar val = 0;
  Scalar dummyx, dummyy, dummyz, ddx, ddy, ddz;
  double nd, sd=0.0;

  switch(dispTyp) 
    {    
    case 0: 
      loc  = c_dsa->locate(node, DofSet::Xdisp);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc]); } 
      break;
    case 1: 
      loc  = c_dsa->locate(node, DofSet::Ydisp);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc]); } 
      break;
    case 2: 
      loc  = c_dsa->locate(node, DofSet::Zdisp);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc]); } 
      break;
    case 3: 
      loc  = c_dsa->locate(node, DofSet::Xrot);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc]); } 
      break;
    case 4: 
      loc  = c_dsa->locate(node, DofSet::Yrot);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc]); } 
      break;
    case 5: 
      loc  = c_dsa->locate(node, DofSet::Zrot);
      if      ( loc >= 0  ) { val = (dblfy? ScalarTypes::d_doublify(sol[loc], grad[loc]) :
				     grad[loc]); } 
      break;
    case 6:
      locx  = c_dsa->locate( node, DofSet::Xdisp);
      loc1x =   dsa->locate( node, DofSet::Xdisp);
      if      ( locx >= 0  ) { dummyx = sol[locx]; ddx = grad[locx]; } 
      else if ( loc1x >= 0 ) { dummyx = bcx[loc1x]; ddx = 0; }
      else                   { dummyx = 0.0; ddx = 0; }
      nd = ScalarTypes::norm(dummyx);
      sd += nd*nd;
      locy  = c_dsa->locate( node, DofSet::Ydisp);
      loc1y =   dsa->locate( node, DofSet::Ydisp);
      if      ( locy >= 0  ) { dummyy = sol[locy]; ddy = grad[locy]; } 
      else if ( loc1y >= 0 ) { dummyy = bcx[loc1y]; ddy = 0; }
      else                   { dummyy = 0.0; ddy = 0;}
      nd = ScalarTypes::norm(dummyy);
      sd += nd*nd;
      locz  = c_dsa->locate( node, DofSet::Zdisp);
      loc1z =   dsa->locate( node, DofSet::Zdisp);
      if      ( locz >= 0  ) { dummyz = sol[locz]; ddz = grad[locz]; } 
      else if ( loc1z>= 0 ) { dummyz = bcx[loc1z]; ddz = 0; }
      else                  { dummyz = 0.0; ddz = 0; }
      nd = ScalarTypes::norm(dummyz);
      sd += nd*nd;
      val = sd>0? ScalarTypes::Real(dummyx*ScalarTypes::conj(ddx) + 
				    dummyy*ScalarTypes::conj(ddy) + 
				    dummyz*ScalarTypes::conj(ddz))/
	sqrt(sd) : 0;
      break;
    default:
      printf("OptcritDisp::evaluate - Wrong Disptype !");
      assert(0);
    }
  return val;
}

//------------------------------------------------------------------------------
template<class Scalar>
void Domain_opt::getGradduNodalDisp(GenVector<Scalar>& sol, GenVector<Scalar>& adj,
				    const Scalar* bcx,
				    int node, int dispTyp,
				    bool dblfy)
{
  int loc, locx, locy, locz, loc1x, loc1y, loc1z;  
  Scalar dummyx, dummyy, dummyz, val = 0;
  double nd, sd=0.0;
  
  switch(dispTyp) 
    {    
    case 0: 
      loc  = c_dsa->locate(node, DofSet::Xdisp);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>()); }
      break;
    case 1: 
      loc  = c_dsa->locate(node, DofSet::Ydisp);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>()); }
      break;
    case 2: 
      loc  = c_dsa->locate(node, DofSet::Zdisp);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>()); }
      break;
    case 3: 
      loc  = c_dsa->locate(node, DofSet::Xrot);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>()); }
      break;
    case 4: 
      loc  = c_dsa->locate(node, DofSet::Yrot);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>()); }
      break;
    case 5: 
      loc  = c_dsa->locate(node, DofSet::Zrot);
      if      ( loc >= 0  ) { adj[loc] += ScalarTypes::conj(dblfy? ScalarTypes::unify(sol[loc]) :
							    ScalarTypes::id<Scalar>()); }
      break;
    case 6:
      locx  = c_dsa->locate( node, DofSet::Xdisp);
      loc1x =   dsa->locate( node, DofSet::Xdisp);
      if      ( locx >= 0  ) { dummyx = sol[locx];  } 
      else if ( loc1x >= 0 ) { dummyx = bcx[loc1x]; }
      else                   { dummyx = 0.0; }
      nd = ScalarTypes::norm(dummyx);
      sd += nd*nd;
      locy  = c_dsa->locate( node, DofSet::Ydisp);
      loc1y =   dsa->locate( node, DofSet::Ydisp);
      if      ( locy >= 0  ) { dummyy = sol[locy];  } 
      else if ( loc1y >= 0 ) { dummyy = bcx[loc1y]; }
      else                   { dummyy = 0.0; }
      nd = ScalarTypes::norm(dummyy);
      sd += nd*nd;
      locz  = c_dsa->locate( node, DofSet::Zdisp);
      loc1z =   dsa->locate( node, DofSet::Zdisp);
      if      ( locz >= 0  ) { dummyz = sol[locz];  } 
      else if ( loc1z>= 0 ) { dummyz = bcx[loc1z]; }
      else                  { dummyz = 0.0; }
      nd = ScalarTypes::norm(dummyz);
      sd += nd*nd;
      val = sqrt(sd);
      if(sd > 0)
	{
	  if(locx>=0) { adj[locx] += ScalarTypes::conj(dummyx)/val; }
	  if(locy>=0) { adj[locy] += ScalarTypes::conj(dummyy)/val; }
	  if(locz>=0) { adj[locz] += ScalarTypes::conj(dummyz)/val; }
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
Scalar Domain_opt::getAverageDisp(GenVector<Scalar>& sol,
				  const Scalar* bcx,
				  const int* nodes, const int* dtypes,
				  int size)
{
  Scalar result = 0;
  for(int i=0; i<size; ++i)
    {
      result += getNodalDisp(sol, bcx, nodes[i], dtypes[i], false);
    }
  return result/static_cast<Scalar>(size);
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar Domain_opt::getGradAverageDisp(GenVector<Scalar>& sol, GenVector<Scalar>& grad,
				      const Scalar* bcx,
				      const int* nodes, const int* dtypes,
				      int size)
{
  Scalar result = 0;
  for(int i=0; i<size; ++i)
    {
      result += getGradNodalDisp(sol, grad, bcx, nodes[i], dtypes[i], false);
    }
  return result/static_cast<Scalar>(size);
}

//------------------------------------------------------------------------------
template<class Scalar>
void Domain_opt::getGradduAverageDisp(GenVector<Scalar>& sol, GenVector<Scalar>& adj,
				      const Scalar* bcx,
				      const int* nodes, const int* dtypes,
				      int size)
{
  for(int i=0; i<size; ++i)
    {
      getGradduNodalDisp(sol, adj, bcx, nodes[i], dtypes[i], false);
    }
  adj *= 1.0/size;
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void Domain_opt::buildDStiffnesDSmDU(GenVector<Scalar>& rhs, 
				     GenVector<Scalar>& sol, 
				     const Scalar* bcx, int flag)
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
      int locIele = (flag==0)? actInf->getEntry(actelm) : actelm;
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
	   { (*getElDisp<Scalar>())[k] = bcx[(*allDOFs)[locIele][k]]; }
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
Scalar Domain_opt::DStiffnesDSmDUAdj(GenVector<Scalar>& sol, GenVector<Scalar>& adj,
				     const Scalar* bcx)
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
      int locIele = actInf->getEntry(actelm);
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
	   { (*getElDisp<Scalar>())[k] = bcx[(*allDOFs)[locIele][k]]; }
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
void Domain_opt::buildDRHSforceDs(GenVector<Scalar>& pseudoforce)
{
  // compute derivative of external force
  for(int i=0; i < numNeuman; ++i) 
    {
      int dof  = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum));
      if(dof < 0) { continue; }
      pseudoforce[dof] += gradnbc[i].val;
    }

  // add gravity forces
  if (gravityFlag()) 
    { fprintf(stderr,"WATCH: Derivatives of Gravity Force neglected\n"); }

  // add pressure load
//   if (geoSource->pressureFlag()) 
//     buildDPressureForceDs(pseudoforce, gs);

  // add thermal load if from input file (i.e. constant)
  // or in the case of adjoint method
//   if (sinfo.thermalLoadFlag) {
//     buildDThermalForceDs(pseudoforce,gs);
//   }
//   if (sinfo.thermoeFlag >= 0 && flag == 2) {
//     buildDThermalForceDs(pseudoforce,gs);
//   }
}

//------------------------------------------------------------------------------
template<class Scalar>
double Domain_opt::getDisplevel(GenVector<Scalar>& sol,	const Scalar* bcx,
				int aFlag, int quoFlag, int size, 
				const int* nodeid, const int* distyp, 
				const double* refVal, const double* difVal, double powFac)
{
  Scalar aveDisp;
  boost::scoped_array<Scalar> difVal1(new Scalar[size]);
  switch(aFlag)
    {
    case 0:
      std::copy(difVal, difVal+size, difVal1.get());
      break;
    case 1:
      aveDisp = getAverageDisp(sol, bcx, nodeid, distyp, size);
      std::fill(difVal1.get(), difVal1.get() + size, aveDisp);
      break;
    case 2:
      for(int i=0; i<size; ++i) 
	{ difVal1[i] = getNodalDisp(sol, bcx, static_cast<int>(difVal[i])-1, distyp[i], false); }
      break;
    default:
      assert(0);      
    }

  double sumVal = 0.0;  
  for(int i=0; i<size; ++i)
    {
      int node   = nodeid[i];
      int typ    = distyp[i];
      Scalar disVal = (getNodalDisp(sol, bcx, node, typ, false) - difVal1[i])/refVal[i];
      sumVal += pow(ScalarTypes::doublify(disVal), powFac);
    }
  double quoVal = quoFlag ? size : 1.0;
  return pow(sumVal/quoVal, 1.0/powFac);  
}

//------------------------------------------------------------------------------
template<class Scalar> 
double Domain_opt::getGradDisplevel(GenVector<Scalar>& sol, GenVector<Scalar>& grad, 
				    const Scalar* bcx,
				    int aFlag, int quoFlag, int size, 
				    const int* nodeid, const int* distyp, 
				    const double* refVal, const double* difVal, double powFac)
{
    Scalar aveDisp, gradAveDisp;
  boost::scoped_array<Scalar> difVal1(new Scalar[size]);
  boost::scoped_array<Scalar> gradDifVal1(new Scalar[size]);

  switch(aFlag)
    {
    case 0:
      std::copy(difVal, difVal+size, difVal1.get());
      std::fill(gradDifVal1.get(), gradDifVal1.get() + size, 0);
      break;
    case 1:
      aveDisp = getAverageDisp(sol, bcx, nodeid, distyp, size);
      std::fill(difVal1.get(), difVal1.get() + size, aveDisp);
      gradAveDisp = getGradAverageDisp(sol, grad, bcx, nodeid, distyp, size);
      std::fill(gradDifVal1.get(), gradDifVal1.get() + size, gradAveDisp);
      break;
    case 2:
      for(int i=0; i<size; ++i) 
	{ 
	  difVal1[i]     = getNodalDisp(sol, bcx, static_cast<int>(difVal[i])-1, distyp[i], false);
	  gradDifVal1[i] = getGradNodalDisp(sol, grad, bcx, static_cast<int>(difVal[i])-1, distyp[i], false);
	}
      break;
    default:
      assert(0);      
    }

  double sumVal = 0.0;  
  for(int i=0; i<size; ++i)
    {
      int node    = nodeid[i];
      int typ     = distyp[i];
      Scalar disVal  = (getNodalDisp(sol, bcx, node, typ, false) - difVal1[i])/refVal[i];
      Scalar dDisVal = (getGradNodalDisp(sol, grad, bcx, node, typ, false) - gradDifVal1[i])/refVal[i];
      sumVal += powFac*pow(ScalarTypes::doublify(disVal), powFac-1) * 
	ScalarTypes::d_doublify(disVal, dDisVal);
    }

  double quoVal = quoFlag ? size : 1.0;
  double dl = getDisplevel(sol, bcx, aFlag, quoFlag, size, nodeid, distyp, refVal, difVal, powFac);
  return dl == 0?
    0 :
    pow(dl, 1.0 - powFac)/powFac*sumVal/quoVal;
}

//------------------------------------------------------------------------------
template<class Scalar> 
void Domain_opt::getGradduDisplevel(GenVector<Scalar>& sol, GenVector<Scalar>& adj,
				    const Scalar* bcx,
				    int aFlag, int quoFlag, int size, 
				    const int* nodeid, const int* distyp, 
				    const double* refVal, const double* difVal, double powFac)
{
  Scalar aveDisp;
  GenVector<Scalar> difAdj1(sol.size());
  GenVector<Scalar> difAdj2(sol.size());

  switch(aFlag)
    {
    case 0:
      for(int i=0; i<size; ++i) 
	{	  
	  Scalar disVal  = (getNodalDisp(sol, bcx, nodeid[i], distyp[i], false) - difVal[i])/refVal[i];
	  difAdj1.zero();
	  getGradduNodalDisp(sol, difAdj1, bcx, nodeid[i], distyp[i], false);
   	  scalAdjVec(difAdj1.size(), difAdj1.data(), 
   		     powFac*pow(ScalarTypes::doublify(disVal), powFac-1)/refVal[i]*
   		     unify(disVal, distyp[i]));
	  adj += difAdj1;
	}
      break;
    case 1:
      assert(typeid(Scalar)!=typeid(DComplex)); // not implemented for complex
      aveDisp = getAverageDisp(sol, bcx, nodeid, distyp, size);
      difAdj2.zero();
      getGradduAverageDisp(sol, difAdj2, bcx, nodeid, distyp, size);
      for(int i=0; i<size; ++i) 
	{
	  Scalar disVal  = (getNodalDisp(sol, bcx, nodeid[i], distyp[i], false) - aveDisp)/refVal[i];
	  difAdj1.zero();
	  getGradduNodalDisp(sol, difAdj1, bcx, nodeid[i], distyp[i], false);
	  difAdj1 -= difAdj2;
	  scalAdjVec(difAdj1.size(), difAdj1.data(), 
		     powFac*pow(ScalarTypes::doublify(disVal), powFac-1)/refVal[i]*
		     ScalarTypes::unify(disVal));
	  adj += difAdj1;
	}
      break;
    case 2:
      for(int i=0; i<size; ++i) 
	{ 
	  Scalar disVal  = (getNodalDisp(sol, bcx, nodeid[i], distyp[i], false) - 
			    getNodalDisp(sol, bcx, static_cast<int>(difVal[i])-1, distyp[i], false))/refVal[i];
	  difAdj1.zero();
	  getGradduNodalDisp(sol, difAdj1, bcx, nodeid[i], distyp[i], false);
	  difAdj2.zero();
	  int refNode = static_cast<int>(difVal[i])-1;
	  assert(refNode >= 0);
	  getGradduNodalDisp(sol, difAdj2, bcx, refNode, distyp[i], false);
	  difAdj1 -= difAdj2;
	  scalAdjVec(difAdj1.size(), difAdj1.data(), 
		     powFac*pow(ScalarTypes::doublify(disVal), powFac-1)/refVal[i]*
		     unify(disVal, distyp[i]));
	  adj += difAdj1;
	}
      break;
    default:
      assert(0);
      break;
    }
  double quoVal = quoFlag ? size : 1.0;
  double dl = getDisplevel(sol, bcx, aFlag, quoFlag, size, nodeid, distyp, refVal, difVal, powFac);
  if(dl == 0)
    { adj.zero(); }
  else
    { adj *= pow(dl, 1.0 - powFac)/powFac/quoVal; }
  return;  
}

//------------------------------------------------------------------------------
template<>
inline GenVector<double>* Domain_opt::getElDisp<double>() { return this->elDisp; }
template<>
inline GenVector<DComplex>* Domain_opt::getElDisp<DComplex>() { return this->c_elDisp.get(); }
template<>
inline GenVector<double>* Domain_opt::getElKelD<double>() { return this->elKelD.get(); }
template<>
inline GenVector<DComplex>* Domain_opt::getElKelD<DComplex>() { return this->c_elKelD.get(); }
template<>
inline GenVector<double>* Domain_opt::getElGrad<double>() { return this->elGrad.get(); }
template<>
inline GenVector<DComplex>* Domain_opt::getElGrad<DComplex>() { return this->c_elGrad.get(); }
template<>
inline GenVector<double>* Domain_opt::getElAdj<double>() { return this->elAdj.get(); }
template<>
inline GenVector<DComplex>* Domain_opt::getElAdj<DComplex>() { return this->c_elAdj.get(); }
template<>
inline GenVector<double>* Domain_opt::getElRHS<double>() { return this->elRHS.get(); }
template<>
inline GenVector<DComplex>* Domain_opt::getElRHS<DComplex>() { return this->c_elRHS.get(); }


//------------------------------------------------------------------------------
template<class Scalar> 
void Domain_opt::densProjectVector(const DensityProjData& dpd, Scalar* vec)
{
  /* do nothing  */
  return;
  if (dpd.densProjFlag > 0)
    {
      for(int i=0; i<this->numUncon(); ++i)
	{ vec[i] *= densProj[i]; }
    }  
  return;
}

//------------------------------------------------------------------------------
template<class Scalar> 
void Domain_opt::densProjectVectorInv(const DensityProjData& dpd, Scalar* vec)
{
  /* do nothing */
  return;
  if (dpd.densProjFlag > 0)
    {
      for(int i=0; i<this->numUncon(); ++i)
	{ vec[i] /= densProj[i]; }
    }  
  return;
}
