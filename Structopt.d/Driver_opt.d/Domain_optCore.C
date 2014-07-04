#ifdef STRUCTOPT

#include <Structopt.d/Structopt_sd.h>
#include <Structopt.d/Driver_opt.d/Domain_opt.h>
#include <Structopt.d/Element_opt.d/Element_opt.h>
#include <Math.d/FullSquareMatrix.h>

#include <map>
#include <algorithm>
#include <utility>
#include <cfloat>

Domain_opt::Domain_opt() : 
  Domain(), analysisList(std::auto_ptr<SolverInfo>(0)), structoptFlag(0), reliabilityFlag(0), numAnalysis(0),
  numStaticAna(0), optgradFlag(0), actLC(0), activeAnalysis(0), actInf(0), semiSAflag(0),
  gradnbc(0), gradcnbc(0)
{}

Domain_opt::~Domain_opt() {}


//----------------------------------------------------------------------------
void 
Domain_opt::activateAnalysis (int sact, int nbFlag) 
{ 
  // check if analysis exists

  if (sact >= numAnalysis ) {
    fprintf(stderr," *** ERROR: Analysis Id %d > number of Analysis %d\n",
                    sact+1,numAnalysis);
    exit(-1);
  } 

  // update analysis id

  sinfo = *(analysisList[sact]); 
  
  activeAnalysis = sact; 
  
  // set Neuman boundary conditions to first load case by default

  if (nbFlag) setActiveLC(0);
} 

//----------------------------------------------------------------------------
int Domain_opt::getNumLC()
{ 
  return numAnaNmLC[activeAnalysis];
}

//----------------------------------------------------------------------------

int Domain_opt::nNeumann(int anaID, int anaLC)
{ 
 if (anaID >= numAnalysis) {
   fprintf(stderr," *** ERROR: analysis id larger > number of analysis\n"); 
   exit(-1);
 }  

 if ( anaLC >= numAnaNmLC[anaID] ) {
   fprintf(stderr," *** ERROR: load case id %d > number of load cases of analysis %d\n",
                    anaLC+1,anaID+1);
   exit(-1);
 }

 int lc = mapAnaNmLC[anaID][anaLC];

 if (lc >= numLC ) {
   fprintf(stderr," *** ERROR: active loadcase to be set > number of global LCs\n");
   exit(-1);
 }
   
 return numNmLC[lc]; 
}


//----------------------------------------------------------------------------

BCond* Domain_opt::getNBC(int anaID, int anaLC)
{ 
 if (anaID >= numAnalysis) {
   fprintf(stderr," *** ERROR: analysis id larger > number of analysis\n"); 
   exit(-1);
 }  

 if ( anaLC >= numAnaNmLC[anaID] ) {
   fprintf(stderr," *** ERROR: load case id %d > number of load cases of analysis %d\n",
                    anaLC+1,anaID+1);
   exit(-1);
 }

 int lc = mapAnaNmLC[anaID][anaLC];

 if (lc >= numLC ) {
   fprintf(stderr," *** ERROR: active loadcase to be set > number of global LCs\n");
   exit(-1);
 }

 return nbcLC[lc]; 
}



//----------------------------------------------------------------------------
void 
Domain_opt::setActiveLC(int anaLC, int *bc, double *bcx)
{ 

 // set free so far loaded DOFS


 int i;
 if (bc != 0 && bcx != 0) {
   for(i=0; i<numNeuman; ++i) {
     int dof  = dsa->locate(nbc[i].nnum, 1 << nbc[i].dofnum);
     int numdof = dsa->size();
     if(dof >= numdof) {
       fprintf(stderr," *** ERROR: check input, found error for FORCE"
        	      " (node %d, dof %d)\n",nbc[i].nnum+1,nbc[i].dofnum+1);
       exit(-1);
     }
     if(dof < 0) {
       fprintf(stderr," *** WARNING: dof does not exist: node %d dof %d\n",
     	    	        nbc[i].nnum+1,nbc[i].dofnum+1);
       continue;
     }
     bc[dof]  = BCFREE;
     bcx[dof] = 0.0;
   }
 }
 
 // check if local load case exists
 
 if ( anaLC >= numAnaNmLC[activeAnalysis] ) {
   fprintf(stderr," *** ERROR: load case id %d > number of load cases of analysis %d\n",
                    anaLC+1,activeAnalysis+1);
   fprintf(stderr,"             number of load cases of this analysis: %d\n",numAnaNmLC[activeAnalysis]);		    
   exit(-1);
 }
 
 int lc = mapAnaNmLC[activeAnalysis][anaLC];
 
 // check if global load case exists

 if (lc >= numLC ) {
   fprintf(stderr," *** ERROR: active loadcase to be set > number of global LCs\n");
   exit(-1);
 }

 // update for new load case

 actLC     = lc; 
 numNeuman = numNmLC[lc];
 nbc       = nbcLC[lc];

 // set new load case
 if (bc != 0 && bcx != 0) {
   for(i=0; i<numNeuman; ++i) {
      int dof  = dsa->locate(nbc[i].nnum, 1 << nbc[i].dofnum);
      int numdof = dsa->size();
      if(dof >= numdof) {
        fprintf(stderr," *** ERROR: check input, found error for FORCE"
         	      " (node %d, dof %d)\n",nbc[i].nnum+1,nbc[i].dofnum+1);
        exit(-1);
      }
      if(dof < 0) {
        fprintf(stderr," *** WARNING: dof does not exist: node %d dof %d\n",
     	    	        nbc[i].nnum+1,nbc[i].dofnum+1);
        continue;
      }
      if(bc[dof] == BCLOAD) {
        fprintf(stderr," *** WARNING: check input, repeated FORCE node %d"
     		      ", dof %d for load case %d\n",
                       nbc[i].nnum+1,nbc[i].dofnum+1,lc+1);
      }
      bc[dof]  = BCLOAD;
      bcx[dof] = nbc[i].val;
    }
  }

#ifdef STRUCTOPT
  // set gradient of Neuman boundary conditions
  setGradActiveLC(actLC);  
#endif
}

double
Domain_opt::getNodalStressStrain(Vector &sol, double *bcx,
				 int node, int stressIndex, int surface, 
				 int eleFlag)
{
  // Allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  // check for elemental evaluation; Number of Connected Elements
  int nconele;
  if (eleFlag >= 0) 
    nconele = 1;
  else
    nconele = nodeToElem->num(node);

  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

  // Watch out for this, as the maximum number of nodes may increase!
  int maxNodesPerElement = 20;

  double stress=0.0;
  double weight=0.0;

  int icon,k,iele;
  for(icon=0; icon<nconele; ++icon) {
 
     elDisp->zero();
     elstress->zero();
     elweight->zero();
     //elNodTemp->zero();

     if (eleFlag >= 0) 
       iele=eleFlag;
     else
       iele = (*nodeToElem)[node][icon];

     int NodesPerElement = elemToNode->num(iele);

     packedEset[iele]->nodes(nodeNumbers);

     // ... DETERMINE ELEMENT DISPLACEMENT VECTOR

     for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elDisp)[k] = sol[cn];
        else
          (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
     }

     // ... DETERMINE NODAL TEMPERATURE VECTOR

     int iNode;
     /**
     if (nodalTemperatures) {
       for (iNode = 0; iNode < NodesPerElement; ++iNode) {
         if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
           (*elNodTemp)[iNode] = packedEset[iele]->prop->Ta;
         else
           (*elNodTemp)[iNode] = nodalTemperatures[nodeNumbers[iNode]];
       }
     }
     **/

     // ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT

     packedEset[iele]->getVonMises(*elstress, *elweight, nodes,
                                   *elDisp, stressIndex, surface);
     //elNodTemp->data());

     // ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
     // ... eleFlag <  0: average stress in node
     // ... eleFlag >= 0: average stress of particular element

     if (eleFlag < 0) {
       for(k=0; k<NodesPerElement; ++k) { 
         int actnod = (*elemToNode)[iele][k];
         if (actnod == node) { 
           stress += (*elstress)[k]*(*elweight)[k];
           weight += (*elweight)[k];
         }
       }
     }
     else {
       for(k=0; k<NodesPerElement; ++k) { 
         stress += (*elstress)[k]*(*elweight)[k];
         weight += (*elweight)[k];
       }
     }
  }

  delete [] nodeNumbers;
  if(sinfo.thermalLoadFlag) delete [] nodalTemperatures;

  return weight != 0 ? (stress/weight) : 0.0;  
}


double
Domain_opt::getStressIntegral(Vector &sol, double *bcx, double& refStress, 
                          double& powFac, int* eleList, int listSize, 
                          int stressIndex, int areaFlag)
{ 
  // allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  //if(elDisp == 0)    elDisp    = new Vector(maxNumDOFs,0.0);
  //if(elNodTemp == 0) elNodTemp = new Vector(maxNumNodes,0.0);

  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

  int    iele, k;
  double vmint, vol;

  //int    elemax;
  //double valmax;

  // WATCH: valvmint and valvol are global variables and are used for sensitivity analysis
  valvmint  = 0.0;
  valvol    = 0.0;
  //valmax    = 0.0;

  for(iele=0; iele<numele; ++iele) {

     if ( ! checkList(iele,eleList,listSize) ) continue;

     int NodesPerElement = elemToNode->num(iele);

     packedEset[iele]->nodes(nodeNumbers);

     // ... DETERMINE ELEMENT DISPLACEMENT VECTOR

     for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elDisp)[k] =  sol[cn];
        else
          (*elDisp)[k] =  bcx[(*allDOFs)[iele][k]];
     }

     /*
     int iNode;
     if (nodalTemperatures) {
       for (iNode = 0; iNode < NodesPerElement; ++iNode) {
         if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
           (*elNodTemp)[iNode] = packedEset[iele]->prop->Ta;
         else
           (*elNodTemp)[iNode] = nodalTemperatures[nodeNumbers[iNode]];
       }
     }
     */

     // ... calculate integral (sig_vm/sigbar)^fac over ele-volume

     vmint = vol = 0.0;

     packedEset[iele]->getVonMisesInt(nodes, *elDisp, refStress, powFac, areaFlag, 
                                      vmint, vol); //, elNodTemp->data());

     //if ( abs(vmint) > valmax) {
     //  valmax = abs(vmint);
     //  elemax = iele;
     //}

     valvmint += vmint;
     valvol   += vol;
  }

  delete [] nodeNumbers;
  if(sinfo.thermalLoadFlag) delete [] nodalTemperatures;

  if (!areaFlag) valvol=1.0;

  //fprintf(stderr,"strint: max. element stress for typ %d : # %d  =  %e\n",
  //        stressIndex,elemax,valmax); 

  return (valvmint != 0) ? pow(valvmint/valvol,1.0/powFac) : 0;   
}

double
Domain_opt::getDCmass(Vector &sol, double *bcx, double& powFac, int* eleList, int listSize)
{
  int iele, k;

  double totDCmass = 0.0;

  //if(elDisp == 0) elDisp = new Vector(maxNumDOFs,0.0);

  for(iele=0; iele<numele; ++iele) {

     if ( ! checkList(iele,eleList,listSize) ) continue;

     // ... DETERMINE ELEMENT DISPLACEMENT VECTOR

     for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elDisp)[k] =  sol[cn];
        else
          (*elDisp)[k] =  bcx[(*allDOFs)[iele][k]];
     }

     double eleDCmass = packedEset[iele]->getDCmass(nodes, *elDisp, powFac);

     totDCmass += eleDCmass;
  }      

  return totDCmass;
}

/*****
void
Domain_opt::getElementForces(Vector &sol, double *bcx, int fileNumber, 
                         int forceIndex, double time)
{
  // allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  // ...REMARK:  FORCES ARE CALCULATED FOR BARS AND BEAMS CURRENTLY 

  if(elDisp == 0)    elDisp    = new Vector(maxNumDOFs,0.0);
  if(elNodTemp == 0) elNodTemp = new Vector(maxNumNodes,0.0);

  // ... Watch out for this as the max number of nodes may increase
  int NodesPerElement = 2;

  FullM forces(numele, NodesPerElement);

  int i, iele;
  // note, we are reusing elstress here for an elemental force vector
  if(elstress == 0) {
    int NodesPerElement, maxNodesPerElement=0;
    for(iele=0; iele<numele; ++iele) {
      NodesPerElement = elemToNode->num(iele);
      maxNodesPerElement = std::max(maxNodesPerElement, NodesPerElement);
    }
    elstress = new Vector(maxNodesPerElement, 0.0);
  }

  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

  for(iele=0; iele<numele; ++iele) {
    packedEset[iele]->nodes(nodeNumbers);

    elDisp->zero();
    elNodTemp->zero();

// ... DETERMINE ELEMENT DISPLACEMENT VECTOR

    int k;
    for(k=0; k<allDOFs->num(iele); ++k) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
      if(cn >= 0)
        (*elDisp)[k] = sol[cn];
      else
        (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
    }

    // ... CALCULATE INTERNAL FORCE VALUE FOR EACH ELEMENT

    // ... taking into account the new temperature in case of thermal coupling
    if(nodalTemperatures) {
      int iNode;
      for(iNode=0; iNode<NodesPerElement; ++iNode) {
        if(nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
     	  (*elNodTemp)[iNode] = packedEset[iele]->prop->Ta;
        else
    	  (*elNodTemp)[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }
    }

    packedEset[iele]->getIntrnForce(*elstress,nodes,elDisp->data(),
                                    forceIndex,elNodTemp->data());

// ... COPY ELEMENT'S NODAL FORCES INTO A TOTAL FORCE MATRIX

    for(i=0; i<NodesPerElement; ++i)
      forces[iele][i] = (*elstress)[i];

  }

// ... PRINT THE ELEMENT FORCES TO A FILE
   geoSource->outputElemVectors(fileNumber, forces.data(), numele, time);

// ... DELETE LOCAL ARRAY
   delete [] nodeNumbers;
   if(sinfo.thermalLoadFlag) delete [] nodalTemperatures;
}

void
Domain_opt::getElasticStiff(FullSquareMatrix *kelArray) 
{
 int iele;
 for(iele=0; iele < numele; ++iele) {
    kelArray[iele] = packedEset[iele]->stiffness(nodes,kelArray[iele].data());
  }
}
**/

/**
void
Domain_opt::getElasticForces(Vector &dsp, double *bcx, Vector &ext_f, double eta,
                         FullSquareMatrix *kelArray)
{

  int size    = sizeof(double)*maxNumDOFs*maxNumDOFs;
  double *kel = (double *) dbg_alloca(size);

  int iele;
  for(iele=0; iele<numele; ++iele) {
       
     int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs      = allDOFs->num(iele);
     
     Vector elDisp   (numEleDOFs,0.0);
     Vector elForce  (numEleDOFs,0.0);
    
     FullSquareMatrix karray(numEleDOFs,kel);  
     
     int k;
     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          elDisp[k] = dsp[cn]; }
        else {
          elDisp[k] = bcx[(*allDOFs)[iele][k]]; }
     }

     if (kelArray) {
       karray=kelArray[iele];
     }
     else {
       karray=packedEset[iele]->stiffness(nodes,karray.data());
     }

     int i,j;
     for (i=0;i<numEleDOFs;i++) {
       for (j=0;j<numEleDOFs;j++) {
         elForce[i] += karray[i][j]*elDisp[j] ;
       }
     }

     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          ext_f[cn] += eta * elForce[k]; }
     }
  }
}
**/

double
Domain_opt::getStructureMass(int* eleList, int listSize)
{
  double totmas = 0.0; // total mass of structure

  for(int iele=0; iele < numele; ++iele) 
    {
      if ( ! checkList(iele,eleList,listSize) ) { continue; }
      double elementMass = packedEset[iele]->getMass(nodes);
      totmas += elementMass;
    }

  // ADD DISCRETE MASS
  typedef std::pair<int, double> CMtp;
  typedef std::map<int,  CMtp> CMmap;
  CMmap nodeCM;
  for(DMassData *current = firstDiMass; current != 0; current = current->next) 
    {
      CMmap::iterator i = nodeCM.lower_bound(current->node);
      if(i == nodeCM.end() || i->first != current->node)
	{
	  // insert a new one
	  nodeCM.insert(i, CMmap::value_type(current->node, CMtp(1, current->diMass)));
	}
      else
	{
	  CMtp& cm = i->second;
	  ++cm.first;
	  cm.second  += current->diMass;
	}
    }

  for(CMmap::const_iterator i = nodeCM.begin();
      i != nodeCM.end(); ++i)
    {
      int nc = std::max(i->second.first, 6);
      double dm = i->second.second;
      totmas += dm/nc;
  }

  return totmas;
}

double 
Domain_opt::getMomOfInertia(int* eleList, int listSize, int momFlag)
{
  double Ixx    = 0.0;
  double Iyy    = 0.0;
  double Izz    = 0.0;

  int*nodeNumbers= new int[maxNumNodes];

  int iele;
  for(iele=0; iele < numele; ++iele) {

    if ( ! checkList(iele,eleList,listSize) ) continue;

    double elementMass = packedEset[iele]->getMass(nodes);

    int numNodesPerElement = packedEset[iele]->numNodes();
    packedEset[iele]->nodes( nodeNumbers );

    double massPerNode = elementMass/numNodesPerElement;

    int i;
    for(i=0; i<numNodesPerElement; ++i) {
       Node &node = nodes.getNode( nodeNumbers[i] );

       double x = node.x;
       double y = node.y;
       double z = node.z;

       Ixx += massPerNode*sqrt(y*y + z*z);
       Iyy += massPerNode*sqrt(x*x + z*z);
       Izz += massPerNode*sqrt(x*x + y*y);
    }
  }

  delete [] nodeNumbers;

  // Add discrete masses

  DMassData *current = firstDiMass;

  while(current != 0) {

    int n = current->node;

    // check if the node is in the model

    if(nodes.exist(n)) { 

      Node &node = nodes.getNode(n);

      switch (current->dof) {
        case 0:
          Iyy += (current->diMass)*sqrt(node.z*node.z);
          Izz += (current->diMass)*sqrt(node.y*node.y);
          break;
        case 1:
          Ixx += (current->diMass)*sqrt(node.z*node.z);
          Izz += (current->diMass)*sqrt(node.x*node.x);
          break;
        case 2:
          Ixx += (current->diMass)*sqrt(node.y*node.y);
          Iyy += (current->diMass)*sqrt(node.x*node.x);
          break;
        case 3:
          Ixx += (current->diMass); 
          break;
        case 4:
          Iyy += (current->diMass); 
          break;
        case 5:
          Izz += (current->diMass); 
          break;
        default:
          fprintf(stderr,"ERROR: Dof for Dimass does not exists\n");
          exit(-1);
      }
    }
    current = current->next;
  }

  switch (momFlag) {
    case 3:
      return Ixx;
      break;
    case 4:
      return Iyy;
      break;
    case 5:
      return Izz;
      break;
  }

  return 0;
}

double Domain_opt::getVariationNodalPos(int nodeId, int dir)
{

  Node& orgnod = nodescopy->getNode(nodeId);
  
  double val;
  double dx,dy,dz;
  
  switch (dir) {
    
    case 0: 
      val = nodes[nodeId]->x - orgnod.x ;
      break;
    case 1:
      val = nodes[nodeId]->y - orgnod.y ;
      break;
    case 2:
      val = nodes[nodeId]->z - orgnod.z ;
      break;
    case 6:
      dx  = nodes[nodeId]->x - orgnod.x ;
      dy  = nodes[nodeId]->y - orgnod.y ;
      dz  = nodes[nodeId]->z - orgnod.z ;
      val = dx*dx + dy*dy + dz*dz;
      val = sqrt(val);
      break;
    default:
      printf("Variation of Nodal Position: evaluate - Wrong Type !");
      exit(-1);
  }
  
  return val;   
}

double
Domain_opt::getStrainEnergy(Vector &sol, double *bcx, SparseMatrix * gStiff,
                        FullSquareMatrix *kelArray, int* eleList, int listSize)
{
  double energy=0.0;

  // Evaluation of Strain Energy if global Stiffness Matrix is given
  
  if (gStiff && ! eleList) {
      
    Vector tmpVec(c_dsa->size());
    
    gStiff->mult(sol,tmpVec);

    energy = 0.5 * (sol * tmpVec);
  
    return energy;
  }
  
  // Evaluation of Strain Energy if local Stiffness Matices are given
  // or has to be determined

  int iele;
  for(iele=0; iele<numele; ++iele) {

     if ( ! checkList(iele,eleList,listSize) ) continue;
       
     int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs      = allDOFs->num(iele);
     
     elDisp->zero();
     elKelD->zero();
   
     FullSquareMatrix karray(numEleDOFs,kelData.get());
     
     int k;
     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elDisp)[k] = sol[cn];
        else
          (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
     }

     if (kelArray) {
       karray=kelArray[iele];
     }
     else {
       karray=packedEset[iele]->stiffness(nodes,karray.data());
     }

     karray.multiply(*elDisp,*elKelD, 1.0); 
				   
     energy += 0.5*(*elKelD * *elDisp);
  }

  return energy;
}

double
Domain_opt::getNodalInternalForce(Vector &dsp, double *bcx, int node, int typ)
{
  // Initialize nodal force
  double nodalForce = 0;

  // Number of connected elements
  int nconele = nodeToElem->num(node);

  if(!elDisp)  elDisp = new Vector(maxNumDOFs,0.0);
  if(!elForce.get()) elForce.reset(new Vector(maxNumDOFs,0.0));
  if (!kelData.get()) kelData.reset(new double[maxNumDOFs*maxNumDOFs]);

  // Extract the relevant DOF
  int dof = dsa->locate(node, (1 << typ));
  if (dof < 0) {
    fprintf(stderr,"Error: Force at node % can not be computed\n",node+1);
    exit(-1);
  }

  // Loop over all connected elements
  int icon;
  for(icon=0; icon<nconele; ++icon) {

     int iele = (*nodeToElem)[node][icon];
       
     int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs      = allDOFs->num(iele);
     
     elDisp->zero(); 
     elForce->zero();
    
     FullSquareMatrix karray(numEleDOFs,kelData.get());  
     
     int k;
     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          (*elDisp)[k] = dsp[cn]; }
        else {
          (*elDisp)[k] = bcx[(*allDOFs)[iele][k]]; }
     }

     karray=packedEset[iele]->stiffness(nodes,karray.data());

     int i,j;
     for (i=0;i<numEleDOFs;i++) {
       for (j=0;j<numEleDOFs;j++) {
         (*elForce)[i] += karray[i][j]*(*elDisp)[j] ;
       }
     }

     for(k=0; k<numEleDOFs; ++k) {
       if ( (*allDOFs)[iele][k] == dof ) nodalForce += (*elForce)[k];
     }
  }
  return nodalForce;

}

/****
void
Domain_opt::getNodalForceMoment(Vector &dsp, double *bcx, int fileNumber, int typ)
{

  // Initialize nodal force
  Vector nodalForce(numnodes,0.0);

  if(!elDisp)  elDisp    = new Vector(maxNumDOFs,0.0);
  if(!elForce) elForce   = new Vector(maxNumDOFs,0.0);

  if (!kelData) kelData  = new double[maxNumDOFs*maxNumDOFs];

  // Loop over all elements
  int iele;
  for(iele=0; iele<numele; ++iele) {

     int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs      = allDOFs->num(iele);
     
     elDisp->zero(); 
     elForce->zero();
    
     FullSquareMatrix karray(numEleDOFs,kelData);  
     
     int k;
     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          (*elDisp)[k] = dsp[cn]; }
        else {
          (*elDisp)[k] = bcx[(*allDOFs)[iele][k]]; }
     }

     karray=packedEset[iele]->stiffness(nodes,karray.data());

     int i,j;
     for (i=0;i<numEleDOFs;i++) {
       for (j=0;j<numEleDOFs;j++) {
         (*elForce)[i] += karray[i][j]*(*elDisp)[j] ;
       }
     }

     for (i=0; i<NodesPerElement; ++i) {
        int node = (*elemToNode)[iele][i];
        for(j=0; j<numEleDOFs; ++j) {
          if ( (*allDOFs)[iele][j] == dsa->locate(node, (1 << typ)) )
            nodalForce[node] += (*elForce)[j];
	}
     }
  }

  geoSource->outputNodeScalars(fileNumber, nodalForce.data(), numnodes, sinfo.dt);
}
***/

//------------------------------------------------------------------------------
/****
void
Domain_opt::buildGlobalK ( SparseMatrix *K,  FullSquareMatrix *kelArray)
{
 K->zeroAll();

 if (!kelData)  kelData  = new double[maxNumDOFs*maxNumDOFs]; 

 FullSquareMatrix karray(1,kelData);

 int iele;
 for(iele=0; iele < numele; ++iele) {
   if(kelArray) {
     karray = kelArray[iele];
   } else {
     karray = packedEset[iele]->stiffness(nodes, karray.data());
   }
   K->add(karray,(*allDOFs)[iele]);  
 }
}
***/

//-----------------------------------------------------------------------------

void
Domain_opt::buildDPressureForceDs(Vector& pseudoforce, GeomState *gs)
{
  Vector ElementDfDs(maxNumDOFs,0.0);

  int numActelm = actInf->size();

  int cflg = geoSource->consistentPFlag();

  int i,iele,actelm;
  for(actelm=0; actelm < numActelm; ++actelm) {

    iele = actInf->getEntry(actelm);

    if ( ! packedEset[iele]->getPressure() ) continue;
  
    ElementDfDs.zero();

    dynamic_cast<Element_opt*>(packedEset[iele])->computeDfDsPressure(nodes, *gradnodes, ElementDfDs,gs,cflg);

    // Assemble element forces into global psuedoforce vector

    for(i=0; i<allDOFs->num(iele); ++i) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][i]);
      if(cn >= 0)
        pseudoforce[cn] += ElementDfDs[i];           
    }
  }
}

//-----------------------------------------------------------------------------
/***
void
Domain_opt::buildGradNodalDensity(double *gnodDen)
{

  //fill projector with appropriate values
  int inode,jele,nconele,j;
  int nodDensFunc      = sinfo.nodDensFunc;
  double nodDensParam1 = sinfo.nodDensParam1;
  double nodDensParam2 = sinfo.nodDensParam2;
  double  dummy,gummy,tempVal;
  int count;

  switch (nodDensFunc) {  
  
  case 1:
   // nodalDen = (sum_i rho_i^p )^(1/p)
   for(inode=0;inode < numnodes;inode++){
      nconele = nodeToElem->num(inode);
      dummy = 0.0;
      gummy = 0.0;
      count = 0;
      for (j=0;j<nconele;j++){
        jele = (*nodeToElem)[inode][j];
        tempVal = packedEset[jele]->getDensity();
        if (tempVal > 0.0){
           dummy  += pow(tempVal,nodDensParam1);
           gummy  += nodDensParam1 * packedEset[jele]->getGradDensity() *
	             pow(tempVal,nodDensParam1-1);
           count++;
        }
      }
      if (count) gnodDen[inode]= gummy * pow(dummy,(1/nodDensParam1)-1) / nodDensParam1;    
    }
    break;
  case 2:
   // nodalDen = 1/N*(sum_i rho_i^p )^(1/p)
   for(inode=0;inode < numnodes;inode++){
      nconele = nodeToElem->num(inode);
      dummy = 0.0;
      gummy = 0.0;
      count = 0;
      for (j=0;j<nconele;j++){
        jele = (*nodeToElem)[inode][j];
        tempVal = packedEset[jele]->getDensity();
        if (tempVal > 0.0){
           dummy += pow(tempVal,nodDensParam1);
	   gummy += nodDensParam1 * packedEset[jele]->getGradDensity() *
	            pow(tempVal,nodDensParam1-1);
           count++;
        }   
      }
      if (count) gnodDen[inode]= gummy * pow(dummy,(1/nodDensParam1)-1)/(nodDensParam1*count);    
    }
    break;
     fprintf(stderr,"\n ERROR Grad Nodal Density Type not defined\n");
     exit(-1);
  }
}
***/

 /***
//-----------------------------------------------------------------------------
double
Domain_opt::getgradNodalTemp(Vector& gradTemp, double* bcx, int node)
{

  int tloc  = c_dsa->locate( node, DofSet::Temp);
  int tloc1 =   dsa->locate( node, DofSet::Temp);
  double res  = (tloc >= 0) ? gradTemp[tloc] : 0.0;
  
  return res;
}

//-------------------------------------------------------------------------

void
Domain_opt::getgradduNodalTemp(Vector& temp, Vector& adj, double *bcx, int node)
{
  int tloc  = c_dsa->locate( node, DofSet::Temp);
  adj[tloc] = 1.0;
  
}
 ***/

//--------------------------------------------------------------------------
/**
void
Domain_opt::computeTHadjLoad(Vector &adj,Vector &load,GeomState *gs, 
                         Corotator** allCorot,double convLambda)
{  
  load.zero();
  int iele,inode,NodesPerElement,globalNode,cn,row,col,k; 

  // Allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  //element temperature information
  if(elNodTemp == 0) elNodTemp = new Vector(maxNumNodes,0.0);
  if(gradelNodTemp == 0) gradelNodTemp = new Vector(maxNumNodes,0.0);

  gradelNodTemp->zero();

  double *nodalTemperatures = 0;
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

  //elemetal load vector 
  if (elForce == 0) elForce = new Vector(maxNumDOFs,0.0); 
  if (elAdj == 0) elAdj   = new Vector(maxNumDOFs,0.0); 
  if (elGrad == 0) elGrad   = new Vector(maxNumDOFs,0.0);

  
  for (iele=0;iele<numele;iele++){ 
     
     NodesPerElement = elemToNode->num(iele);
     packedEset[iele]->nodes(nodeNumbers);
     elAdj->zero();
     
     // ... DETERMINE ELEMENT ADJOINT VECTOR
     for(k=0; k<allDOFs->num(iele); ++k) {
       cn = c_dsa->getRCN((*allDOFs)[iele][k]);
       if(cn >= 0)
    	 (*elAdj)[k] = adj[cn];
       else
    	 (*elAdj)[k] = 0.0;
     } //k=0
  
     //... DETERMINE ELMENTAL TEMPERATURE VECTOR
     if (nodalTemperatures) { 
       for (inode = 0; inode < NodesPerElement; ++inode) {
         if (nodalTemperatures[nodeNumbers[inode]] == defaultTemp)
           (*elNodTemp)[inode] = packedEset[iele]->prop->Ta;
         else
           (*elNodTemp)[inode] = nodalTemperatures[nodeNumbers[inode]];
       } //inode
     }

     //compute d(thermalforce)/dT*adj
     //need to differentiate between linear and nonlinear cases
     for (inode = 0; inode < NodesPerElement; ++inode) {    
        (*gradelNodTemp)[inode] = 1.0;
        elForce->zero();
        elGrad->zero();
        
        //the linear case
        if (probType() != SolverInfo::NonLinDynam){

           packedEset[iele]->computeDfDsThermal(nodes, *gradnodes, 
                                                *elNodTemp,*gradelNodTemp,
                                                *elGrad,0);
            }
        else { // the nonliner case
	  	     
             if (sinfo.resstrFlag == 1){
	   
                 //portion from temperatures
                 packedEset[iele]->getThermalForce(nodes, *elNodTemp, 
		                                   *elForce, 1, gs,2);
                 
		 (*elForce) *= convLambda;
       
                 packedEset[iele]->computeDfDsThermal(nodes, *gradnodes, *elNodTemp, *gradelNodTemp,
		     		                      *elGrad,1,gs,2);
                 (*elGrad) *= convLambda;			      
       
                 //portion from residual stress
                 packedEset[iele]->getThermalForce(nodes, *elNodTemp , 
		                                   *elDisp, 1, gs,1);

                 (*elForce) += *elDisp;
       
                 packedEset[iele]->computeDfDsThermal(nodes, *gradnodes, 
		                                      *elNodTemp, *gradelNodTemp,
		     			              *elDisp,1,gs,1);
                 (*elGrad) += *elDisp;
                }
             else {

                //compute the thermal linear thermal force on the element level
                packedEset[iele]->getThermalForce(nodes, *elNodTemp,
		                                  *elForce, 1, gs);
     
                //scale elemetal force by lambda
                (*elForce) *= convLambda;
                
                packedEset[iele]->computeDfDsThermal(nodes, *gradnodes, 
						     *elNodTemp, *gradelNodTemp,
		    				     *elGrad,1,gs);
		    		     
                (*elGrad) *= convLambda;
               }
	    	
              allCorot[iele]->getDExternalForceDs(nodes, *gradnodes, *gs, 
		                                  elForce->data(), elGrad->data());			 
             } 
        
        //mulitiply elGrad*elForce and add into load
        load[nodeNumbers[inode]] += (*elGrad)*(*elAdj);
        
        (*gradelNodTemp)[inode] = 0.0;
	
     } //inode



  } //iele
 
 delete [] nodeNumbers; 
 if(sinfo.thermalLoadFlag) delete [] nodalTemperatures;
}
***/

//-------------------------------------------------------------------------
 /***
void
Domain_opt::buildDThermalForceDs(Vector& pseudoforce, GeomState *gs,int flag)
{
  // allocate local memory
  int *nodeNumbers = (int*) dbg_alloca(sizeof(int)*maxNumNodes);
    
  // get nodal temperatures

  double *nodalTemperatures = 0;
  double *gradnodalTemperatures = 0;
  if(sinfo.thermalLoadFlag){
      nodalTemperatures     = getNodalTemperatures();
      gradnodalTemperatures = getGradNodalTemperatures();
  }
  if(sinfo.thermoeFlag >=0) {
       nodalTemperatures     = temprcvd;
       gradnodalTemperatures = gradtemprcvd;
  }
  
  int numActelm,NodesPerElement;
  if (flag ==0) numActelm = actInf->size();
  else          numActelm = numele;
  
  int i,iele,actelm,iNode;
  for(actelm=0; actelm < numActelm; ++actelm) {

    if (flag ==0 ) iele = actInf->getEntry(actelm);
    else           iele = actelm;

    elForce->zero();
    elNodTemp->zero();
    gradelNodTemp->zero();

    NodesPerElement = elemToNode->num(iele);
    packedEset[iele]->nodes(nodeNumbers);

    iNode;
    for (iNode = 0; iNode < NodesPerElement; ++iNode) {
       if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp){
         (*elNodTemp)[iNode] = packedEset[iele]->prop->Ta;
	 (*gradelNodTemp)[iNode] = 0.0;
	 }
       else{
     	 (*elNodTemp)[iNode] = nodalTemperatures[nodeNumbers[iNode]];
         (*gradelNodTemp)[iNode] = gradnodalTemperatures[nodeNumbers[iNode]];
         }
    }

    packedEset[iele]->computeDfDsThermal(nodes, *gradnodes, *elNodTemp, *gradelNodTemp,
                                         *elForce,0,gs);

    // Assemble element forces into global psuedoforce vector

    for (i=0; i<allDOFs->num(iele); ++i) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][i]);
      if(cn >= 0)
        pseudoforce[cn] += (*elForce)[i];           
    }
  }
  if(sinfo.thermalLoadFlag ) delete [] nodalTemperatures;
  if(sinfo.thermalLoadFlag ) delete [] gradnodalTemperatures;
}
**/

//-----------------------------------------------------------------------------

double
Domain_opt::getGradStructureMass(int* eleList, int listSize)
{
  double gradtotmas = 0.0;

  int numActelm = actInf->size();
  double mratio = geoSource->getMRatio();
  
  int iele,actelm;
  for(actelm=0; actelm < numActelm; ++actelm) {

    iele = actInf->getEntry(actelm);

    if ( ! checkList(iele,eleList,listSize) ) continue;

    double gradelementMass = dynamic_cast<Element_opt*>(packedEset[iele])->
      getGradMass(nodes,*gradnodes);
    gradtotmas += gradelementMass;
}
//
//  WATCH: derivatives w.r.t. variable nodal mass not implemented
//  
//  double *nodeCount = (double *) dbg_alloca(sizeof(double)*numnodes);
//  double *nodeMass  = (double *) dbg_alloca(sizeof(double)*numnodes);

//  int n;
//  for(n=0; n<numnodes; ++n)
//    nodeCount[n] = nodeMass[n] = 0.0;

//  // ADD DISCRETE MASS
//  DMassData *current = firstDiMass;
//  while(current != 0) {
//    int n = current->node;
//    nodeMass[n]  += current->diMass;
//    nodeCount[n] += 1;
//    current = current->next;
//  }

//  for(n=0; n<numnodes; ++n) {
//    if(nodeCount[n] == 0.0) continue;
//    if(nodeCount[n] > 6) nodeCount[n] = 6;
//    totmas += nodeMass[n]/nodeCount[n];
//  }

   return gradtotmas;
}

//-----------------------------------------------------------------------------

double 
Domain_opt::getGradMomOfInertia(int* eleList, int listSize, int momFlag)
{
  double dIxx    = 0.0;
  double dIyy    = 0.0;
  double dIzz    = 0.0;

  int *nodeNumbers = new int[maxNumNodes];

  int numActelm = actInf->size();
  double mratio = geoSource->getMRatio();
  int iele,actelm;
  for(actelm=0; actelm < numActelm; ++actelm) {

    iele = actInf->getEntry(actelm);

    if ( ! checkList(iele,eleList,listSize) ) continue;

    double elementMass     = packedEset[iele]->getMass(nodes);
    double gradelementMass = dynamic_cast<Element_opt*>(packedEset[iele])->
      getGradMass(nodes,*gradnodes);

    int numNodesPerElement = packedEset[iele]->numNodes();
    packedEset[iele]->nodes( nodeNumbers );

    double massPerNode     = elementMass/numNodesPerElement;
    double gradmassPerNode = gradelementMass/numNodesPerElement;
    
    int i;
    for(i=0; i<numNodesPerElement; ++i) {

       Node &node  = nodes.getNode( nodeNumbers[i] );
       Node &dnode = gradnodes->getNode( nodeNumbers[i] );

       double x = node.x;
       double y = node.y;
       double z = node.z;

       double dx = dnode.x;
       double dy = dnode.y;
       double dz = dnode.z;

       dIxx += gradmassPerNode*sqrt(y*y + z*z);
       dIyy += gradmassPerNode*sqrt(x*x + z*z);
       dIzz += gradmassPerNode*sqrt(x*x + y*y);

       if (sqrt(y*y + z*z) > 1e-16) dIxx += massPerNode*(y*dy+z*dz)/sqrt(y*y + z*z);
       if (sqrt(x*x + z*z) > 1e-16) dIyy += massPerNode*(x*dx+z*dz)/sqrt(x*x + z*z);
       if (sqrt(x*x + y*y) > 1e-16) dIzz += massPerNode*(x*dx+y*dy)/sqrt(x*x + y*y);

    }
  }

  delete [] nodeNumbers;

  // Add discrete masses

  DMassData *current = firstDiMass;

  while(current != 0) {

    int n = current->node;

    // check if the node is in the model

    if(nodes.exist(n)) { 

      Node &node  = nodes.getNode(n);
      Node &dnode = gradnodes->getNode(n);

      double x = node.x;
      double y = node.y;
      double z = node.z;

      double dx = dnode.x;
      double dy = dnode.y;
      double dz = dnode.z;

      if(current->dof == 0 ) {
        if (sqrt(z*z) > 1e-16) dIyy += (current->diMass)*(z*dz)/sqrt(z*z);
        if (sqrt(y*y) > 1e-16) dIzz += (current->diMass)*(y*dy)/sqrt(y*y);
      }

      if(current->dof == 1) {
        if (sqrt(z*z) > 1e-16) dIxx += (current->diMass)*(z*dz)/sqrt(z*z);
        if (sqrt(x*x) > 1e-16) dIzz += (current->diMass)*(x*dx)/sqrt(x*x);
      }
  
      if(current->dof == 2) {
        if (sqrt(y*y) > 1e-16) dIxx += (current->diMass)*(y*dy)/sqrt(y*y);
        if (sqrt(x*x) > 1e-16) dIyy += (current->diMass)*(x*dx)/sqrt(x*x);
      }
    }

    current = current->next;
  }

  switch (momFlag) {
    case 3:
      return dIxx;
      break;
    case 4:
      return dIyy;
      break;
    case 5:
      return dIzz;
      break;
  }

  return 0;
}

//-----------------------------------------------------------------------------

double Domain_opt::getGradVariationNodalPos(int nodeId, int dir)
{

  Node &orgnod = nodescopy->getNode(nodeId);

  double grad = 0.0;
  
  double dx,dy,dz,val;
  double ddx,ddy,ddz;
  
  switch (dir) {
    
    case 0: 
      grad = (*gradnodes)[nodeId]->x;
      break;
    case 1:
      grad = (*gradnodes)[nodeId]->y;
      break;
    case 2:
      grad = (*gradnodes)[nodeId]->z;
      break;
    case 6:
       dx = nodes[nodeId]->x - orgnod.x ;
       dy = nodes[nodeId]->y - orgnod.y ;
       dz = nodes[nodeId]->z - orgnod.z ;
      ddx = (*gradnodes)[nodeId]->x;
      ddy = (*gradnodes)[nodeId]->y;
      ddz = (*gradnodes)[nodeId]->z;
      val = dx*dx + dy*dy + dz*dz;
      val = sqrt(val);
      if ( val > 1.0e-12 ) grad = (dx*ddx+dy*ddy+dz*ddz)/val;
      break;
    default:
      printf("Derivative of Variation of Nodal Position: Wrong Type !");
      exit(-1);
  }
  
  return grad;   
}

//-----------------------------------------------------------------------------

double
Domain_opt::getGradStrainEnergy(Vector &sol, Vector &grad, double *bcx,
                            FullSquareMatrix *kelArray, 
                            int* eleList, int listSize)
{
  double gradenergy = 0.0;

  int iele,actelm,i;

  for(iele=0; iele < numele; ++iele) { // loop over all elements ("el" means element)

     if ( ! checkList(iele,eleList,listSize) ) continue;

     int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs      = allDOFs->num(iele);
     
     FullSquareMatrix karray(numEleDOFs,kelData.get());
     FullSquareMatrix gradkarray(numEleDOFs,gkelData.get());

     // get element displacements
     
     elDisp->zero();

     for (i=0;i<numEleDOFs;i++) {
       int cn = c_dsa->getRCN((*allDOFs)[iele][i]);
       if(cn >= 0) 
     	(*elDisp)[i] = sol[cn];
       else 
     	(*elDisp)[i] = bcx[(*allDOFs)[iele][i]];
     }       
     
     // get contribution of Kel*elGrad
     elKelD->zero();
     elGrad->zero();

     for (i=0;i<numEleDOFs;i++) {
       int cn = c_dsa->getRCN((*allDOFs)[iele][i]);
       if(cn >= 0) (*elGrad)[i] = grad[cn]; 
     }

     if (kelArray)
       karray=kelArray[iele];
     else
       karray=packedEset[iele]->stiffness(nodes,karray.data());

     karray.multiply(*elGrad,*elKelD,2.0);
          
     // get contribution of dKel*elDisp

     actelm = actInf->check(iele); 

     if(actelm) {
       
       dynamic_cast<Element_opt*>(packedEset[iele])->gradstiffness(nodes,*gradnodes,gradkarray);

       gradkarray.multiply(*elDisp,*elKelD, 1.0);

     }
				   
     gradenergy += 0.5*((*elKelD)*(*elDisp));
  }

  return gradenergy;
}


//-----------------------------------------------------------------------------

double
Domain_opt::getGradPartStrainEnergy(Vector &sol, double *bcx,
                                FullSquareMatrix *kelArray, 
                                int* eleList, int listSize)
{ 
  double gradenergy = 0.0;

  int numActelm = actInf->size();
  
  int iele,actelm,i;

  for(actelm=0; actelm < numActelm; ++actelm) {

     iele = actInf->getEntry(actelm);

     if ( ! checkList(iele,eleList,listSize) ) continue;

     int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs      = allDOFs->num(iele);
     
     FullSquareMatrix gradkarray(numEleDOFs,gkelData.get());

     // get element displacements
     elDisp->zero();

     for (i=0;i<numEleDOFs;i++) {
       int cn = c_dsa->getRCN((*allDOFs)[iele][i]);
       if(cn >= 0) 
     	(*elDisp)[i] = sol[cn];
       else 
     	(*elDisp)[i] = bcx[(*allDOFs)[iele][i]];
     }       
     
     // get contribution of dKel*elDisp
       
     dynamic_cast<Element_opt*>(packedEset[iele])->gradstiffness(nodes,*gradnodes,gradkarray);

     elKelD->zero();

     gradkarray.multiply(*elDisp,*elKelD, 1.0);
				   
     gradenergy += 0.5*((*elKelD)*(*elDisp));
  }

  return gradenergy;
}

//------------------------------------------------------------------------------

void
Domain_opt::getGradduStrainEnergy(Vector &sol, Vector &adj, double *bcx,
                              FullSquareMatrix *kelArray, 
                              int* eleList, int listSize)
{
  for(int iele=0; iele<numele; ++iele) 
    {
      if ( ! checkList(iele,eleList,listSize) ) 
	{ continue; }

      int NodesPerElement = elemToNode->num(iele);
      int numEleDOFs      = allDOFs->num(iele);
      
      FullSquareMatrix karray(numEleDOFs,kelData.get());      
     
      for(int k=0; k<numEleDOFs; ++k) 
	{
	  int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
	  if(cn >= 0) 
	    { (*elDisp)[k] = sol[cn]; }
	  else 
	    { (*elDisp)[k] = bcx[(*allDOFs)[iele][k]]; }
	}
      
      if (kelArray) 
	{ karray=kelArray[iele]; }
      else 
	{ karray=packedEset[iele]->stiffness(nodes,karray.data()); }

      elAdj->zero();
      
      karray.multiply(*elDisp,*elAdj, 1.0);
      
      for(int k=0; k<numEleDOFs; ++k) 
	{
	  int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
	  if(cn >= 0) 
	    { adj[cn] += (*elAdj)[k]; }
	}
    }
}

//------------------------------------------------------------------------------

double
Domain_opt::getGradStressIntegral(Vector &sol, Vector &grad, double *bcx, 
                              double& refStress, double& powFac, 
                              int* eleList, int listSize, 
                              int stressIndex, int areaFlag)
{
  
  int iele, k, i, iNode;

  double vmint, vol, dvmint, dvol;
  double graddvmint, graddvol;

  graddvmint = graddvol = 0.0;

  // Allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  /*
  // Get temperature and derivatives
  double *nodalTemperatures = 0;
  double *gradnodalTemperatures = 0;

  if(sinfo.thermalLoadFlag) {
     nodalTemperatures = getNodalTemperatures();
     gradnodalTemperatures = getGradNodalTemperatures();
  }

  if(sinfo.thermoeFlag >=0) {
     nodalTemperatures = temprcvd;
     gradnodalTemperatures = gradtemprcvd;
  }
  */

  for(iele=0; iele<numele; ++iele) {

     if ( ! checkList(iele,eleList,listSize) ) continue;

     int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs  = allDOFs->num(iele);

     packedEset[iele]->nodes(nodeNumbers);

     elDisp->zero();
     elGrad->zero();

     //elNodTemp->zero();
     //gradelNodTemp->zero();
 
     // ... DETERMINE ELEMENT DISPLACEMENT VECTOR

     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elDisp)[k] =  sol[cn];
        else
          (*elDisp)[k] =  bcx[(*allDOFs)[iele][k]];
     }

     for (i=0;i<numEleDOFs;i++) {
       int cn = c_dsa->getRCN((*allDOFs)[iele][i]);
       if(cn >= 0) (*elGrad)[i] = grad[cn]; 
     }

     /***
     // ... DETERMINE NODAL TEMPERATURE VECTOR
     if (nodalTemperatures) {
       for (iNode = 0; iNode < NodesPerElement; ++iNode) {
         if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp){
           (*elNodTemp)[iNode] = packedEset[iele]->prop->Ta;
	   (*gradelNodTemp)[iNode] = 0.0;
	 }  
         else{
           (*elNodTemp)[iNode] = nodalTemperatures[nodeNumbers[iNode]];
	   (*gradelNodTemp)[iNode] = gradnodalTemperatures[nodeNumbers[iNode]];
	 }  
       }
     }
     ***/

     // ... calculate integral (sig_vm/sigbar)^fac over ele-volume

     vmint = vol = dvmint = dvol = 0.0;

     dynamic_cast<Element_opt*>(packedEset[iele])->
       getGradVonMisesInt(nodes, *gradnodes, *elDisp, *elGrad, 
			  refStress, powFac, areaFlag,
			  vmint, dvmint, vol, dvol);
     //elNodTemp->data(), gradelNodTemp->data());

     graddvmint += dvmint;
     graddvol   +=   dvol;
  }

  delete [] nodeNumbers;
  /*if(sinfo.thermalLoadFlag) delete [] nodalTemperatures;
    if(sinfo.thermalLoadFlag) delete [] gradnodalTemperatures;*/

  if (!areaFlag) {
    valvol   = 1.0;
    graddvol = 0.0;
  }

  if (valvmint != 0) {
    double fac = pow(valvmint/valvol,1.0/powFac-1.0)/powFac;
    return  fac*(graddvmint/valvol - graddvol*valvmint/(valvol*valvol));
  }
  else return 0;
  
}

//------------------------------------------------------------------------------

double
Domain_opt::getGradPartStressIntegral(Vector &sol, Vector &adj, double *bcx,
                                double& refStress, double& powFac, 
                                int* eleList, int listSize, 
                                int stressIndex, int areaFlag)
				
{  
  int k, actelm, iele,iNode;

  double vmint, vol, dvmint, dvol;
  double graddvmint, graddvol;
  
  int numActelm = actInf->size();

  // Allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  /*
  double *nodalTemperatures = 0;
  double *gradnodalTemperatures = 0;

  if(sinfo.thermalLoadFlag){
      nodalTemperatures     = getNodalTemperatures();
      gradnodalTemperatures = getGradNodalTemperatures();
  }

  if(sinfo.thermoeFlag >=0) {
       nodalTemperatures     = temprcvd;
       gradnodalTemperatures = gradtemprcvd;
  }
  */

  graddvmint= graddvol = 0.0;
  elGrad->zero();

 for(actelm=0; actelm < numActelm; ++actelm) {  

     iele = actInf->getEntry(actelm);

     if ( ! checkList(iele,eleList,listSize) ) continue;

     elDisp->zero();
     /*elNodTemp->zero();
       gradelNodTemp->zero();*/

     int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs  = allDOFs->num(iele);

     packedEset[iele]->nodes(nodeNumbers);

     // ... DETERMINE ELEMENT DISPLACEMENT VECTOR

     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elDisp)[k] =  sol[cn];
        else
          (*elDisp)[k] =  bcx[(*allDOFs)[iele][k]];
     }

     /***
     // ... DETERMINE NODAL TEMPERATURE VECTOR
     if (nodalTemperatures) {
       for (iNode = 0; iNode < NodesPerElement; ++iNode) {
         if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp){
           (*elNodTemp)[iNode] = packedEset[iele]->prop->Ta;
	   (*gradelNodTemp)[iNode] = 0.0;
	 }  
         else{
           (*elNodTemp)[iNode] = nodalTemperatures[nodeNumbers[iNode]];
	   (*gradelNodTemp)[iNode] = gradnodalTemperatures[nodeNumbers[iNode]];
	 }  
        }
      }
     ***/

     // ... calculate integral (sig_vm/sigbar)^fac over ele-volume

     vmint = vol = dvmint = dvol = 0.0;

     dynamic_cast<Element_opt*>(packedEset[iele])->getGradVonMisesInt(nodes, *gradnodes, *elDisp, *elGrad, 
                                          refStress, powFac, areaFlag,
                                          vmint, dvmint, vol, dvol);
     //elNodTemp->data(),
     //gradelNodTemp->data());

     graddvmint += dvmint;
     graddvol   +=   dvol;

  }

  delete [] nodeNumbers;
  //if(sinfo.thermalLoadFlag) delete [] nodalTemperatures;
  //if(sinfo.thermalLoadFlag) delete [] gradnodalTemperatures;

  if (!areaFlag) {
    valvol   = 1.0;
    graddvol = 0.0;
  }

  if (valvmint != 0 ) {
    double fac = pow(valvmint/valvol,1.0/powFac-1.0)/powFac;
    return  fac*(graddvmint/valvol - graddvol*valvmint/(valvol*valvol));
  } 
  else return 0;
}

//------------------------------------------------------------------------------

void
Domain_opt::getGradduStressIntegral(Vector &sol, Vector &adj, double *bcx,
                                double& refStress, double& powFac, 
                                int* eleList, int listSize, 
                                int stressIndex, int areaFlag)
{ 
  int k, iele, iNode;
  double vmint, vol, dvmint, dvol;

  // Allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  //double *nodalTemperatures = 0;
  //double *gradnodalTemperatures = 0;

//   if(sinfo.thermalLoadFlag){
//       nodalTemperatures     = getNodalTemperatures();
//       gradnodalTemperatures = getGradNodalTemperatures();
//   }

//   if(sinfo.thermoeFlag >=0) {
//        nodalTemperatures     = temprcvd;
//        gradnodalTemperatures = gradtemprcvd;
//   }

  elGrad->zero();

  for(iele=0; iele<numele; ++iele) {

    if ( ! checkList(iele,eleList,listSize) ) continue;

    elDisp->zero();
//     elNodTemp->zero();
//     gradelNodTemp->zero();

    int NodesPerElement = elemToNode->num(iele);
    int ldofs	       = allDOFs->num(iele);

     packedEset[iele]->nodes(nodeNumbers);

    // ... DETERMINE ELEMENT DISPLACEMENT VECTOR

    for(k=0; k<ldofs; ++k) {
       int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
       if(cn >= 0) {
    	 (*elDisp)[k] = sol[cn]; }
       else {
    	 (*elDisp)[k] = bcx[(*allDOFs)[iele][k]]; }
    }

    // ... DETERMINE NODAL TEMPERATURE VECTOR
//    if (nodalTemperatures){
//      for (iNode = 0; iNode < NodesPerElement; ++iNode) {
//        if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp){
//          (*elNodTemp)[iNode] = packedEset[iele]->prop->Ta;
// 	 (*gradelNodTemp)[iNode] = 0.0;
// 	 }
//        else{
//      	 (*elNodTemp)[iNode] = nodalTemperatures[nodeNumbers[iNode]];
//          (*gradelNodTemp)[iNode] = gradnodalTemperatures[nodeNumbers[iNode]];
//        }
//      }
//    }

    // ... PERTURB EACH ELEMENTAL DOF

    int i;
    for(i=0; i<ldofs; ++i) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][i]);
      if(cn >= 0) {
    	(*elGrad)[i] = 1.0;

    	dvmint= vmint = vol = dvol = 0.0;

    	dynamic_cast<Element_opt*>(packedEset[iele])->
	  getGradVonMisesInt(nodes, *gradnodes, *elDisp,
			     *elGrad, refStress, powFac, areaFlag,
			     vmint, dvmint, vol, dvol);
	//elNodTemp->data(), gradelNodTemp->data());

    	 adj[cn] += dvmint;
    	(*elGrad)[i] = 0.0;
      }
    }
  }

  delete [] nodeNumbers;
//   if(sinfo.thermalLoadFlag) delete [] nodalTemperatures;
//   if(sinfo.thermalLoadFlag) delete [] gradnodalTemperatures;

  double scale;
  if (valvmint != 0) 
      scale = pow(valvmint/valvol,1.0/powFac-1.0)/powFac/valvol; 
  else 
      scale = 0.0;

  int i;
  for(i=0;i<adj.size();i++)  adj[i] *= scale;
}

//------------------------------------------------------------------------------

double
Domain_opt::getGradDCmass(Vector &sol, Vector &grad, double *bcx, double& powFac, int* eleList, int listSize)
{
  int k, iele;

  double gradDCmass = 0.0;

  for(iele=0; iele<numele; ++iele) {

     if ( ! checkList(iele,eleList,listSize) ) continue;

     elDisp->zero();
     elGrad->zero();

     for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          (*elDisp)[k] =  sol[cn];
	  (*elGrad)[k] = grad[cn];
	}
        else
          (*elDisp)[k] =  bcx[(*allDOFs)[iele][k]];
     }

     double edcm = dynamic_cast<Element_opt*>(packedEset[iele])->
       getGradDCmass(nodes, *gradnodes, *elDisp, *elGrad, powFac);

     gradDCmass += edcm;
  }      

  return gradDCmass;
}

//------------------------------------------------------------------------------

double
Domain_opt::getGradPartDCmass(Vector &sol, double *bcx, double& powFac, int* eleList, int listSize)
{
  int k, actelm, iele;

  int numActelm = actInf->size();

  double gradDCmass = 0.0;

  elGrad->zero();

  for(actelm=0; actelm < numActelm; ++actelm) {  

     iele = actInf->getEntry(actelm);

     if ( ! checkList(iele,eleList,listSize) ) continue;

     elDisp->zero();

     for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0)
          (*elDisp)[k] =  sol[cn];
        else
          (*elDisp)[k] =  bcx[(*allDOFs)[iele][k]];
     }

     double edcm = dynamic_cast<Element_opt*>(packedEset[iele])->
       getGradDCmass(nodes, *gradnodes, *elDisp, *elGrad, powFac);

     gradDCmass += edcm;
  }      

  return gradDCmass;
}

//------------------------------------------------------------------------------

void
Domain_opt::getGradduDCmass(Vector &sol, Vector &adj, double *bcx, double& powFac, int* eleList, int listSize)
{
  int i, k, iele;

  elGrad->zero();

  for(iele=0; iele<numele; ++iele) {

     if ( ! checkList(iele,eleList,listSize) ) continue;

     elDisp->zero();
     elGrad->zero();

     int ldofs = allDOFs->num(iele);
 
     for(k=0; k<ldofs; ++k) {

        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);

        if(cn >= 0)
          (*elDisp)[k] =  sol[cn];
        else
          (*elDisp)[k] =  bcx[(*allDOFs)[iele][k]];
     }

     for(i=0; i<ldofs; ++i) {
 
        int cn = c_dsa->getRCN((*allDOFs)[iele][i]);

        if(cn >= 0) {

      	   (*elGrad)[i] = 1.0;

           double edcm = dynamic_cast<Element_opt*>(packedEset[iele])->
	     getGradDCmass(nodes, *gradnodes, *elDisp, *elGrad, powFac);
    	 
	   adj[cn] += edcm;
    	 
	  (*elGrad)[i] = 0.0;
        }
     }
  }      
}

//------------------------------------------------------------------------------

double
Domain_opt::getGradNodalStressStrain(Vector &sol,Vector &grad, double *bcx,
                                 int node, int stressIndex, int surface, 
                                 int eleFlag)
{
  // Allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  // check for elemental evaluation; Number of Connected Elements
  int nconele;
  if (eleFlag >= 0) 
    nconele = 1;
  else
    nconele = nodeToElem->num(node);

//   double *nodalTemperatures = 0;
//   double *gradnodalTemperatures = 0;

//   if(sinfo.thermalLoadFlag){
//       nodalTemperatures     = getNodalTemperatures();
//       gradnodalTemperatures = getGradNodalTemperatures();
//   }

//   if(sinfo.thermoeFlag >=0) {
//        nodalTemperatures     = temprcvd;
//        gradnodalTemperatures = gradtemprcvd;
//   }

  double gradstress = 0.0;
  double weight     = 0.0;

  int icon, iNode, k, iele;
  for(icon=0; icon<nconele; ++icon) {
 
     elDisp->zero();
     elGrad->zero();
//      elNodTemp->zero();
//      gradelNodTemp->zero();
     elstress->zero();
     elweight->zero();

     if (eleFlag >= 0) 
       iele=eleFlag;
     else
       iele = (*nodeToElem)[node][icon];

     int NodesPerElement = elemToNode->num(iele);

     packedEset[iele]->nodes(nodeNumbers);

     // ... DETERMINE ELEMENT DISPLACEMENT VECTOR

     for(k=0; k<allDOFs->num(iele); ++k) {
       int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
       if(cn >= 0) {
          (*elDisp)[k] = sol[cn];
          (*elGrad)[k] = grad[cn]; }
       else {
         (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
         (*elGrad)[k] = 0.0; }
     }

     // ... DETERMINE NODAL TEMPERATURE VECTOR
//      if (nodalTemperatures) {
//        for (iNode = 0; iNode < NodesPerElement; ++iNode) {
//          if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp){
//            (*elNodTemp)[iNode] = packedEset[iele]->prop->Ta;
// 	   (*gradelNodTemp)[iNode] = 0.0;
//          }
//          else{
//      	   (*elNodTemp)[iNode] = nodalTemperatures[nodeNumbers[iNode]];
//            (*gradelNodTemp)[iNode] = gradnodalTemperatures[nodeNumbers[iNode]];
//          }
//        }     
//      }
     
     // ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT

     dynamic_cast<Element_opt*>(packedEset[iele])->
       getGradVonMises(*gradelstress, *elweight, nodes, 
		       *gradnodes, *elDisp, *elGrad, 
		       stressIndex, surface);
     //elNodTemp->data(), gradelNodTemp->data());

     // ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT

     // ... eleFlag <  0: average stress in node
     // ... eleFlag >= 0: average stress of particular element

     if (eleFlag < 0) {
       for(k=0; k<NodesPerElement; ++k) { 
         int actnod = (*elemToNode)[iele][k];
         if (actnod == node) { 
           gradstress += (*gradelstress)[k]*(*elweight)[k];
           weight     += (*elweight)[k];
         }
       }
     }
     else {
       for(k=0; k<NodesPerElement; ++k) { 
         gradstress += (*gradelstress)[k]*(*elweight)[k];
         weight     += (*elweight)[k];
       }
     }
  }

  delete [] nodeNumbers;
//   if(sinfo.thermalLoadFlag) delete [] nodalTemperatures;
//   if(sinfo.thermalLoadFlag) delete [] gradnodalTemperatures;

  return weight != 0 ? gradstress/weight : 0.0;  
}

//------------------------------------------------------------------------------

double
Domain_opt::getGradPartNodalStressStrain(Vector &sol, double *bcx,
                                     int node, int stressIndex, int surface, 
                                     int eleFlag)
{
  // Allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  // check for elemental evaluation; Number of Connected Elements
  int nconele;
  if (eleFlag >= 0) 
    nconele = 1;
  else
    nconele = nodeToElem->num(node);

//   double *nodalTemperatures = 0;
//   double *gradnodalTemperatures = 0;

//   if(sinfo.thermalLoadFlag){
//       nodalTemperatures     = getNodalTemperatures();
//       gradnodalTemperatures = getGradNodalTemperatures();
//   }

//   if(sinfo.thermoeFlag >=0) {
//        nodalTemperatures     = temprcvd;
//        gradnodalTemperatures = gradtemprcvd;
//   }

  double gradstress = 0.0;
  double weight     = 0.0;

  elGrad->zero();

  int icon,k,iNode,iele;
  for(icon=0; icon<nconele; ++icon) {

     elDisp->zero();
//      elNodTemp->zero();
//      gradelNodTemp->zero();
     elstress->zero();
     elweight->zero();
 
     if (eleFlag >= 0) 
       iele=eleFlag;
     else
       iele = (*nodeToElem)[node][icon];

     int NodesPerElement = elemToNode->num(iele);

     packedEset[iele]->nodes(nodeNumbers);

     // ... DETERMINE ELEMENT DISPLACEMENT VECTOR

     for(k=0; k<allDOFs->num(iele); ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) 
          (*elDisp)[k] = sol[cn];
        else 
          (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
     }

     // ... DETERMINE NODAL TEMPERATURE VECTOR       
//      if (nodalTemperatures){
//        for (iNode = 0; iNode < NodesPerElement; ++iNode) {
//          if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp){
//            (*elNodTemp)[iNode] = packedEset[iele]->prop->Ta;
// 	   (*gradelNodTemp)[iNode] = 0.0;
//          }
//          else{
//      	   (*elNodTemp)[iNode] = nodalTemperatures[nodeNumbers[iNode]];
//            (*gradelNodTemp)[iNode] = gradnodalTemperatures[nodeNumbers[iNode]];
//          }
//        }
//      }

     // ... CALCULATE STRESS/STRAIN VALUE FOR EACH NODE OF THE ELEMENT

     dynamic_cast<Element_opt*>(packedEset[iele])->
       getGradVonMises(*gradelstress, *elweight, nodes, 
		       *gradnodes, *elDisp, *elGrad, 
		       stressIndex, surface);
     //elNodTemp->data(), gradelNodTemp->data());

     // ... ASSEMBLE ELEMENT'S NODAL STRESS/STRAIN & WEIGHT
     // ... eleFlag <  0: average stress in node
     // ... eleFlag >= 0: average stress of particular element

     if (eleFlag < 0) {
       for(k=0; k<NodesPerElement; ++k) { 
         int actnod = (*elemToNode)[iele][k];
         if (actnod == node) { 
           gradstress += (*gradelstress)[k]*(*elweight)[k];
           weight += (*elweight)[k];
         }
       }
     }
     else {
       for(k=0; k<NodesPerElement; ++k) { 
         gradstress += (*gradelstress)[k]*(*elweight)[k];
         weight += (*elweight)[k];
       }
     }
  }

  delete [] nodeNumbers;
//   if(sinfo.thermalLoadFlag) delete [] nodalTemperatures;
//   if(sinfo.thermalLoadFlag) delete [] gradnodalTemperatures;

  return weight != 0 ? gradstress/weight : 0.0;  
}


//------------------------------------------------------------------------------

void
Domain_opt::getGradduNodalStressStrain(Vector &sol,Vector &adj, double *bcx,
                                   int node, int stressIndex, int surface, 
                                   int eleFlag)
{
  // Allocate integer array to store node numbers
  int *nodeNumbers = new int[maxNumNodes];

  // check for elemental evaluation; Number of Connected Elements
  int nconele;
  if (eleFlag >= 0) 
    nconele = 1;
  else
    nconele = nodeToElem->num(node);

//   double *nodalTemperatures = 0;
//   double *gradnodalTemperatures = 0;

//   if(sinfo.thermalLoadFlag){
//       nodalTemperatures     = getNodalTemperatures();
//       gradnodalTemperatures = getGradNodalTemperatures();
//   }

//   if(sinfo.thermoeFlag >=0) {
//        nodalTemperatures     = temprcvd;
//        gradnodalTemperatures = gradtemprcvd;
//   }

  int    icount  = 0;
  double cweight = 0.0;

  elGrad->zero();

  int icon,k,iNode,iele;
  for(icon=0; icon<nconele; ++icon) {
 
     elDisp->zero();
//      elNodTemp->zero();
//      gradelNodTemp->zero();
     elstress->zero();
     elweight->zero();

     if (eleFlag >= 0) 
       iele=eleFlag;
     else
       iele = (*nodeToElem)[node][icon];

     int NodesPerElement = elemToNode->num(iele);
     int ldofs           = allDOFs->num(iele);

     packedEset[iele]->nodes(nodeNumbers);

     // ... DETERMINE ELEMENT DISPLACEMENT VECTOR

     for(k=0; k<ldofs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          (*elDisp)[k] = sol[cn]; }
        else {
          (*elDisp)[k] = bcx[(*allDOFs)[iele][k]]; }
     }

     // ... DETERMINE NODAL TEMPERATURE VECTOR

//      if (nodalTemperatures){
//        for (iNode = 0; iNode < NodesPerElement; ++iNode) {
//          if (nodalTemperatures[nodeNumbers[iNode]] == defaultTemp){
//            (*elNodTemp)[iNode] = packedEset[iele]->prop->Ta;
// 	   (*gradelNodTemp)[iNode] = 0.0;
//          }
//          else{
//      	   (*elNodTemp)[iNode] = nodalTemperatures[nodeNumbers[iNode]];
//            (*gradelNodTemp)[iNode] = gradnodalTemperatures[nodeNumbers[iNode]];
//          }
//        }
//      }

     // ... PERTURB EACH ELEMENTAL DOF

     int i;
     for(i=0; i<ldofs; ++i) {
       int cn = c_dsa->getRCN((*allDOFs)[iele][i]);
       if(cn >= 0) {
         (*elGrad)[i] = 1.0;
          
         dynamic_cast<Element_opt*>(packedEset[iele])->
	   getGradVonMises(*gradelstress, *elweight,
			   nodes, *gradnodes, *elDisp, 
			   *elGrad, stressIndex, surface);
//                                   elNodTemp->data(), gradelNodTemp->data());

         (*elGrad)[i] = 0.0;

         // ... eleFlag <  0: average stress in node
         // ... eleFlag >= 0: average stress of particular element

         if (eleFlag < 0) {
           for(k=0; k<NodesPerElement; ++k) { 
             int actnod = (*elemToNode)[iele][k];
             if (actnod == node) { 
               adj[cn] += (*gradelstress)[k]*(*elweight)[k];
               if ( icount == icon ) { 
         	 cweight += (*elweight)[k];
         	 icount++;
               }
             } 
           }
         }
         else {
           for(k=0; k<NodesPerElement; ++k) { 
             adj[cn] += (*gradelstress)[k]*(*elweight)[k];
             if ( icount == icon ) cweight += (*elweight)[k];
           }
	   if ( icount == icon ) icount++;
         }
       }
     }
  }

  if ( cweight != 0.0 ) adj *= (1.0/cweight);

  delete [] nodeNumbers;
//   if(sinfo.thermalLoadFlag) delete [] nodalTemperatures;
//   if(sinfo.thermalLoadFlag) delete [] gradnodalTemperatures;
}


//------------------------------------------------------------------------------
void
Domain_opt::getGradduNodalDisp(Vector &sol, double *bcx, int node, int dispTyp, 
                           double* nodFac, int* nodLoc, int& numNod)
{
    int loc,loc1;
    int loc11,loc22,loc33;
    double d1,d2,d3;
    double sum = 0.0;

    numNod=0;
        
    switch (dispTyp) {
              
    case 0: 
      loc  = c_dsa->locate( node, DofSet::Xdisp);
      if ( loc >= 0  ) { numNod=1; nodFac[0]=1.0; nodLoc[0]=loc; } 
      break;
    case 1: 
      loc  = c_dsa->locate( node, DofSet::Ydisp);
      if ( loc >= 0  ) { numNod=1; nodFac[0]=1.0; nodLoc[0]=loc; } 
      break;
    case 2: 
      loc  = c_dsa->locate( node, DofSet::Zdisp);
      if ( loc >= 0  ) { numNod=1; nodFac[0]=1.0; nodLoc[0]=loc; } 
      break;
    case 3: 
      loc  = c_dsa->locate( node, DofSet::Xrot);
      if ( loc >= 0  ) { numNod=1; nodFac[0]=1.0; nodLoc[0]=loc;  } 
      break;
    case 4: 
      loc  = c_dsa->locate( node, DofSet::Yrot);
      if ( loc >= 0  ) { numNod=1; nodFac[0]=1.0; nodLoc[0]=loc; } 
      break;
    case 5: 
      loc  = c_dsa->locate( node, DofSet::Zrot);
      if ( loc >= 0  ) { numNod=1; nodFac[0]=1.0; nodLoc[0]=loc; } 
      break;
    case 6:
      loc11 = c_dsa->locate( node, DofSet::Xdisp);
      loc1  =   dsa->locate( node, DofSet::Xdisp);
      if      ( loc11 >= 0 ) { d1 = sol[loc11]; } 
      else if ( loc1  >= 0 ) { d1 = bcx[loc1]; }
      else                   { d1 = 0.0; }
      sum  += d1*d1;
      loc22 = c_dsa->locate( node, DofSet::Ydisp);
      loc1  =   dsa->locate( node, DofSet::Ydisp);
      if      ( loc22 >= 0 ) { d2 = sol[loc22]; }  
      else if ( loc1  >= 0 ) { d2 = bcx[loc1]; } 
      else                   { d2 = 0.0; }
      sum  += d2*d2;
      loc33 = c_dsa->locate( node, DofSet::Zdisp);
      loc1  =   dsa->locate( node, DofSet::Zdisp);
      if      ( loc33 >= 0 ) { d3 = sol[loc33]; }  
      else if ( loc1  >= 0 ) { d3 = bcx[loc1]; } 
      else                   { d3 = 0.0; }
      sum += d3*d3;
      sum  = sqrt(sum);
      numNod=0;
      if ( loc11 >= 0 ) { nodFac[numNod]= d1/sum; nodLoc[numNod]=loc11; numNod++; }
      if ( loc22 >= 0 ) { nodFac[numNod]= d2/sum; nodLoc[numNod]=loc22; numNod++; }
      if ( loc33 >= 0 ) { nodFac[numNod]= d3/sum; nodLoc[numNod]=loc33; numNod++; }
      break;
    case 7:
      loc  = c_dsa->locate( node, DofSet::Temp);
      if ( loc >= 0  ) { numNod=1; nodFac[0]=1.0; nodLoc[0]=loc; } 
      break; 
    default:
      printf("OptcritDisp::evaluate - Wrong Disptype !");
    }
}


//------------------------------------------------------------------------------
void Domain_opt::getGradduNodalDisp(ComplexVector &sol, DComplex *bcx, int node, int dispTyp, 
				    DComplex* nodFac, int* nodLoc, int& numNod)
{
  DComplex one(1.0, 1.0);

    int loc,loc1;
    int loc11,loc22,loc33;
    DComplex d1,d2,d3;
    double sum = 0.0;

    numNod=0;
        
    switch (dispTyp) {
              
    case 0: 
      loc  = c_dsa->locate( node, DofSet::Xdisp);
      if ( loc >= 0  ) 
	{ numNod=1; nodFac[0]= ScalarTypes::norm(sol[loc])>0?
	    ScalarTypes::conj(sol[loc])/ScalarTypes::norm(sol[loc]) : 0;
	  nodLoc[0]=loc; }
      break;
    case 1: 
      loc  = c_dsa->locate( node, DofSet::Ydisp);
      if ( loc >= 0  )
	{ numNod=1; nodFac[0]= ScalarTypes::norm(sol[loc])>0?
	    ScalarTypes::conj(sol[loc])/ScalarTypes::norm(sol[loc]) : 0;
	  nodLoc[0]=loc; }
      break;
    case 2: 
      loc  = c_dsa->locate( node, DofSet::Zdisp);
      if ( loc >= 0  )
	{ numNod=1; nodFac[0]= ScalarTypes::norm(sol[loc])>0?
	    ScalarTypes::conj(sol[loc])/ScalarTypes::norm(sol[loc]) : 0;
	  nodLoc[0]=loc; }
      break;
    case 3: 
      loc  = c_dsa->locate( node, DofSet::Xrot);
      if ( loc >= 0  )
	{ numNod=1; nodFac[0]= ScalarTypes::norm(sol[loc])>0?
	    ScalarTypes::conj(sol[loc])/ScalarTypes::norm(sol[loc]) : 0;
	  nodLoc[0]=loc; }
      break;
    case 4: 
      loc  = c_dsa->locate( node, DofSet::Yrot);
      if ( loc >= 0  )
	{ numNod=1; nodFac[0]= ScalarTypes::norm(sol[loc])>0?
	    ScalarTypes::conj(sol[loc])/ScalarTypes::norm(sol[loc]) : 0;
	  nodLoc[0]=loc; }
      break;
    case 5: 
      loc  = c_dsa->locate( node, DofSet::Zrot);
      if ( loc >= 0  )
	{ numNod=1; nodFac[0]= ScalarTypes::norm(sol[loc])>0?
	    ScalarTypes::conj(sol[loc])/ScalarTypes::norm(sol[loc]) : 0;
	  nodLoc[0]=loc; }
      break;
    case 6:
      loc11 = c_dsa->locate( node, DofSet::Xdisp);
      loc1  =   dsa->locate( node, DofSet::Xdisp);
      if      ( loc11 >= 0 ) { d1 = sol[loc11]; } 
      else if ( loc1  >= 0 ) { d1 = bcx[loc1]; }
      else                   { d1 = 0.0; }
      sum  += ScalarTypes::sqNorm(d1);
      loc22 = c_dsa->locate( node, DofSet::Ydisp);
      loc1  =   dsa->locate( node, DofSet::Ydisp);
      if      ( loc22 >= 0 ) { d2 = sol[loc22]; }  
      else if ( loc1  >= 0 ) { d2 = bcx[loc1]; } 
      else                   { d2 = 0.0; }
      sum  += ScalarTypes::sqNorm(d2);
      loc33 = c_dsa->locate( node, DofSet::Zdisp);
      loc1  =   dsa->locate( node, DofSet::Zdisp);
      if      ( loc33 >= 0 ) { d3 = sol[loc33]; }  
      else if ( loc1  >= 0 ) { d3 = bcx[loc1]; } 
      else                   { d3 = 0.0; }
      sum += ScalarTypes::sqNorm(d3);
      sum  = sqrt(sum);
      numNod=0;
      if ( loc11 >= 0 && sum>0) { nodFac[numNod]= ScalarTypes::conj(d1/sum); nodLoc[numNod]=loc11; numNod++; }
      if ( loc22 >= 0 && sum>0) { nodFac[numNod]= ScalarTypes::conj(d2/sum); nodLoc[numNod]=loc22; numNod++; }
      if ( loc33 >= 0 && sum>0) { nodFac[numNod]= ScalarTypes::conj(d3/sum); nodLoc[numNod]=loc33; numNod++; }
      break;
    default:
      printf("OptcritDisp::evaluate - Wrong Disptype !");
    }
}

//------------------------------------------------------------------------------

double
Domain_opt::getGradNodalInternalForce(Vector &dsp, Vector &grad, double *bcx, int node, int typ)
{
  // Initialize nodal force
  double gradnodalForce = 0;

  // Number of connected elements
  int nconele = nodeToElem->num(node);

  // Extract the relevant DOF
  int dof  = dsa->locate(node, (1 << typ));
  
  //if(elGrad.get() == 0)    elGrad.reset(new Vector(maxNumDOFs,0.0));
  
  // Loop over all connected elements
  int icon;
  for(icon=0; icon<nconele; ++icon) {

     int iele = (*nodeToElem)[node][icon];
       
     int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs      = allDOFs->num(iele);
     
     elDisp->zero(); 
     elGrad->zero();
     elForce->zero();
    
     FullSquareMatrix     karray(numEleDOFs,kelData.get());  
     FullSquareMatrix gradkarray(numEleDOFs,gkelData.get());  
     
     int k;
     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          (*elDisp)[k] = dsp[cn];
          (*elGrad)[k] = grad[cn]; }
        else {
         (*elDisp)[k] = bcx[(*allDOFs)[iele][k]];
         (*elGrad)[k] = 0.0; }
     }

     karray=packedEset[iele]->stiffness(nodes,karray.data());

     dynamic_cast<Element_opt*>(packedEset[iele])->
       gradstiffness(nodes,*gradnodes,gradkarray);

     int i,j;
     for (i=0;i<numEleDOFs;i++) {
       for (j=0;j<numEleDOFs;j++) {
         (*elForce)[i] += karray[i][j]*(*elGrad)[j] + gradkarray[i][j]*(*elDisp)[j]  ;
       }
     }

     for(k=0; k<numEleDOFs; ++k) {
       if ( (*allDOFs)[iele][k] == dof ) gradnodalForce += (*elForce)[k];
     }
  }
  return gradnodalForce;

}

//------------------------------------------------------------------------------

double
Domain_opt::getGradPartNodalInternalForce(Vector &dsp, double *bcx, int node, int typ)
{
  // Initialize nodal force
  double gradnodalForce = 0;

  // Number of connected elements
  int nconele = nodeToElem->num(node);

  // Extract the relevant DOF
  int dof = dsa->locate(node, (1 << typ));


  // Loop over all connected elements

  int icon;
  for(icon=0; icon<nconele; ++icon) {

     int iele = (*nodeToElem)[node][icon];
       
     int NodesPerElement = elemToNode->num(iele);
     int numEleDOFs      = allDOFs->num(iele);
     
     elDisp->zero(); 
     elForce->zero();
    
     FullSquareMatrix gradkarray(numEleDOFs,gkelData.get());  
     
     int k;
     for(k=0; k<numEleDOFs; ++k) {
        int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
        if(cn >= 0) {
          (*elDisp)[k] = dsp[cn]; }
        else {
         (*elDisp)[k] = bcx[(*allDOFs)[iele][k]]; }
     }

     dynamic_cast<Element_opt*>(packedEset[iele])->
       gradstiffness(nodes,*gradnodes,gradkarray);

     int i,j;
     for (i=0;i<numEleDOFs;i++) {
       for (j=0;j<numEleDOFs;j++) {
         (*elForce)[i] += gradkarray[i][j]*(*elDisp)[j]  ;
       }
     }

     for(k=0; k<numEleDOFs; ++k) {
       if ( (*allDOFs)[iele][k] == dof ) gradnodalForce += (*elForce)[k];
     }
  }
  return gradnodalForce;

}

//------------------------------------------------------------------------------

void
Domain_opt::getGradduNodalInternalForce(Vector &dsp, Vector &adj, double *bcx, int node, int typ)
{
  // Number of connected elements
  int nconele = nodeToElem->num(node);

  // Extract the relevant DOF
  int dof  = dsa->locate(node, (1 << typ));

  // Loop over all connected elements
  //if(elGrad.get() == 0)  elGrad.reset(new Vector(maxNumDOFs,0.0));
  //else 
  elGrad->zero();
  
  //if(elDisp == 0)  elDisp    = new Vector(maxNumDOFs,0.0);
  //else 
  elDisp->zero();
  

  int icon;
  for(icon=0; icon<nconele; ++icon) {

    int iele = (*nodeToElem)[node][icon];
      
    int NodesPerElement = elemToNode->num(iele);
    int numEleDOFs	= allDOFs->num(iele);
    
    elDisp->zero(); 
    elForce->zero();
    
    FullSquareMatrix karray(numEleDOFs,kelData.get());  
    
    karray=packedEset[iele]->stiffness(nodes,karray.data());

    // ... PERTURB EACH ELEMENTAL DOF

    int i,k;
    for(i=0; i<numEleDOFs; ++i) {
      int cn = c_dsa->getRCN((*allDOFs)[iele][i]);
      if(cn >= 0) {
    	for(k=0; k<numEleDOFs; ++k) {
    	  if ( (*allDOFs)[iele][k] == dof ) adj[cn] += karray[k][i];;
    	} 
      }
    }
  }
}

int
Domain_opt::checkList(int i, int* list, int size)
{
  if (! list) return 1;
 
  int j;
  for(j=0;j<size;j++) if (i==list[j]) return 1;
  return 0;
}


//------------------------------------------------------------------------------
  /****
void
Domain_opt::buildGlobalDK ( DBSparseMatrix *DK )
{
  DK->zeroAll();

  int numActelm = actInf->size();
  
  int iele,actelm,numEleDOFs;
  for(actelm=0; actelm < numActelm; ++actelm) {

    iele       = actInf->getEntry(actelm);
    numEleDOFs = allDOFs->num(iele);

    FullSquareMatrix gradkarray(numEleDOFs,gkelData);

    packedEset[iele]->gradstiffness(nodes,*gradnodes,gradkarray);
    
    DK->add(gradkarray,(*allDOFs)[iele]);
 }
}
  ****/
//---------------------------------------------------------------------------
   /****
double *
Domain_opt::getGradNodalTemperatures()
{
  double *gradnodalTemperatures = new double[numnodes];

  // initialize gradient of nodal temperatures to zero
  int i;
  for(i=0; i<numnodes; ++i)
      gradnodalTemperatures[i] = 0.0;
   
  return gradnodalTemperatures;
}
   ****/
//---------------------------------------------------------------------------
    /****
void
Domain_opt::getTemp(Vector &temp)
{  

  double *nodalTemperatures = 0;
  
  if(sinfo.thermalLoadFlag) nodalTemperatures = getNodalTemperatures();
  if(sinfo.thermoeFlag >=0) nodalTemperatures = temprcvd;

   temp.copy(nodalTemperatures);
   if(sinfo.thermalLoadFlag) delete [] nodalTemperatures;
   
}
    ****/
//---------------------------------------------------------------------------
     /***
void
Domain_opt::getgradTemp(Vector &gtemp)
{
  double *gradnodalTemperatures = 0;
  
  if(sinfo.thermalLoadFlag) gradnodalTemperatures = getGradNodalTemperatures();
  if(sinfo.thermoeFlag >=0) gradnodalTemperatures = gradtemprcvd;

   gtemp.copy(gradnodalTemperatures);
   delete [] gradnodalTemperatures;
}
****/


void Domain_opt::addAnalysis(int sact) 
{  
  if (numAnalysis > 0) 
    {
      analysisList[numAnalysis-1].reset(new SolverInfo(sinfo));
      sinfo = SolverInfo();
    }
  ++numAnalysis;
}


void Domain_opt::updateAnalysis() 
{  
  // make sure that there is at least one analysis  
  if (numAnalysis == 0)  { numAnalysis=1; }
  
  // set solver info for last analysis type  
  analysisList[numAnalysis-1].reset(new SolverInfo(sinfo));
  
  // examine type of analyses for optimization
  if (structoptFlag || reliabilityFlag) 
    {      
    for(int ina=0;ina<numAnalysis; ++ina) 
      {
	activateAnalysis(ina,0);    
	switch ( probType() ) 
	  {    
	  case SolverInfo::Static:
	    numStaticAna++;
	    break;
	    /*
	      case SolverInfo::NonLinStatic:
	      numNlnstcAna++;
	      break;
	      case SolverInfo::Modal:
	      numEigenAna++; 
	      break;
	      case SolverInfo::Dynamic:
	      numDynamAna++; 
	      break;
	      case SolverInfo::NonLinDynam:
	      numNLdynAna++; 
	      break;
	      case SolverInfo::FVibr:
	      numFVibrAna++; 
	      break;
	    */
	  case SolverInfo::Top:
	    /* do nothing */
	    break;
	  default:
	    assert(0);
	    break;
	  }
      }    
    }

  // check consistency of solver parameters
  // 
  // 1. renumbering  
  int renumRef   = analysisList[0]->renum;
  int renumError = 0;  
  for(int ina=1;ina<numAnalysis; ++ina) 
    {
      if(analysisList[ina]->renum != renumRef) { renumError++; }
    }
  if (renumError) 
    {
      fprintf(stderr," *** ERROR: renumber not consistent among analysis\n");
      exit(-1);
    }
  return;
}

void Domain_opt::setUpData()
{
  Domain::setUpData();
  // update list of analysis types
  updateAnalysis();
  return;
}

int Domain_opt::setNeuman(int _numNmLC, BCond* _nbcLC)
{
  // this function implicitly assumes that there is only one load case and
  // one analysis
  numLC = 1;
  numAnaNmLC.reset(new int[1]);  numAnaNmLC[0] = 1;
  mapAnaNmLC.reset(new boost::scoped_array<int>[1]); 
  mapAnaNmLC[0].reset(new int[1]); mapAnaNmLC[0][0] = 0;  
  numNmLC.reset(new int[1]); numNmLC[0] = _numNmLC;   
  nbcLC.reset(new BCond*[1]); nbcLC[0] =  _nbcLC;
  
  return Domain::setNeuman(_numNmLC, _nbcLC);
}


extern "C"      
{
  void  _FORTRAN(dscal)(const int& n, const double& alpha, double* x, const int& incx);
}

//------------------------------------------------------------------------------
template<>
void Domain_opt::scalAdjVec<DComplex>(int len, DComplex* v, DComplex s)
{
  int incx = 2;
  double rscal = ScalarTypes::Real(s);
  double iscal = ScalarTypes::Imag(s);

  _FORTRAN(dscal)(len, rscal, reinterpret_cast<double*>(v),   incx);
  _FORTRAN(dscal)(len, iscal, reinterpret_cast<double*>(v)+1, incx);
}

//------------------------------------------------------------------------------
template<>
void Domain_opt::scalAdjVec<double>(int len, double* v, double s)
{
  int incx = 1;
  _FORTRAN(dscal)(len, s, v, incx);
}

//------------------------------------------------------------------
void Domain_opt::buildNodalDensity(const NodalDensityData& ndd)
{
  if (!nodalDens) { nodalDens.reset(new double[this->numnodes]); }
  switch(ndd.nodDensFunc)
    {  
    case 1:
      {
	// nodalDens = (sum_i rho_i^p )^(1/p)
	for(int inode=0; inode<numnodes; ++inode)
	  {
	    const int nconele = nodeToElem->num(inode);
	    int count = 0;
	    nodalDens[inode] = 0.0;
	    for(int j=0; j<nconele; ++j)
	      {
		const int jele = (*nodeToElem)[inode][j];
		const double dummy = (packedEset[jele]->getProperty() != 0)? packedEset[jele]->getProperty()->rho : 0.0;
		if(dummy > 0.0)
		  {
		    nodalDens[inode] += pow(dummy, ndd.nodDensParam1);
		    ++count;
		  }
	      }
	    if(count) 
	      { nodalDens[inode] = pow(nodalDens[inode],(1/ndd.nodDensParam1)); }
	  }
      }
      break;
    case 2:
      {
	// nodalDen = 1/N*(sum_i rho_i^p )^(1/p)
	// not implemented
	assert(0);
      }
      break;    
    case 3:
      {
	// nodalDens = max(connected nodes)
	for(int inode=0; inode < numnodes; ++inode)
	  {
	    const int nconele = nodeToElem->num(inode);
	    nodalDens[inode] = 0.0;
	    for(int j=0; j<nconele; ++j)
	      {
		const int jele = (*nodeToElem)[inode][j];
		const double dummy = (packedEset[jele]->getProperty() != 0)? packedEset[jele]->getProperty()->rho : 0.0;
		nodalDens[inode] = std::max(nodalDens[inode], dummy);
	      }
	  }
      }
      break;
    default:
      assert(0);
      break;
    }
  return;
}

//--------------------------------------------------------------------------
static void updDensMinMax(double nodalDens, 
			  double densProjParam1, double densProjParam2, double densProjParam3,
			  double& densProj, double& minVal, double& maxVal)
{
  if (nodalDens == 0.0) 
    { densProj = 1.0; }
  else 
    { densProj = 1 + sqrt(densProjParam1)*pow(densProjParam3/nodalDens,densProjParam2); }
  minVal = std::min(minVal, densProj);
  maxVal = std::max(maxVal, densProj);      
  return;
}

//--------------------------------------------------------------------------
void Domain_opt::buildDensProj(const DensityProjData& dpd)
{
  if (dpd.densProjFlag <= 0) 
    { return; }  
  
  if (!densProj) 
    { 
      densProj.reset(new double[this->numUncon()]); 
      std::fill(densProj.get(), densProj.get()+this->numUncon(), 1.0);
    }
  assert(nodalDens != 0);
  
  densProj_minVal    =  DBL_MAX;
  densProj_maxVal    = -DBL_MAX;
  int loc = 0;
    
  for(int inode=0; inode<numnodes; ++inode) 
    {     
      switch (dpd.densProjFlag) 
	{     
	case 1:  
	  //x-translational
	  loc  = c_dsa->locate(inode, DofSet::Xdisp);  
	  if(loc>=0)
	    {
	      updDensMinMax(nodalDens[inode], 
			    dpd.densProjParam1, dpd.densProjParam2, dpd.densProjParam3,
			    densProj[loc], densProj_minVal, densProj_maxVal);
	    }
	  //y-translational
	  loc  = c_dsa->locate(inode, DofSet::Ydisp);  
	  if(loc>=0)
	    {
	      updDensMinMax(nodalDens[inode], 
			    dpd.densProjParam1, dpd.densProjParam2, dpd.densProjParam3,
			    densProj[loc], densProj_minVal, densProj_maxVal);
	    }
	  //z-translational
	  loc  = c_dsa->locate(inode, DofSet::Zdisp);  
	  if(loc>=0)
	    {
	      updDensMinMax(nodalDens[inode], 
			    dpd.densProjParam1, dpd.densProjParam2, dpd.densProjParam3,
			    densProj[loc], densProj_minVal, densProj_maxVal);
	    }
	  //x-rotational
	  loc  = c_dsa->locate(inode, DofSet::Xrot);  
	  if(loc>=0)
	    {
	      updDensMinMax(nodalDens[inode], 
			    dpd.densProjParam1, dpd.densProjParam2, dpd.densProjParam3,
			    densProj[loc], densProj_minVal, densProj_maxVal);
	    }
	  //y-rotational
	  loc  = c_dsa->locate(inode, DofSet::Yrot);  
	  if(loc>=0)
	    {
	      updDensMinMax(nodalDens[inode], 
			    dpd.densProjParam1, dpd.densProjParam2, dpd.densProjParam3,
			    densProj[loc], densProj_minVal, densProj_maxVal);
	    }
	  //z-rotational
	  loc  = c_dsa->locate(inode, DofSet::Zrot);  
	  if(loc>=0)
	    {
	      updDensMinMax(nodalDens[inode], 
			    dpd.densProjParam1, dpd.densProjParam2, dpd.densProjParam3,
			    densProj[loc], densProj_minVal, densProj_maxVal);
	    }
	  break;
	case 2: 
	  // not implemented
	  assert(0);
	  break;
	default:
	  assert(0);
	  break;
    }  
  }
  return;
}


//------------------------------------------------------------------------------
void Domain_opt::scaleDensProj()
{
  assert(densProj != 0 && nodalDens != 0);
  //make it so maximum node has no multiplier
  if (densProj_minVal < DBL_MAX)
    {
      for (int inode=0; inode<this->numUncon(); ++inode) 
	{ densProj[inode] /= densProj_minVal; }
      densProj_maxVal /= densProj_minVal;
      densProj_minVal /= densProj_minVal;
    }
  
  //filePrint(stderr,"\n Maximum value on density projector is %9.5e\n", densProj_maxVal); 
  //filePrint(stderr,"\n Minimum value on density projector is %9.5e\n", densProj_minVal);
}

//------------------------------------------------------------------------------
void Domain_opt::densProjectStiffness(GenFullSquareMatrix<double>& kel, int nel)
{
  // do nothing temporarily
  return;
  //assert(!(packedEset[nel]->isMpcElement() || packedEset[nel]->isPhantomElement()));
  const DensityProjData& dpd = static_cast<GeoSource_opt*>(geoSource)->getDensityProjData();
  if (dpd.densProjFlag > 0)
    {
      const int numEleDOFs = allDOFs->num(nel); 
      for(int jdof=0; jdof<numEleDOFs; ++jdof)
	{
	  const int cn = c_dsa->getRCN((*allDOFs)[nel][jdof]);
	  if (cn < 0) { continue; }
	  const double fac = densProj[cn];
	  for(int idof=0; idof<numEleDOFs; ++idof)
	    {
	      if (c_dsa->getRCN((*allDOFs)[nel][idof])<0) { continue; }
	      kel[idof][jdof] *= fac;
	      kel[jdof][idof] *= fac;
	    }      
	}
    }
  return;
}

//------------------------------------------------------------------------------
void Domain_opt::densProjectStiffnessC(GenFullSquareMatrix<DComplex>& kel, int nel)
{
  // do nothing temporarily
  return;
  const DensityProjData& dpd = static_cast<GeoSource_opt*>(geoSource)->getDensityProjData();
  if (dpd.densProjFlag > 0)
    {
      const int numEleDOFs = allDOFs->num(nel); 
      for(int jdof=0; jdof<numEleDOFs; ++jdof)
	{
	  const int cn = c_dsa->getRCN((*allDOFs)[nel][jdof]);
	  if (cn < 0) { continue; }
	  const double fac = densProj[cn];
	  for(int idof=0; idof<numEleDOFs; ++idof)
	    {
	      if (c_dsa->getRCN((*allDOFs)[nel][idof])<0) { continue; }
	      kel[idof][jdof] *= fac;
	      kel[jdof][idof] *= fac;
	    }      
	}
    }
  return;
}

//------------------------------------------------------------------------------
double Domain_opt::densProjCoeff(int dof)
{ 
  const DensityProjData& dpd = 
    static_cast<GeoSource_opt*>(geoSource)->getDensityProjData();
  return (dpd.densProjFlag > 0)? densProj[dof]: 1.0;
}



#endif
