#include <stdio.h>
#include <Utils.d/dbg_alloca.h>
#include <iostream>

#include <Driver.d/Domain.h>
#include <Math.d/VectorSet.h>
#include <Utils.d/Connectivity.h>
#include <Driver.d/SubDomain.h>
#include <Corotational.d/DistrGeomState.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <Utils.d/pstress.h>
#include <Utils.d/BinFileHandler.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/CornerMaker.h>
#include <Feti.d/Feti.h>
#include <Driver.d/Mpc.h>
#include <Driver.d/DynamProbType.h>
#include <Solvers.d/DiagParallelSolver.h>
#include <Paral.d/MDDynam.h>
#include <Paral.d/GenMS.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>
#include <Paral.d/DomainGroupTask.h>
#ifdef USE_MPI
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#endif
extern Connectivity *procMpcToMpc;

#include <Driver.d/SubDomainFactory.h>

extern const char* problemTypeMessage[];

template<class Scalar>
GenDecDomain<Scalar>::GenDecDomain(Domain *d)
 : mt(d->getTimers())
{
  domain = d;
  initialize(); 
#ifdef USE_MPI
  communicator = new FSCommunicator(structCom);
#else
  communicator = new FSCommunicator;
#endif
}

template<class Scalar>
void 
GenDecDomain<Scalar>::initialize()
{
  subToElem  = 0;
  elemToNode = 0;
  subToNode  = 0;
  nodeToSub  = 0;
  subToSub   = 0;
  mpcToSub_dual = 0;
  mpcToSub_primal = 0;
  mpcToMpc   = 0;
  cpuToSub   = 0;
  subDomain  = 0;
  primalFile = 0;
  stress     = 0;
  weight     = 0;
  globalStress = 0;
  globalWeight = 0;
  grToSub    = 0;
  elemToSub = 0;
  wetInterfaceNodes = 0;
  numWetInterfaceNodes = 0;
  outFreqCount = 0;
  outEigCount = 0;
  localSubToGl = 0; 
  glSubToLocal = 0;
  communicator = 0;
  cpuToCPU = 0; 
  subToCPU = 0; 
  mpcToCpu = 0;
  numPrimalMpc = 0;
  numDualMpc = 0;
  firstOutput = true;
  soweredInput = false;
  nodeVecInfo = 0;
} 

template<class Scalar>
GenDecDomain<Scalar>::~GenDecDomain()
{
  if(elemToNode) { delete elemToNode; elemToNode = 0; }
  if(subToNode) { delete subToNode; subToNode = 0; }
  if(nodeToSub) { delete nodeToSub; nodeToSub = 0; }
  if(subToSub) { delete subToSub; subToSub = 0; }
  if(mpcToSub_dual) { delete mpcToSub_dual; mpcToSub_dual = 0; }
  if(mpcToSub_primal) { delete mpcToSub_primal; mpcToSub_primal = 0; }
  if(mpcToMpc) { delete mpcToMpc; mpcToMpc = 0; }
  if(subDomain) {
    for(int i=0; i<numSub; ++i)
      if(subDomain[i]) { delete subDomain[i]; subDomain[i] = 0; }
     delete [] subDomain; subDomain = 0;
  }  
  if(primalFile) { delete primalFile; primalFile = 0; }
  if(stress) { delete stress; stress = 0; }
  if(weight) { delete weight; weight = 0; }
  if(globalStress) { delete [] globalStress; globalStress = 0;}
  if(globalWeight) {delete [] globalWeight; globalWeight = 0;}
  if(grToSub) { delete grToSub; grToSub = 0; }

  if(elemToSub) { delete elemToSub; elemToSub = 0; }
  if(wetInterfaceNodes) { delete [] wetInterfaceNodes; }
  if(!cpuToSub && localSubToGl) { delete [] localSubToGl; localSubToGl = 0; }
  if(communicator) { delete communicator; communicator = 0; }
  if(cpuToCPU) { delete cpuToCPU; cpuToCPU = 0; }
  if(subToCPU) { delete [] subToCPU; subToCPU = 0; }
  if(glSubToLocal) { delete [] glSubToLocal; glSubToLocal = 0; }
  if(mpcToCpu) { delete mpcToCpu; mpcToCpu = 0; }
  if(subToElem) { delete subToElem; subToElem = 0; }
  if(nodeVecInfo) delete nodeVecInfo;
}

template<class Scalar>
void GenDecDomain<Scalar>::addBMPCs()
{
  // step 1. create all the candidate bmpcs 
  vector<LMPCons *> bmpcs; // candidate bmpcs
  int nbmpc = 0; // number of bmpcs
  int nsn; // number of subdomains sharing node
  for(int i=0; i<nodeToSub->csize(); ++i) {
    if((nsn = nodeToSub->num(i)) < 2) continue; 
    for(int dof=0; dof<DofSet::max_known_nonL_dof; ++dof) { // x,y,z translations and rotations, temperature, helmholtz and intpress
      for(int j=0; j<nsn; ++j) {
        for(int k=j+1; k<nsn; ++k) {
          LMPCons *bmpc = new LMPCons(nbmpc, 0.0, new LMPCTerm(i,dof,0.0)); // dummy mpc, set coef later in distributeMPCs
          bmpc->type = 2;
          int subA = (*nodeToSub)[i][j];
          int subB = (*nodeToSub)[i][k];
          bmpc->psub = (subA > subB) ? subA : subB;
          bmpc->nsub = (subA > subB) ? subB : subA;
          bmpcs.push_back(bmpc);
          nbmpc++;
        }
      }
    }
  }
  //cerr << "number of candidate bmpcs = " << nbmpc << endl;

  // step 2. check which of the candidates are active (ie both subdomains share an active "r" dof)
  int *pstatus = new int[2*nbmpc];
  int *nstatus = pstatus+nbmpc;
  for(int i=0; i<nbmpc; ++i) pstatus[i] = nstatus[i] = 0;
  paralApply(numSub, subDomain, &GenSubDomain<Scalar>::bmpcQualify, &bmpcs, pstatus, nstatus);
#ifdef DISTRIBUTED
  communicator->globalSum(2*nbmpc, pstatus);
#endif

  // step 3. now add the active bmpcs to the global mpc collection
  // (check for and eliminate redundant bmpcs)
  nbmpc = 0; // reset counter
  int total = 0;
  for(int i=0; i<nodeToSub->csize(); ++i) {
    if((nsn = nodeToSub->num(i)) < 2) continue;
    for(int dof=0; dof<DofSet::max_known_nonL_dof; ++dof) { // x,y,z translations and rotations, temperature, helmholtz and intpress
      bool got_enough = false;
      for(int j=0; j<nsn; ++j) {
        int count = 0;
        for(int k=j+1; k<nsn; ++k) {
          //if(pstatus[nbmpc] && nstatus[nbmpc] && (*nodeToSub)[j][0] >= (*nodeToSub)[k][0]) {
          if(pstatus[nbmpc] && nstatus[nbmpc] && !got_enough) {
            bmpcs[nbmpc]->lmpcnum = 0; // shouldn't be used 
            domain->addLMPC(bmpcs[nbmpc],false);
            count++; total++;
          }
          nbmpc++;
        }
        if(count > 0) { 
          got_enough = true;
          //cerr << "added " << count << " bmpcs at node " << i << " dof " << dof << endl;
        }
      }
    }
  }
  //cerr << "final number of bmpcs = " << total << endl;
}

template<class Scalar>
void GenDecDomain<Scalar>::getSharedNodes()
{
  // Start timer
  startTimerMemory(mt.makeConnectivity, mt.memoryConnect);

  // PJSA 9-1-04 for coupled_dph all the wet interface nodes need to be included in the sharedNodes list in SComm
  bool coupled_dph = false;
  int *wetInterfaceNodeMap = 0;
  if(domain->solInfo().isCoupled) {
    coupled_dph = true;
    //wetInterfaceNodes = domain->getAllWetInterfaceNodes(numWetInterfaceNodes);
    wetInterfaceNodeMap = new int[nodeToSub->csize()];
    for(int i=0; i<nodeToSub->csize(); ++i) wetInterfaceNodeMap[i] = -1;
    for(int i=0; i<numWetInterfaceNodes; ++i) wetInterfaceNodeMap[wetInterfaceNodes[i]] = i;
  }

  // create each subdomain's interface lists
  int iSub, jSub, subJ, iNode;
  int totConnect;
  int *nConnect = new int[subToNode->csize()];
  int *flag = new int[subToNode->csize()];
  for(iSub = 0; iSub < subToNode->csize(); ++iSub) {
     flag[iSub] = -1;
     nConnect[iSub] = 0;
  }

  // Count connectivity
  totConnect = 0;
  for(iSub = 0; iSub < subToNode->csize(); ++iSub) {
    for(iNode = 0; iNode < subToNode->num(iSub); ++iNode) { // loop over the nodes
      int thisNode = (*subToNode)[iSub][iNode];
      bool isWetInterfaceNode = (coupled_dph && (wetInterfaceNodeMap[thisNode] != -1)) ? true : false;
      for(jSub = 0; jSub < nodeToSub->num(thisNode); ++jSub) {
        // loop over the subdomains connected to this node
        subJ = (*nodeToSub)[thisNode][jSub];
        bool neighbWithSelf = ((subJ == iSub) && isWetInterfaceNode && (domain->solInfo().fetiInfo.fsi_corner == 0)) ? true : false;
        if((subJ > iSub) || neighbWithSelf) {
          // only deal with connection to highered numbered subdomains, guarantees symmetry of lists
          if(flag[subJ] != iSub) {
            flag[subJ] = iSub;
            nConnect[iSub]++;
            if(!neighbWithSelf) {
              nConnect[subJ]++;
              totConnect += 2;
            }
            else totConnect += 1;
          }
        }
      }
    }
  }

  int **nodeCount = new int*[subToNode->csize()];
  int **connectedDomain = new int*[subToNode->csize()];
  int **remoteID = new int*[subToNode->csize()];

  // Allocate memory for list of connected subdomains
  // (de-allocated in ~SComm)
  for(iSub = 0; iSub < subToNode->csize(); ++iSub) {
    int size = nConnect[iSub];
    connectedDomain[iSub] = new int[size];
    remoteID[iSub] = new int[size];
    nodeCount[iSub] = new int[size];
    flag[iSub] = -1;
    nConnect[iSub] = 0;
  }

  int *whichLocal  = new int[subToNode->csize()];
  int *whichRemote = new int[subToNode->csize()];

  for(iSub=0; iSub < subToNode->csize(); ++iSub) {
    for(iNode = 0; iNode < subToNode->num(iSub); ++iNode) {
      int nd = (*subToNode)[iSub][iNode];
      bool isWetInterfaceNode = (coupled_dph && (wetInterfaceNodeMap[nd] != -1)) ? true : false;
      for(jSub = 0; jSub < nodeToSub->num(nd); ++jSub) {
        subJ = (*nodeToSub)[nd][jSub];
        bool neighbWithSelf = ((subJ == iSub) && isWetInterfaceNode && (domain->solInfo().fetiInfo.fsi_corner == 0)) ? true : false;
        if((subJ > iSub) || neighbWithSelf) {
          if(flag[subJ] != iSub) { // attribute location for this sub
            flag[subJ] = iSub;
            connectedDomain[subJ][nConnect[subJ]] = iSub;
            connectedDomain[iSub][nConnect[iSub]] = subJ;
            remoteID[subJ][nConnect[subJ]] = nConnect[iSub];
            remoteID[iSub][nConnect[iSub]] = nConnect[subJ];
            if(!neighbWithSelf) whichLocal[subJ] = nConnect[iSub]++;
            else whichLocal[subJ] = nConnect[iSub];
            whichRemote[subJ] = nConnect[subJ]++;
            nodeCount[iSub][whichLocal[subJ]] = 1;
            nodeCount[subJ][whichRemote[subJ]] = 1;
          }
          else {
            nodeCount[iSub][whichLocal[subJ]]++;
            if(!neighbWithSelf) nodeCount[subJ][whichRemote[subJ]]++;
          }
        }
      }
    }
  }

  // allocate memory for interface node lists
  Connectivity **interfNode = new Connectivity *[subToNode->csize()];
  for(iSub=0; iSub < subToNode->csize(); ++iSub)
    interfNode[iSub] = new Connectivity(nConnect[iSub], nodeCount[iSub]);

  // fill the lists
  for(iSub = 0; iSub < subToNode->csize(); ++iSub) {
     flag[iSub]     = -1;
     nConnect[iSub] =  0;
  }

  for(iSub = 0; iSub < subToNode->csize(); ++iSub) {
    for(iNode = 0; iNode < subToNode->num(iSub); ++iNode) {
      int nd = (*subToNode)[iSub][iNode];
      bool isWetInterfaceNode = (coupled_dph && (wetInterfaceNodeMap[nd] != -1)) ? true : false;
      for(jSub = 0; jSub < nodeToSub->num(nd); ++jSub) {
        subJ = (*nodeToSub)[nd][jSub];
        bool neighbWithSelf = ((subJ == iSub) && isWetInterfaceNode && (domain->solInfo().fetiInfo.fsi_corner == 0)) ? true : false;
        if((subJ > iSub) || neighbWithSelf) {
          if(flag[subJ] != iSub) { // attribute location for this sub
            flag[subJ] = iSub;
            if(!neighbWithSelf) whichLocal[subJ] = nConnect[iSub]++;
            else whichLocal[subJ] = nConnect[iSub];
            whichRemote[subJ] = nConnect[subJ]++;
            (*interfNode[iSub])[whichLocal[subJ]][0] = nd;
            (*interfNode[subJ])[whichRemote[subJ]][0] = nd;
            nodeCount[iSub][whichLocal[subJ]]=1;
            nodeCount[subJ][whichRemote[subJ]]=1;
          }
          else {
            int il = nodeCount[iSub][whichLocal[subJ]]++;
            (*interfNode[iSub])[whichLocal[subJ]][il] = nd;
            if(!neighbWithSelf) {
              int jl = nodeCount[subJ][whichRemote[subJ]]++;
              (*interfNode[subJ])[whichRemote[subJ]][jl] = nd;
            }
          }
        }
      }
    }
  }
  delete [] flag;
  delete [] whichLocal;
  delete [] whichRemote;
  for(iSub = 0; iSub < subToNode->csize(); ++iSub) delete [] nodeCount[iSub];
  delete [] nodeCount;
  delete [] wetInterfaceNodeMap;

  for(iSub = 0; iSub < numSub; ++iSub) {
    int subI = (localSubToGl) ? localSubToGl[iSub] : iSub;
    GenSubDomain<Scalar> **subds = new GenSubDomain<Scalar> * [nConnect[subI]];
    for(jSub = 0; jSub < nConnect[subI]; ++jSub) {
       int subJ = glSubToLocal[connectedDomain[subI][jSub]];
       subds[jSub] = (subJ >= 0) ? subDomain[subJ] : NULL;
    }
    SComm *sc = new SComm(nConnect[subI], connectedDomain[subI], subds,
                          remoteID[subI], interfNode[subI]);
    sc->locSubNum = iSub;
    sc->glSubToLocal = glSubToLocal;
    subDomain[iSub]->setSComm(sc);
    delete [] subds;
  }

  delete [] remoteID;
  delete [] connectedDomain;
  delete [] nConnect;
  delete [] interfNode;

  paralApply(numSub, subDomain, &GenSubDomain<Scalar>::renumberSharedNodes);

  stopTimerMemory(mt.makeConnectivity, mt.memoryConnect);
}

template<class Scalar>
void
GenDecDomain<Scalar>::preProcessBCsEtc()
{
  if(!soweredInput) {
    distributeBCs();
    distributeControlLawData();
    distributeDiscreteMass();
    paralApply(numSub, subDomain, &GenSubDomain<Scalar>::renumberBCsEtc); 
    paralApply(numSub, subDomain, &GenSubDomain<Scalar>::renumberControlLaw); 
  }
  else {
    if (geoSource->binaryInputControlLeft) {
      distributeControlLawData();
      paralApply(numSub, subDomain, &GenSubDomain<Scalar>::renumberControlLaw);
    }
  }

  paralApply(numSub, subDomain, &BaseSub::makeCDSA); 
  paralApply(numSub, subDomain, &Domain::makeNodeToNode_sommer);
}

template<class Scalar>
void
GenDecDomain<Scalar>::preProcessMPCs()
{
#ifdef SOWER_SURFS
  if(soweredInput) {
    //HB compute mortar LMPCs
    domain->SetMortarPairing();
    domain->SetUpSurfaces(); //domain->SetUpSurfaces(geoSource->sower_nodes);
    if(domain->solInfo().newmarkBeta != 0.0) { // not for explicit dynamics
      domain->ComputeMortarLMPC();
      domain->computeMatchingWetInterfaceLMPC();
      domain->CreateMortarToMPC();
    }
#ifdef MORTAR_DEBUG
    domain->PrintSurfaceEntities();
    domain->PrintMortarConds();
#endif
  }
#endif
  if(domain->solInfo().fetiInfo.bmpc) addBMPCs();
  if(domain->getNumLMPC() > 0) {
    if(verboseFlag) filePrint(stderr, " ... Applying the Multi-Point Constraints");
    // check for mpcs involving bad nodes and constrained DOFs
    if(nodeToSub) domain->checkLMPCs(nodeToSub);
    // select which mpcs are to be included in the coarse problem
    domain->setPrimalLMPCs(numDualMpc, numPrimalMpc);
    // distribute mpcs
    execParal(numSub, this, &GenDecDomain<Scalar>::extractSubDomainMPCs);
    makeMpcToSub(); // new version works for distributed data
    // renumber local mpcs
    paralApplyToAll(numSub, subDomain, &GenSubDomain<Scalar>::renumberMPCs);
    // locate and store dsa, c_dsa and cc_dsa numbering for future reference
    paralApply(numSub, subDomain, &GenSubDomain<Scalar>::locateMpcDofs);
    // make mpcToMpc connectivity
    makeMpcToMpc(); // this should be in the Feti-DP constructor
    if(verboseFlag) filePrint(stderr, " ...\n");
  }
  else domain->solInfo().getFetiInfo().mpcflag = 0;
}

template<class Scalar>
void
GenDecDomain<Scalar>::deleteMPCs()
{
  paralApply(numSub, subDomain, &GenSubDomain<Scalar>::deleteMPCs);
  if(mpcToSub_dual) delete mpcToSub_dual; mpcToSub_dual = 0;
  if(mpcToMpc) delete mpcToMpc; mpcToMpc = 0;
  if(mpcToCpu) delete mpcToCpu; mpcToCpu = 0;
  numDualMpc = 0;
}

template<class Scalar>
void
GenDecDomain<Scalar>::reProcessMPCs()
{
  deleteMPCs();
  preProcessMPCs();
  getSharedMPCs();
  paralApply(numSub, subDomain, &BaseSub::mergeInterfaces);
  paralApply(numSub, subDomain, &GenSubDomain<Scalar>::applyMpcSplitting);
}

template<class Scalar>
void
GenDecDomain<Scalar>::extractSubDomainMPCs(int iSub)
{
  if(numDualMpc) subDomain[iSub]->extractMPCs(domain->numLMPC, domain->lmpc);
  if(numPrimalMpc) subDomain[iSub]->extractMPCs_primal(domain->numLMPC, domain->lmpc);
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeSubToSubEtc()
{
  if(soweredInput) {
    subToSub = geoSource->getSubToSub();
    subToNode = geoSource->getSubToNode();
    subToNode->sortTargets(); // PJSA 11-16-2006

    mt.memoryNodeToSub -= memoryUsed();
    nodeToSub = subToNode->reverse();
    mt.memoryNodeToSub += memoryUsed();

    geoSource->setNumNodes(nodeToSub->csize());
    geoSource->computeClusterInfo(localSubToGl[0]);
  }
  else {
    mt.memoryElemToNode -= memoryUsed();
    elemToNode = new Connectivity(&domain->packedEset);
    mt.memoryElemToNode += memoryUsed();

    mt.memorySubToNode -= memoryUsed();
    subToNode = subToElem->transcon(elemToNode);
    subToNode->sortTargets(); // PJSA 11-16-2006
    mt.memorySubToNode += memoryUsed();

    mt.memoryNodeToSub -= memoryUsed();
    nodeToSub = subToNode->reverse();
    mt.memoryNodeToSub += memoryUsed();
    domain->setNumNodes(nodeToSub->csize());

    mt.memorySubToNode -= memoryUsed();
    subToSub = subToNode->transcon(nodeToSub);
    mt.memorySubToNode += memoryUsed();

    //if(domain->solInfo().isAcoustic()) {
    if(domain->numSSN() || domain->solInfo().isCoupled) { // sommerfeld, scatter, wet, distributed neum
      mt.memoryNodeToElem -= memoryUsed();
      domain->nodeToElem = elemToNode->reverse();
      mt.memoryNodeToElem += memoryUsed();

      mt.memoryElemToSub -= memoryUsed();
      elemToSub = subToElem->reverse();
      mt.memoryElemToSub += memoryUsed();
      if(domain->numSSN() > 0) domain->checkSommerTypeBC(domain, elemToNode, domain->nodeToElem); // flip normals if necessary
    }
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeSubDomains() 
{
  makeSubDMaps();
  makeSubToSubEtc();
  subDomain = new GenSubDomain<Scalar> *[numSub];

  startTimerMemory(mt.makeSubDomains, mt.memorySubdomain);
  if(soweredInput) {
    for(int iSub = 0; iSub < this->numSub; iSub++) { 
      subDomain[iSub] = geoSource->template readDistributedInputFiles<Scalar>(iSub, localSubToGl[iSub]);  
    }
#ifdef SOWER_SURFS
    geoSource->template readDistributedSurfs<Scalar>(localSubToGl[0]); //pass dummy sub number
#endif
  }
  else {
    execParal(numSub, this, &GenDecDomain<Scalar>::constructSubDomains);
    if(domain->solInfo().isCoupled) { // JLchange
      for(int iSub = 0; iSub < numSub; ++iSub) subDomain[iSub]->setnodeToSubConnectivity(nodeToSub); // JLchange
      addFsiElements();
    }
    paralApply(numSub, subDomain, &GenSubDomain<Scalar>::renumberElements); 
  }
  
  paralApply(numSub, subDomain, &BaseSub::makeDSA); 
  stopTimerMemory(mt.makeSubDomains, mt.memorySubdomain);
/*
  //added JF to check the memory allocated
  fprintf(stderr,"in Driver.d/DecDomain.C\n");
  for(int iSub = 0; iSub < this->numSub;iSub++) {
    fprintf(stderr,"SubDomain # = %d\n",iSub);
    fprintf(stderr,"   Number of nodes in the nodeSet = %d\n",subDomain[iSub]->numNodes());
    fprintf(stderr,"   Number of elements in the elemSet = %d\n",subDomain[iSub]->numElements());
    fprintf(stderr,"  size of the memory allocated:\n");
    fprintf(stderr,"   name of geoSource pointer: %d\n",geoSource);
    fprintf(stderr,"   nodes = %d\n",geoSource->numNode());
    fprintf(stderr,"   elements = %d\n",geoSource->numElem());
  }
*/
}

template<class Scalar>
void
GenDecDomain<Scalar>::distributeDiscreteMass()
{
  // Distribute masses to the subdomains
  DMassData *cmass = domain->firstDiMass;
  while(cmass != 0) {
   int node = cmass->node;
   for(int iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
     int tSub = glSubToLocal[(*nodeToSub)[node][iSub]];
     if(tSub >= 0)
        subDomain[tSub]->addDMass(node, cmass->dof, cmass->diMass);
   }
   cmass = cmass->next;
  }
}

template<class Scalar>
GenFetiSolver<Scalar> *
GenDecDomain<Scalar>::getFetiSolver()
{
 FetiInfo *finfo = &domain->solInfo().getFetiInfo();
 if(finfo->version == FetiInfo::fetidp) {
   GenSolver<Scalar> **sysMatrices = 0;
   GenSparseMatrix<Scalar> **sysMat = 0;
   Rbm **rbms = 0;
   bool computeRbms = (domain->solInfo().isStatic() || domain->probType() == SolverInfo::Modal) && !geoSource->isShifted();
   return new GenFetiDPSolver<Scalar>(numSub, subDomain, subToSub, finfo,
                                      communicator, glSubToLocal, mpcToSub_dual, mpcToSub_primal,
                                      mpcToMpc, mpcToCpu, cpuToSub, grToSub, sysMatrices, sysMat, rbms, 0, computeRbms);
 }
 else {
   return new GenFetiSolver<Scalar>(numSub, subDomain, subToSub, finfo, communicator,
                                    glSubToLocal, mpcToSub_dual, cpuToSub);
 }
}
template<class Scalar>
GenFetiSolver<Scalar> *
GenDecDomain<Scalar>::getDynamicFetiSolver(GenDomainGroupTask<Scalar> &dgt)
{
 FetiInfo *finfo = &domain->solInfo().getFetiInfo();
 bool computeRbms = (domain->solInfo().isStatic() || domain->probType() == SolverInfo::Modal) && !geoSource->isShifted(); 
 if(finfo->version == FetiInfo::fetidp) {
   return new GenFetiDPSolver<Scalar>(numSub, subDomain, subToSub, finfo, communicator, glSubToLocal,
                                      mpcToSub_dual, mpcToSub_primal, mpcToMpc, mpcToCpu, cpuToSub, grToSub,
                                      dgt.dynMats, dgt.spMats, dgt.rbms, 0, computeRbms);
 }
 else {
   return new GenFetiSolver<Scalar>(numSub, subDomain, subToSub, finfo, communicator,
                                    glSubToLocal, mpcToSub_dual, cpuToSub,
                                    dgt.dynMats, dgt.spMats, dgt.rbms);
 }
}

template<class Scalar>
DiagParallelSolver<Scalar> *
GenDecDomain<Scalar>::getDiagSolver(int nSub, GenSubDomain<Scalar> **sd,
              GenSolver<Scalar> **sol)//GenDomainGroupTask<Scalar> &dgt)
{
 return new DiagParallelSolver<Scalar>(nSub, sd, sol, cpuToSub, communicator);
}

template<class Scalar>
void
GenDecDomain<Scalar>::getCPUMap()
{
  mt.memoryCPUMAP -= memoryUsed();
#ifdef DISTRIBUTED
  char *mapName = geoSource->getCpuMapFile(); // PJSA
  FILE *f = fopen(mapName,"r");
  numCPU = geoSource->getCPUMap(f, globalNumSub);
  subToCPU = geoSource->getSubToCPU();
  cpuToCPU = geoSource->getCpuTOCPU();
  if(f) fclose(f);
#else
  numCPU = 1; myCPU = 0;
  geoSource->createSingleCpuToSub(globalNumSub);
#endif
  cpuToSub = geoSource->getCpuToSub();
  mt.memoryCPUMAP += memoryUsed();
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeSubDMaps()
{
  glSubToLocal = new int[globalNumSub];
  if(cpuToSub) {
    numSub = cpuToSub->num(myCPU);
    for(int iSub = 0; iSub < globalNumSub; ++iSub)
      glSubToLocal[iSub] = -1;
    for(int iSub = 0; iSub < numSub; ++iSub) 
      glSubToLocal[ (*cpuToSub)[myCPU][iSub] ] = iSub;
    localSubToGl = (*cpuToSub)[myCPU];
  }
  else { // shared memory, all subds on single cpu
    for(int i=0; i<globalNumSub; ++i) glSubToLocal[i] = i;
    numSub = globalNumSub;
    localSubToGl = new int[numSub];
    for(int i=0; i<numSub; ++i) localSubToGl[i] = i;
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::preProcess()
{
 filePrint(stderr, "%s", problemTypeMessage[domain->probType()]); 
 soweredInput = geoSource->binaryInput;

 if(verboseFlag) filePrint(stderr, " ... Reading Decomposition File     ...\n");
 subToElem = geoSource->getDecomposition();
 subToElem->sortTargets(); // PJSA 11-16-2006
 globalNumSub = subToElem->csize();

 getCPUMap();

 if(verboseFlag) filePrint(stderr, " ... Making the Subdomains          ...\n");
 makeSubDomains();

 preProcessBCsEtc();

 preProcessFSIs();// FLuid-Structure Interaction

 getSharedNodes();

 makeCorners();// Corners for FETI-DP
 
 getSharedDOFs();

 preProcessMPCs();//Multi-Point Constraint

 getSharedFSIs();

 getSharedMPCs();

 paralApply(numSub, subDomain, &BaseSub::mergeInterfaces);
 paralApply(numSub, subDomain, &GenSubDomain<Scalar>::applySplitting);
 //paralApply(numSub, subDomain, &GenSubDomain<Scalar>::initSrc);

 makeInternalInfo();
 makeInternalInfo2();
 makeNodeInfo();

#ifdef DISTRIBUTED
 geoSource->setNumNodalOutput();
 if(geoSource->getNumNodalOutput()) {
   for(int i=0; i<numSub; ++i)
     geoSource->distributeOutputNodesX(subDomain[i], nodeToSub);
 }
#endif

 // free up some memory
 //delete nodeToSub; nodeToSub = 0;
 delete elemToSub; elemToSub = 0;
 if(!geoSource->elemOutput() && elemToNode) { delete elemToNode; elemToNode = 0; }
}

template<class Scalar>
void
GenDecDomain<Scalar>::scaleDisp(GenDistrVector<Scalar> &u)
{
  execParal(numSub, this, &GenDecDomain<Scalar>::scaleSubDisp, u);
}

template<class Scalar>
void
GenDecDomain<Scalar>::scaleSubDisp(int iSub, GenDistrVector<Scalar> &u)
{
  subDomain[iSub]->scaleDisp(u.subData(iSub));
}

template<class Scalar>
void
GenDecDomain<Scalar>::postProcessing(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &f,
                                     double eigV, GenDistrVector<Scalar> *aeroF, int x,
                                     GenMDDynamMat<Scalar> *dynOps, SysState<GenDistrVector<Scalar> > *distState, int ndflag)
{
  int numOutInfo = geoSource->getNumOutInfo();
  if(numOutInfo == 0) return;

  if(verboseFlag && numOutInfo && x == 0 && ndflag == 0 && !domain->solInfo().isDynam())
    filePrint(stderr," ... Postprocessing                 ...\n");

  Scalar *globVal = 0;  
  int numNodes = geoSource->numNode(); 
  int i, j, iSub, inode;

  // initialize and merge displacements from subdomains into global array
  Scalar (*mergedDis)[11] = new Scalar[numNodes][11];//DofSet::max_known_nonL_dof
  for(i = 0; i < numNodes; ++i)
    for(j=0; j<11; j++) mergedDis[i][j] = 0.0;//DofSet::max_known_nonL_dof
  for(iSub = 0; iSub < numSub; ++iSub)
    subDomain[iSub]->mergeAllDisp(mergedDis, u.subData(iSub));

  // intialize and merge aeroelastic forces from subdomains into global array
  Scalar (*mergedAeroF)[6] = 0;
  if(aeroF) {
    mergedAeroF = new Scalar[numNodes][6];
    for(i = 0; i < numNodes; ++i)
      for(j=0; j<6; ++j) mergedAeroF[i][j] = 0.0;
    for(iSub = 0; iSub < numSub; ++iSub) 
      subDomain[iSub]->mergeForces(mergedAeroF, aeroF->subData(iSub));
  }

  // initialize and merge velocities and accelerations from subdomains into global array
  GenDistrVector<Scalar> *v_n = 0, *a_n = 0;
  Scalar (*mergedVel)[11] = 0, (*mergedAcc)[11] = 0;//DofSet::max_known_nonL_dof
  if(distState) {
    v_n = &distState->getVeloc();
    a_n = &distState->getAccel();
    mergedVel = new Scalar[numNodes][11];//DofSet::max_known_nonL_dof
    mergedAcc = new Scalar[numNodes][11];//DofSet::max_known_nonL_dof
    for(i = 0; i < numNodes; ++i)
      for(j=0; j<11; ++j) mergedVel[i][j] = mergedAcc[i][j] = 0.0; //DofSet::max_known_nonL_dof
    for(iSub = 0; iSub < numSub; ++iSub) {
      subDomain[iSub]->mergeAllVeloc(mergedVel, v_n->subData(iSub));
      subDomain[iSub]->mergeAllAccel(mergedAcc, a_n->subData(iSub));
    }
  }

  // compute current time (or frequency in the case of a helmholtz problem)
  double time;
  if(geoSource->isShifted() && domain->probType() != SolverInfo::Modal) {
    time = domain->getFrequencyOrWavenumber();
    if(domain->solInfo().doFreqSweep) x = outFreqCount++;
  } else if(domain->probType() == SolverInfo::Modal) {
    time = eigV;
    if(domain->solInfo().doEigSweep) x = outEigCount++; 
  }
  else time = x*domain->solInfo().getTimeStep();

  // get output information
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // open output files
  if(x == 0 && firstOutput) geoSource->openOutputFiles();


  Scalar Wext = 0, Waero = 0, Wela = 0, Wkin = 0, Wdmp = 0;

  for(i = 0; i < numOutInfo; i++) { 
    if(oinfo[i].ndtype != ndflag) continue;
    if(ndflag !=0 && oinfo[i].type != OutputInfo::Disp6DOF && oinfo[i].type !=  OutputInfo::Displacement) continue;
    if(oinfo[i].interval != 0 && x % oinfo[i].interval == 0) {
      switch(oinfo[i].type) {
        case OutputInfo::EigenPair:
        case OutputInfo::FreqRespModes:
        case OutputInfo::Displacement:
          getPrimalVector(i, mergedDis, numNodes, 3, time);
          break;
        case OutputInfo::Velocity:
          if(distState) getPrimalVector(i, mergedVel, numNodes, 3, time);
          break;
        case OutputInfo::Acceleration:
          if(distState) getPrimalVector(i, mergedAcc, numNodes, 3, time);
          break;
        case OutputInfo::Disp6DOF:
          getPrimalVector(i, mergedDis, numNodes, 6, time);
          break;
        case OutputInfo::Velocity6:
          if(distState) getPrimalVector(i, mergedVel, numNodes, 6, time);
          break;
        case OutputInfo::Accel6:
          if(distState) getPrimalVector(i, mergedAcc, numNodes, 6, time);
          break;
        case OutputInfo::Temperature:
          getPrimalScalar(i, mergedDis, numNodes, 6, time);
          break;
        case OutputInfo::TemperatureFirstTimeDerivative:
          if(distState) getPrimalScalar(i, mergedVel, numNodes, 6, time);
          break;
        case OutputInfo::EigenPressure:
        case OutputInfo::HelmholtzModes:
        case OutputInfo::Helmholtz:
        case OutputInfo::EigenSlosh:
          getPrimalScalar(i, mergedDis, numNodes, 7, time);
          break;
        case OutputInfo::PressureFirstTimeDerivative:
          if(distState) getPrimalScalar(i, mergedVel, numNodes, 7, time);
          break;
        case OutputInfo::PressureSecondTimeDerivative:
          if(distState) getPrimalScalar(i, mergedAcc, numNodes, 7, time);
          break;
        case OutputInfo::StressXX:
          getStressStrain(u, i, SXX, time);
          break;
        case OutputInfo::StressYY:
          getStressStrain(u, i, SYY, time);
          break;
        case OutputInfo::StressZZ:
          getStressStrain(u, i, SZZ, time);
          break;
        case OutputInfo::StressXY:
          getStressStrain(u, i, SXY, time);
          break; 
        case OutputInfo::StressYZ:
          getStressStrain(u, i, SYZ, time);
          break;
        case OutputInfo::StressXZ:
          getStressStrain(u, i, SXZ, time);
          break;
        case OutputInfo::StrainXX:
          getStressStrain(u, i, EXX, time);
          break;
        case OutputInfo::StrainYY:
          getStressStrain(u, i, EYY, time);
          break;
        case OutputInfo::StrainZZ:
          getStressStrain(u, i, EZZ, time);
          break;
        case OutputInfo::StrainXY:
          getStressStrain(u, i, EXY, time);
          break;
        case OutputInfo::StrainYZ:
          getStressStrain(u, i, EYZ, time);
          break;
        case OutputInfo::StrainXZ:
          getStressStrain(u, i, EXZ, time);
          break;
        case OutputInfo::StressVM:
          getStressStrain(u, i, VON, time);
          break;
        case OutputInfo::StrainVM:
          getStressStrain(u,i,STRAINVON, time);
          break;
        case OutputInfo::ContactPressure:
          getStressStrain(u, i, CONPRESS, time);
          break;
        case OutputInfo::Damage:
          getStressStrain(u, i, DAMAGE, time);
          break;
        case OutputInfo::StressPR1:
          getPrincipalStress(u, i, PSTRESS1, time);
          break;
        case OutputInfo::StressPR2:
          getPrincipalStress(u, i, PSTRESS2, time);
          break;
        case OutputInfo::StressPR3:
          getPrincipalStress(u, i, PSTRESS3, time);
          break;
        case OutputInfo::StrainPR1:
          getPrincipalStress(u, i, PSTRAIN1, time);
          break;
        case OutputInfo::StrainPR2:
          getPrincipalStress(u, i, PSTRAIN2, time);
          break;
        case OutputInfo::StrainPR3:
          getPrincipalStress(u, i, PSTRAIN3, time);
          break;
        case OutputInfo::StressPR1Direc:
          getPrincipalStress(u, i, PSTRESS1DIREC, time);
          break;
        case OutputInfo::StressPR2Direc:
          getPrincipalStress(u, i, PSTRESS2DIREC, time);
          break;
        case OutputInfo::StressPR3Direc:
          getPrincipalStress(u, i, PSTRESS3DIREC, time);
          break;
        case OutputInfo::StrainPR1Direc:
          getPrincipalStress(u, i, PSTRAIN1DIREC, time);
          break;
        case OutputInfo::StrainPR2Direc:
          getPrincipalStress(u, i, PSTRAIN2DIREC, time);
          break;
        case OutputInfo::StrainPR3Direc:
          getPrincipalStress(u, i, PSTRAIN3DIREC, time);
          break;
        case OutputInfo::InXForce:
          getElementForce(u, i, INX, time);
          break;
        case OutputInfo::InYForce:
          getElementForce(u, i, INY, time);
          break;
        case OutputInfo::InZForce:
          getElementForce(u, i, INZ, time);
          break;
        case OutputInfo::AXMoment:
          getElementForce(u, i, AXM, time);
          break;
        case OutputInfo::AYMoment:
          getElementForce(u, i, AYM, time);
          break;
        case OutputInfo::AZMoment:
          getElementForce(u, i, AZM, time);
          break;
        case OutputInfo::DispX:
          getPrimalScalar(i, mergedDis, numNodes, 0, time);
          break;
        case OutputInfo::DispY:
          getPrimalScalar(i, mergedDis, numNodes, 1, time);
          break;
        case OutputInfo::DispZ:
          getPrimalScalar(i, mergedDis, numNodes, 2, time);
          break;
        case OutputInfo::RotX:
          getPrimalScalar(i, mergedDis, numNodes, 3, time);
          break;
        case OutputInfo::RotY:
          getPrimalScalar(i, mergedDis, numNodes, 4, time);
          break;
        case OutputInfo::RotZ:
          getPrimalScalar(i, mergedDis, numNodes, 5, time);
          break;
        case OutputInfo::DispMod:
          if(oinfo[i].nodeNumber == -1) {
            if(!globVal) globVal = new Scalar[numNodes];
            for(inode=0; inode<numNodes; ++inode) {
              globVal[inode] = ScalarTypes::sqrt(mergedDis[inode][0]*mergedDis[inode][0] +
                                                 mergedDis[inode][1]*mergedDis[inode][1] +
                                                 mergedDis[inode][2]*mergedDis[inode][2]);
            }
            geoSource->outputNodeScalars(i, globVal, numNodes, time);
          }
          else {
            inode = oinfo[i].nodeNumber;
            Scalar dm = ScalarTypes::sqrt(mergedDis[inode][0]*mergedDis[inode][0] +
                                          mergedDis[inode][1]*mergedDis[inode][1] +
                                          mergedDis[inode][2]*mergedDis[inode][2]);
            geoSource->outputNodeScalars(i, &dm, 1, time);
          }
          break;
	case OutputInfo::RotMod:
          if(oinfo[i].nodeNumber == -1) {
            if(!globVal) globVal = new Scalar[numNodes];
            for(inode=0; inode<numNodes; ++inode) {
              globVal[inode] = ScalarTypes::sqrt(mergedDis[inode][3]*mergedDis[inode][3] +
                                                 mergedDis[inode][4]*mergedDis[inode][4] +
                                                 mergedDis[inode][5]*mergedDis[inode][5]);
            }
            geoSource->outputNodeScalars(i, globVal, numNodes, time);
          }
          else {
            inode = oinfo[i].nodeNumber;
            Scalar rm = ScalarTypes::sqrt(mergedDis[inode][3]*mergedDis[inode][3] +
                                          mergedDis[inode][4]*mergedDis[inode][4] +
                                          mergedDis[inode][5]*mergedDis[inode][5]);
            geoSource->outputNodeScalars(i, &rm, 1, time);
          }
          break;
        case OutputInfo::TotMod:
          if(oinfo[i].nodeNumber == -1) {
            if(!globVal) globVal = new Scalar[numNodes];
            for(inode=0; inode<numNodes; ++inode) {
              globVal[inode] = ScalarTypes::sqrt(mergedDis[inode][0]*mergedDis[inode][0] +
                                                 mergedDis[inode][1]*mergedDis[inode][1] +
                                                 mergedDis[inode][2]*mergedDis[inode][2] +
                                                 mergedDis[inode][3]*mergedDis[inode][3] +
                                                 mergedDis[inode][4]*mergedDis[inode][4] +
                                                 mergedDis[inode][5]*mergedDis[inode][5]);
            }
            geoSource->outputNodeScalars(i, globVal, numNodes, time);
          }
          else {
            inode = oinfo[i].nodeNumber;
            Scalar tm = ScalarTypes::sqrt(mergedDis[inode][0]*mergedDis[inode][0] +
                                          mergedDis[inode][1]*mergedDis[inode][1] +
                                          mergedDis[inode][2]*mergedDis[inode][2] +
                                          mergedDis[inode][3]*mergedDis[inode][3] +
                                          mergedDis[inode][4]*mergedDis[inode][4] +
                                          mergedDis[inode][5]*mergedDis[inode][5]);
            geoSource->outputNodeScalars(i, &tm, 1, time);
          }
          break;
/*
        case OutputInfo::Rigid: {
#ifdef DISTRIBUTED
          int glNumRBM = (fetiSolver) ? fetiSolver->numRBM() : 0;
          if (glNumRBM) {
            int rbmSize = glNumRBM*subDomain[0]->numUncon();
            Scalar *localRBMs = new Scalar[rbmSize];
            fetiSolver->getRBMs(localRBMs);
            VectorSet globalRBMs(glNumRBM,numNodes*6, 0.0);
            subDomain[0]->expandRBM(localRBMs,globalRBMs);
            int iRBM;
            for (iRBM=0; iRBM<glNumRBM; ++iRBM)
              communicator->globalSum(numNodes*6,globalRBMs[iRBM].data());

            for (iRBM=0; iRBM<glNumRBM; ++iRBM) {
              filePrint(oinfo[i].filptr,"%e\n",0.0);
              for (iNode=0; iNode<numNodes; ++iNode)
                filePrint(oinfo[i].filptr,"%f %f %f\n",
                   globalRBMs[iRBM][6*iNode+0], globalRBMs[iRBM][6*iNode+1],
                                                globalRBMs[iRBM][6*iNode+2]);
            }
          }
#endif
          }
          break;
*/
        case OutputInfo::Energies: {
            Wext = f * u;
            if(aeroF) Waero = (*aeroF) * u;
            Wela = 0.5 * Wext;
            Scalar error = Wext+Wela+Wkin+Wdmp;
            geoSource->outputEnergies(i, time, Wext, Waero, Wela, Wkin, Wdmp, error);
          }
          break;
        case OutputInfo::Farfield: 
          domain->nffp = oinfo[i].interval;
          buildFFP(u,oinfo[i].filptr);
          break;
        case OutputInfo::AeroXForce:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 0, time);
          break;
        case OutputInfo::AeroYForce:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 1, time);
          break;
        case OutputInfo::AeroZForce:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 2, time);
          break;
        case OutputInfo::AeroXMom:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 3, time);
          break;
        case OutputInfo::AeroYMom:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 4, time);
          break;
        case OutputInfo::AeroZMom:
          if(aeroF) getAeroForceScalar(i, mergedAeroF, numNodes, 5, time);
          break;
        case OutputInfo::YModulus:
          this->getElementAttr(i,YOUNG, time);
          break;
        case OutputInfo::MDensity:
          this->getElementAttr(i,MDENS, time);
          break;
        case OutputInfo::Thicknes:
          this->getElementAttr(i,THICK, time);
          break;
        case OutputInfo::TDEnforcement: {
          double *plot_data = new double[numNodes]; 
          for(int iNode=0; iNode<numNodes; ++iNode) plot_data[iNode] = 0.0;
          for(int iMortar=0; iMortar<domain->GetnMortarConds(); iMortar++) {
            domain->GetMortarCond(iMortar)->get_plot_variable(oinfo[i].tdenforc_var,plot_data);
          }
          if(oinfo[i].nodeNumber == -1) 
            geoSource->outputNodeScalars(i, plot_data, numNodes, time);
          else 
            geoSource->outputNodeScalars(i, &plot_data[oinfo[i].nodeNumber], 1, time);
          delete [] plot_data;
        } break;
        default:
          filePrint(stderr," *** WARNING: Output case %d not implemented \n", i);
          break;
      }
    }
    if(globVal) { delete [] globVal; globVal = 0; }
  }
  firstOutput = false; 
  // ... CALCULATE STRUCTURE MASS IF REQUESTED
  if(domain->solInfo().massFlag)  {
    double mass = domain->computeStructureMass();
    filePrint(stderr," ... Structure mass = %10.4f    ...\n",mass);
    filePrint(stderr," --------------------------------------\n");
  }

  if(mergedDis) delete [] mergedDis;
  if(aeroF) delete [] mergedAeroF;
  if(distState) { delete [] mergedVel; delete [] mergedAcc; }
  if(globVal) delete [] globVal; 

}

template<class Scalar>
void
GenDecDomain<Scalar>::getPrimalVector(int fileNumber, Scalar (*xyz)[11], int numNodes,
                                      int ndof, double time)//DofSet::max_known_nonL_dof
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];

  int inode;
  if(oinfo.nodeNumber == -1) {
    if(ndof == 6) geoSource->outputNodeVectors6(fileNumber, xyz, numNodes, time);
    else geoSource->outputNodeVectors(fileNumber, xyz, numNodes, time);
  }
  else {
    inode = oinfo.nodeNumber;
    if(ndof == 6) geoSource->outputNodeVectors6(fileNumber, xyz+inode, 1, time);
    else geoSource->outputNodeVectors(fileNumber, xyz+inode, 1, time);
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::getPrimalScalar(int fileNumber, Scalar (*xyz)[11], int numNodes, 
                                      int dof, double time)//DofSet::max_known_nonL_dof
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber]; 

  int inode;
  if(oinfo.nodeNumber == -1) {
    Scalar *xyz_dof = new Scalar[numNodes];
    for(inode = 0; inode < numNodes; ++inode) xyz_dof[inode] = xyz[inode][dof];
    geoSource->outputNodeScalars(fileNumber, xyz_dof, numNodes, time);
  }
  else {
    inode = oinfo.nodeNumber;
    geoSource->outputNodeScalars(fileNumber, xyz[inode]+dof, 1, time);
  }
} 

template<class Scalar>
void
GenDecDomain<Scalar>::getAeroForceScalar(int fileNumber, Scalar (*mergedAeroF)[6],
                                         int numNodes, int dof, double time)
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];

  int inode;
  if(oinfo.nodeNumber == -1) {
    geoSource->outputNodeScalars(fileNumber, (double*)0, 0, time); // output time
    for(inode = 0; inode < numNodes; ++inode)
      geoSource->outputNodeScalars(fileNumber, mergedAeroF[inode]+dof, 1);
  }
  else {
    inode = oinfo.nodeNumber;
    geoSource->outputNodeScalars(fileNumber, mergedAeroF[inode]+dof, 1, time);
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::computeSubdElemForce(int iSub, Scalar *globForce,
                                          GenDistrVector<Scalar> *u, int Findex)
{
  Scalar *locForce = new Scalar[subDomain[iSub]->countElemNodes()];
  subDomain[iSub]->computeElementForce(u->subData(iSub), Findex, locForce);
  subDomain[iSub]->mergeElemStress(locForce, globForce, elemToNode);
  delete [] locForce;
}

template<class Scalar>
void GenDecDomain<Scalar>::getElementForce(GenDistrVector<Scalar> &u, int fileNumber,
                                           int Findex, double time)
{
  int numElemNodes = elemToNode->numConnect();
  Scalar *globForce = new Scalar[numElemNodes];
  execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdElemForce,
            globForce, &u, Findex);
  geoSource->outputElemStress(fileNumber, globForce, elemToNode->csize(),
                              elemToNode->ptr(), time);
  delete [] globForce;
}

template<class Scalar>
void
GenDecDomain<Scalar>::computeSubdStress(int iSub, GenDistrVector<Scalar> *globStress, 
                                        GenDistrVector<Scalar> *globWeight, GenDistrVector<Scalar> *u, 
                                        int fileNumber, int Findex)
{
  GenStackVector<Scalar> gstress(globStress->subData(iSub),globStress->subLen(iSub));
  GenStackVector<Scalar> gweight(globWeight->subData(iSub),globWeight->subLen(iSub));
  subDomain[iSub]->computeStressStrain(fileNumber, u->subData(iSub), 
	                               Findex, gstress.data(), gweight.data());
}

template<class Scalar>
void GenDecDomain<Scalar>::computeSubdElemStress(int iSub, Scalar *glElemStress,
                                                 GenDistrVector<Scalar> *u, int fileNumber, int Findex) 
{
  Scalar *locStress = new Scalar[subDomain[iSub]->countElemNodes()];
  subDomain[iSub]->computeStressStrain(fileNumber, u->subData(iSub),
                                       Findex, locStress);
  subDomain[iSub]->mergeElemStress(locStress, glElemStress, elemToNode);
  delete [] locStress;
}

template<class Scalar>
void
GenDecDomain<Scalar>::computeSubdStress(int iSub, GenDistrVector<Scalar> *globStress,
                                        GenDistrVector<Scalar> *globWeight, DistrGeomState *u, 
                                        Corotator ***allCorot, int *fileNumber, int *Findex)
{
  // Non-linear version of computeSubdStress
  subDomain[iSub]->computeStressStrain((*u)[iSub], allCorot[iSub],
                                       *fileNumber, *Findex, globStress->subData(iSub),
                                       globWeight->subData(iSub));
}

template<class Scalar>
void GenDecDomain<Scalar>::computeSubdElemStress(int iSub, Scalar *glElemStress,
                                                 DistrGeomState *u, Corotator ***allCorot, 
                                                 int fileNumber, int Findex)
{
  // Non-linear version of computeSubdElemStress
  Scalar *locStress = new Scalar[subDomain[iSub]->countElemNodes()];
  subDomain[iSub]->computeStressStrain((*u)[iSub], allCorot[iSub], fileNumber,
                                       Findex, locStress);
  subDomain[iSub]->mergeElemStress(locStress, glElemStress, elemToNode);
  delete [] locStress;
}

template<class Scalar>
void GenDecDomain<Scalar>::getElementStressStrain(DistrGeomState *gs, Corotator ***allCorot, 
                                                  int fileNumber, int Findex, double time)
{
  // Non-linear version of getElementStressStrain
  int numElemNodes = elemToNode->numConnect();
  Scalar *globStress = new Scalar[numElemNodes];
  execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdElemStress,
            globStress, gs, allCorot, fileNumber, Findex);
  geoSource->outputElemStress(fileNumber, globStress, elemToNode->csize(),
                              elemToNode->ptr(), time);
  delete [] globStress;
}

template<class Scalar>
void
GenDecDomain<Scalar>::getStressStrain(DistrGeomState *gs, Corotator ***allCorot,
                                      int fileNumber, int Findex, double time) 
{
 // Non-linear version of getStressStrain
 OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
 if(oinfo.averageFlg == 0) {
   getElementStressStrain(gs, allCorot, fileNumber, Findex, time);
   return;
 }

 // Allocate a distributed vector and initialize it to zero
 // if it hasn't already been allocated.
 if(stress == 0) stress = new GenDistrVector<Scalar>(nodeInfo);
 if(weight == 0) weight = new GenDistrVector<Scalar>(nodeInfo);

 stress->zero();
 weight->zero();

 // each subdomain computes its stress vector
 execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdStress,
           stress, weight, gs, allCorot, &fileNumber, &Findex);

 int numNodes = geoSource->getNumGlobNodes();
 if(globalStress == 0) globalStress = new Scalar[numNodes]; 
 if(globalWeight == 0) globalWeight = new Scalar[numNodes];

 int i;
 for(i=0; i<numNodes; ++i)
   globalStress[i] = globalWeight[i] = 0.0;

 int iSub;
 for(iSub=0; iSub<numSub; ++iSub)
   subDomain[iSub]->mergeStress(stress->subData(iSub), weight->subData(iSub),
                                globalStress, globalWeight);

 for(i=0; i < numNodes; ++i)  {
   if(globalWeight[i] == 0.0)
     globalStress[i] = 0.0;
   else
     globalStress[i] = globalStress[i]/globalWeight[i];
 }

 geoSource->outputNodeScalars(fileNumber, globalStress, numNodes, time);

 delete [] globalWeight; globalWeight=0;
 delete [] globalStress; globalStress=0;

 delete stress; stress = 0;
 delete weight; weight = 0;
}

template<class Scalar>
void GenDecDomain<Scalar>::setsizeSfemStress(int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  if(avgnum == 1)  sizeSfemStress = geoSource->getNumGlobNodes();  // node-based output
  else if(avgnum == 0) {  // element-based output
   sizeSfemStress = 0;
/*   Connectivity *elemToNode = new Connectivity(domain->getEset());
   int numele = geoSource->getNumAttributes();  // number of elements; another option domain->numElements();
   for(int iele=0; iele<numele; ++iele)   {
//     cerr << "number of nodes in this element  = " << elemToNode->num(iele) << endl;
     sizeSfemStress = sizeSfemStress + elemToNode->num(iele); // add number of nodes for each element
   }*/
  }
  else {
   cerr << "avgnum = " << avgnum << " not implemented in Domain::setsizeSfemStress()" << endl;
   sizeSfemStress = 0;
  }
}


template<class Scalar>
void GenDecDomain<Scalar>::updateSfemStress(Scalar* str, int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;
  int numNodes = geoSource->getNumGlobNodes();
  if(avgnum == 1)  for (int i=0;i<numNodes;++i) globalStress[i] = str[i];
  else if(avgnum == 0) cerr << "updateSfemStress for element not yet implemented" << endl; // for (int i=0;i<stressAllElems->size();++i) (*stressAllElems) = str[i]; // YYY DG
  else {cerr << "avgnum = " << avgnum << " not implemented in Domain::updateSfemStress()" << endl;}
}


template<class Scalar>
void GenDecDomain<Scalar>::getElementStressStrain(GenDistrVector<Scalar> &u, int fileNumber,
                                                  int Findex, double time, int printFlag)  
{
  // allocate arrays
  int numElemNodes = elemToNode->numConnect();
  Scalar *globStress = new Scalar[numElemNodes];
  // each subdomain computes its stress vector
  execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdElemStress,
            globStress, &u, fileNumber, Findex);
  geoSource->outputElemStress(fileNumber, globStress, elemToNode->csize(),
                              elemToNode->ptr(), time);
}

template<class Scalar>
void GenDecDomain<Scalar>::getStressStrain(GenDistrVector<Scalar> &u, int fileNumber, 
           				   int Findex, double time, int printFlag)  
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementStressStrain(u, fileNumber, Findex, time, printFlag);
    return;
  }

  // Allocate a distributed vector and initialize it to zero
  // if it hasn't already been allocated.
  if(stress == 0) stress = new GenDistrVector<Scalar>(nodeInfo);
  if(weight == 0) weight = new GenDistrVector<Scalar>(nodeInfo);

  stress->zero();
  weight->zero();
  

  if(printFlag != 2) {
    // each subdomain computes its stress vector
    if(Findex != 16)
      execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdStress,
                stress, weight, &u, fileNumber, Findex); // YYY DG Does it print ? Probably not
  }


  // allocate global stress and weight arrays 
  int numNodes = geoSource->getNumGlobNodes();
  if(globalStress == 0) globalStress = new Scalar[numNodes]; 
  if(globalWeight == 0) globalWeight = new Scalar[numNodes];

  if (printFlag != 2) { 
    int i;
    for (i = 0; i < numNodes; ++i)
      globalStress[i] = globalWeight[i] = 0.0;
 

    int iSub;
    for(iSub=0; iSub < numSub; ++iSub) {
      if(Findex != 16) {
        subDomain[iSub]->mergeStress(stress->subData(iSub), weight->subData(iSub),
                                     globalStress, globalWeight);
      }
      else {
        subDomain[iSub]->computeContactPressure(globalStress, globalWeight); // PJSA
      }
    }
    for(i = 0; i < numNodes; ++i)  {
      if(globalWeight[i] == 0.0)
        globalStress[i] = 0.0;
      else
        globalStress[i] = globalStress[i]/globalWeight[i];
    }

  }

  if(printFlag != 1) geoSource->outputNodeScalars(fileNumber, globalStress, numNodes, time); // YYY DG Probably the only printing function

  delete stress; stress = 0;
  delete weight; weight = 0;
}

template<class Scalar>
void
GenDecDomain<Scalar>::getPrincipalStress(DistrGeomState *gs, Corotator ***allCorot,
                                         int fileNumber, int strIndex, double time)
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementPrincipalStress(gs, allCorot, fileNumber, strIndex, time);
    return;
  }

  // set stress VS. strain for element subroutines (put this in a function!!)
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];
  bool direction;

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i) strDir[i] = i;
    direction = false;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i) strDir[i] = i+7;
    direction = false;
  }
  else if((strIndex==6) || (strIndex==7) || (strIndex==8)) {
    strInd = strIndex-6;
    stressORstrain = 0;
    for(i=0; i<6; ++i) strDir[i] = i;
    direction = true;
  }
  else if((strIndex==9) || (strIndex==10) || (strIndex==11)) {
    strInd = strIndex-9;
    stressORstrain = 1;
    for(i=0; i<6; ++i) strDir[i] = i+7;
    direction = true;
  }
  else {
    filePrint(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // Allocate a distributed vector for stress if it hasn't already been allocated
  if(stress == 0) stress = new GenDistrVector<Scalar>(nodeInfo);
  if(weight == 0) weight = new GenDistrVector<Scalar>(nodeInfo);

  // stress storage
  int n = domain->numnodes;
  Scalar (*globalAllStress)[6] = new Scalar [n][6];

  // Compute Each Required Stress (all 6) using same routines as for 
  // individual stresses
  int str_loop;
  int Findex;
  for (str_loop = 0; str_loop < 6; ++str_loop) {
    // get current stress/strain index
    Findex = strDir[str_loop];
    // Initialize distributed vector to zero
    stress->zero();
    weight->zero();
    // each subdomain computes its stress vector
    execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdStress,
             stress, weight, gs, allCorot, &fileNumber, &Findex);	
    Scalar *globalStress = new Scalar[domain->numnodes]; 
    Scalar *globalWeight = new Scalar[domain->numnodes];
    for(i = 0; i < domain->numnodes; ++i)
      globalStress[i] = globalWeight[i] = 0.0;
    int iSub;
    for(iSub=0; iSub<numSub; ++iSub)
      subDomain[iSub]->mergeStress(stress->subData(iSub),
                                   weight->subData(iSub),
                                   globalStress,globalWeight);
    for(i = 0; i < domain->numnodes; ++i)
      if(globalWeight[i] != 0.0)
        globalAllStress[i][str_loop] = globalStress[i]/globalWeight[i];
      else
        globalAllStress[i][str_loop] = 0.0;
    delete [] globalWeight; globalWeight=0;
    delete [] globalStress; globalStress=0;
  }

  // ... CALCULATE PRINCIPALS AT EACH NODE
  // PJSA 3-24-05: modified to compute principal direction if required
  Scalar svec[6], pvec[3];
  Scalar *globalPVec = 0; 
  Scalar (*globalPDir)[11] = 0;//DofSet::max_known_nonL_dof 
  Scalar *pdir = 0;
  int numNodes = domain->numNode();
  if(direction) {
    globalPDir = new Scalar[numNodes][11];//DofSet::max_known_nonL_dof 
    pdir = (Scalar *) dbg_alloca(sizeof(Scalar)*11);
  }
  else globalPVec = new Scalar[numNodes];

  for(i = 0; i < numNodes; ++i) {
    for(j = 0; j < 6; ++j)
      svec[j] = globalAllStress[i][j];
    // Convert Engineering to Tensor Strains
    if(stressORstrain != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec, pvec, pdir);
    if(direction) for(j=0;j<3;++j) globalPDir[i][j] = pdir[strInd*3+j];
    else globalPVec[i] = pvec[strInd];
  }
  if(direction) {
    if(oinfo.nodeNumber == -1) 
       geoSource->outputNodeVectors(fileNumber, globalPDir, numNodes, time);
    else 
      geoSource->outputNodeVectors(fileNumber, globalPDir+oinfo.nodeNumber, 1, time);
  }
  else {
    if(oinfo.nodeNumber == -1)
      geoSource->outputNodeScalars(fileNumber, globalPVec, numNodes, time);
    else
      geoSource->outputNodeScalars(fileNumber, globalPVec+oinfo.nodeNumber, 1, time);
  }

  delete [] globalAllStress;
  if(globalPVec) delete [] globalPVec;
  if(globalPDir) delete [] globalPDir;
}

template<class Scalar>
void GenDecDomain<Scalar>::getElementPrincipalStress(GenDistrVector<Scalar> &u, int fileNumber,
                                                      int strIndex, double time)
{
  // set stress VS. strain for element subroutines
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];
                                                                                                                               
  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  else {
    filePrint(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // allocate arrays
  int numElemNodes = elemToNode->numConnect();
  Scalar **globAllStress = new Scalar * [6];
  for(i=0; i<6; ++i) globAllStress[i] = new Scalar[numElemNodes];

  // Compute Each Required Stress (all 6) using same routines as for
  // individual stresses
  int str_loop;
  int Findex;
  for(str_loop = 0; str_loop < 6; ++str_loop) {
    // get current stress/strain index
    Findex = strDir[str_loop];

    // each subdomain computes its stress vector
    execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdElemStress,
              globAllStress[str_loop], &u, fileNumber, Findex);
  }

  // ... CALCULATE PRINCIPALS AT EACH NODE
  Scalar svec[6], pvec[3];
  Scalar *globalPVec = new Scalar[numElemNodes];
  for(i = 0; i < numElemNodes; ++i) {
    for(j = 0; j < 6; ++j)
      svec[j] = globAllStress[j][i];
    // Convert Engineering to Tensor Strains
    if(stressORstrain != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec, pvec);
    globalPVec[i] = pvec[strInd];
  }
  geoSource->outputElemStress(fileNumber, globalPVec, elemToNode->csize(),
                              elemToNode->ptr(), time);

  for(i=0; i<6; ++i) delete [] globAllStress[i];
  delete [] globAllStress;
  delete [] globalPVec;
}

template<class Scalar>
void
GenDecDomain<Scalar>::getElementPrincipalStress(DistrGeomState *gs, Corotator ***allCorot,
                                                int fileNumber, int strIndex, double time)
{
  // PJSA 3-23-05 Non-linear version of getElementPrincipalStress
  // set stress VS. strain for element subroutines
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  else {
    filePrint(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // allocate arrays
  int numElemNodes = elemToNode->numConnect();
  Scalar **globAllStress = new Scalar * [6];
  for(i=0; i<6; ++i) globAllStress[i] = new Scalar[numElemNodes];

  // Compute Each Required Stress (all 6) using same routines as for
  // individual stresses
  int str_loop;
  int Findex;
  for(str_loop = 0; str_loop < 6; ++str_loop) {
    // get current stress/strain index
    Findex = strDir[str_loop];

    // each subdomain computes its stress vector
    execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdElemStress,
              globAllStress[str_loop], gs, allCorot, fileNumber, Findex);
  }

  // ... CALCULATE PRINCIPALS AT EACH NODE
  Scalar svec[6], pvec[3];
  Scalar *globalPVec = new Scalar[numElemNodes];
  for(i = 0; i < numElemNodes; ++i) {
    for(j = 0; j < 6; ++j)
      svec[j] = globAllStress[j][i];
    // Convert Engineering to Tensor Strains
    if(stressORstrain != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec, pvec);
    globalPVec[i] = pvec[strInd];
  }
  geoSource->outputElemStress(fileNumber, globalPVec, elemToNode->csize(),
                              elemToNode->ptr(), time);

  for(i=0; i<6; ++i) delete [] globAllStress[i];
  delete [] globAllStress;
  delete [] globalPVec;
}


template<class Scalar>
void
GenDecDomain<Scalar>::getPrincipalStress(GenDistrVector<Scalar> &u, int fileNumber, int strIndex, 
			     	         double time)  
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementPrincipalStress(u, fileNumber, strIndex, time);
    return;
  }

  // set stress VS. strain for element subroutines (put this in a function!!)
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];
  bool direction;

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i) strDir[i] = i;
    direction = false;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i) strDir[i] = i+7;
    direction = false;
  }
  else if((strIndex==6) || (strIndex==7) || (strIndex==8)) {
    strInd = strIndex-6;
    stressORstrain = 0;
    for(i=0; i<6; ++i) strDir[i] = i;
    direction = true;
  }
  else if((strIndex==9) || (strIndex==10) || (strIndex==11)) {
    strInd = strIndex-9;
    stressORstrain = 1;
    for(i=0; i<6; ++i) strDir[i] = i+7;
    direction = true;
  }
  else {
    filePrint(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // Allocate a distributed vector for stress if it hasn't already been allocated
  if(stress == 0) stress = new GenDistrVector<Scalar>(nodeInfo);
  if(weight == 0) weight = new GenDistrVector<Scalar>(nodeInfo);

  // stress storage
  int n = domain->numnodes;
  Scalar (*globalAllStress)[6] = new Scalar [n][6];

  // Compute Each Required Stress (all 6) using same routines as for 
  // individual stresses
  int str_loop;
  int Findex;
  for (str_loop = 0; str_loop < 6; ++str_loop) {
    // get current stress/strain index
    Findex = strDir[str_loop];
    // Initialize distributed vector to zero
    stress->zero();
    weight->zero();
    // each subdomain computes its stress vector
    execParal(numSub, this, &GenDecDomain<Scalar>::computeSubdStress,
              stress, weight, &u, fileNumber, Findex);
    Scalar *globalStress = new Scalar[domain->numnodes]; 
    Scalar *globalWeight = new Scalar[domain->numnodes];
    for(i = 0; i < domain->numnodes; ++i)
      globalStress[i] = globalWeight[i] = 0.0;
    int iSub;
    for(iSub=0; iSub<numSub; ++iSub)
      subDomain[iSub]->mergeStress(stress->subData(iSub),
                                   weight->subData(iSub),
                                   globalStress,globalWeight);
    for(i = 0; i < domain->numnodes; ++i)
      if(globalWeight[i] != 0.0)
        globalAllStress[i][str_loop] = globalStress[i]/globalWeight[i];
      else
        globalAllStress[i][str_loop] = 0.0;
    delete [] globalWeight; globalWeight=0;
    delete [] globalStress; globalStress=0;
  }

  // ... CALCULATE PRINCIPALS AT EACH NODE
  // PJSA 3-24-05: modified to compute principal direction if required
  Scalar svec[6], pvec[3];
  Scalar *globalPVec = 0; 
  Scalar (*globalPDir)[11] = 0;//DofSet::max_known_nonL_dof 
  Scalar *pdir = 0;
  int numNodes = domain->numNode();
  if(direction) {
    globalPDir = new Scalar[numNodes][11];//DofSet::max_known_nonL_dof 
    pdir = (Scalar *) dbg_alloca(sizeof(Scalar)*11);
  }
  else globalPVec = new Scalar[numNodes];

  for(i = 0; i < numNodes; ++i) {
    for(j = 0; j < 6; ++j)
      svec[j] = globalAllStress[i][j];
    // Convert Engineering to Tensor Strains
    if(stressORstrain != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec, pvec, pdir);
    if(direction) for(j=0;j<3;++j) globalPDir[i][j] = pdir[strInd*3+j];
    else globalPVec[i] = pvec[strInd];
  }
  if(direction) {
    if(oinfo.nodeNumber == -1) 
       geoSource->outputNodeVectors(fileNumber, globalPDir, numNodes, time);
    else 
      geoSource->outputNodeVectors(fileNumber, globalPDir+oinfo.nodeNumber, 1, time);
  }
  else {
    if(oinfo.nodeNumber == -1)
      geoSource->outputNodeScalars(fileNumber, globalPVec, numNodes, time);
    else
      geoSource->outputNodeScalars(fileNumber, globalPVec+oinfo.nodeNumber, 1, time);
  }

  delete [] globalAllStress;
  if(globalPVec) delete [] globalPVec;
  if(globalPDir) delete [] globalPDir;
}

// to output primal error
template<class Scalar>
void
GenDecDomain<Scalar>::outputPrimal(GenDistrVector<Scalar> &primal, int iter)
{
 int numNode = domain->numnodes;

 // Open primal error file
 if(primalFile == 0) {
   if((primalFile=fopen("primalError","w"))==(FILE *) 0 )
     filePrint(stderr," *** ERROR: Cannot open %s ...\n", "primalError");
   fflush(primalFile);

   filePrint(primalFile,"Scalar PRIMAL under FETI for %s\n%d\n",
             geoSource->getCheckFileInfo()->nodeSetName, numNode);
 }

 Scalar *error = (Scalar*) dbg_alloca(numNode*sizeof(Scalar));
 int i;
 for(i=0; i<numNode; ++i)
   error[i] = 0.0;

 int isub;
 for(isub=0; isub<numSub; ++isub) {
   subDomain[isub]->mergePrimalError(error, primal.subData(isub));
 }

 // print iteration number
 filePrint(primalFile,"%d\n",iter);

 // print primal error at each node
 for(i=0; i<numNode; ++i)
   filePrint(primalFile,"%e\n",error[i]);
}

// -----------------------------
// Nonlinear DecDomain functions
// -----------------------------

template<class Scalar>
void
GenDecDomain<Scalar>::postProcessing(DistrGeomState *geomState, Corotator ***allCorot,
                                     double x, SysState<GenDistrVector<Scalar> > *distState)  
{
  // NOTE: for dynamic runs, x represents the time
  //       for static runs, x represents the load parameter, lambda
  int numOutInfo = geoSource->getNumOutInfo();
  if(verboseFlag && numOutInfo && x == 0)
    filePrint(stderr," ... Postprocessing                 ...\n");

  int numNodes = geoSource->numNode();
  Scalar (*xyz)[11] = new Scalar[numNodes][11];//DofSet::max_known_nonL_dof
  Scalar *globVal = 0;  // for output

  int i,j,iSub;
  for(i = 0; i < numNodes; ++i)
    for (j = 0 ; j < 11 ; j++)
    xyz[i][j] = 0.0;

  int isub;
  for(isub = 0; isub < numSub; ++isub)
    subDomain[isub]->mergeDisp(xyz,(*geomState)[isub]);

  // for nonlinear dynamics: initialize and merge velocities and accelerations from subdomains into global array
  GenDistrVector<Scalar> *v_n = 0, *a_n = 0;
  Scalar (*mergedVel)[11] = 0, (*mergedAcc)[11] = 0;//DofSet::max_known_nonL_dof
  if(distState) {
    v_n = &distState->getVeloc();
    a_n = &distState->getAccel();
    mergedVel = new Scalar[numNodes][11];//DofSet::max_known_nonL_dof
    mergedAcc = new Scalar[numNodes][11];//DofSet::max_known_nonL_dof
    for(i = 0; i < numNodes; ++i)
      for(j=0; j<11; ++j) mergedVel[i][j] = mergedAcc[i][j] = 0.0; //DofSet::max_known_nonL_dof
    for(iSub = 0; iSub < numSub; ++iSub) {
      subDomain[iSub]->mergeAllVeloc(mergedVel, v_n->subData(iSub));
      subDomain[iSub]->mergeAllAccel(mergedAcc, a_n->subData(iSub));
    }
  }

  if(x == 0) {
    geoSource->openOutputFiles();
  }

  int inode;
  OutputInfo *oinfo = geoSource->getOutputInfo();
  for(i = 0; i < numOutInfo; i++) {
   int step = (domain->solInfo().isDynam()) ? int(x/domain->solInfo().getTimeStep()+0.5) : int(x/domain->solInfo().getNLInfo().dlambda+0.5);
   //cerr << "i = " << i << ", x = " << x << ", step = " << step << ", interval = " << oinfo[i].interval << endl;
   if(oinfo[i].interval != 0 && step % oinfo[i].interval == 0) {
    // int dof = -1;
    switch(oinfo[i].type) {
     case OutputInfo::FreqRespModes:
     case OutputInfo::Displacement:
       getPrimalVector(i, xyz, numNodes, 3, x);
       break;
     case OutputInfo::Velocity:
       if(distState) getPrimalVector(i, mergedVel, numNodes, 3, x);
       break;
     case OutputInfo::Acceleration:
       if(distState) getPrimalVector(i, mergedAcc, numNodes, 3, x);
       break;
     case OutputInfo::Disp6DOF:
       getPrimalVector(i, xyz, numNodes, 6, x);
       break;
     case OutputInfo::Velocity6:
       if(distState) getPrimalVector(i, mergedVel, numNodes, 6, x);
       break;
     case OutputInfo::Accel6:
       if(distState) getPrimalVector(i, mergedAcc, numNodes, 6, x);
       break;
     case OutputInfo::Temperature:
       getPrimalScalar(i, xyz, numNodes, 0, x);
       break;
     case OutputInfo::TemperatureFirstTimeDerivative:
       if(distState) getPrimalScalar(i, mergedVel, numNodes, 6, x);
       break;
     case OutputInfo::PressureFirstTimeDerivative:
       if(distState) getPrimalScalar(i, mergedVel, numNodes, 7, x);
       break;
     case OutputInfo::PressureSecondTimeDerivative:
       if(distState) getPrimalScalar(i, mergedAcc, numNodes, 7, x);
       break;
     case OutputInfo::StressXX:
       getStressStrain(geomState, allCorot, i, SXX, x);
       break;
     case OutputInfo::StressYY:
       getStressStrain(geomState, allCorot, i, SYY, x);
       break;
     case OutputInfo::StressZZ:
       getStressStrain(geomState, allCorot, i, SZZ, x);
       break;
     case OutputInfo::StressXY:
       getStressStrain(geomState, allCorot, i, SXY, x);
       break;
     case OutputInfo::StressYZ:
       getStressStrain(geomState, allCorot, i, SYZ, x);
       break;
     case OutputInfo::StressXZ:
       getStressStrain(geomState, allCorot, i, SXZ, x);
       break;
     case OutputInfo::StrainXX:
       getStressStrain(geomState, allCorot, i, EXX, x);
       break;
     case OutputInfo::StrainYY:
       getStressStrain(geomState, allCorot, i, EYY, x);
       break;
     case OutputInfo::StrainZZ:
       getStressStrain(geomState, allCorot, i, EZZ, x);
       break;
     case OutputInfo::StrainXY:
       getStressStrain(geomState, allCorot, i, EXY, x);
       break;
     case OutputInfo::StrainYZ:
       getStressStrain(geomState, allCorot, i, EYZ, x);
       break;
     case OutputInfo::StrainXZ:
       getStressStrain(geomState, allCorot, i, EXZ, x);
       break;
     case OutputInfo::StressVM:
       getStressStrain(geomState, allCorot, i, VON, x);
       break;
     case OutputInfo::StrainVM:
       getStressStrain(geomState, allCorot,i,STRAINVON, x);
       break;
     case OutputInfo::StressPR1:
       getPrincipalStress(geomState, allCorot, i, PSTRESS1, x);
       break;
     case OutputInfo::StressPR2:
       getPrincipalStress(geomState, allCorot, i, PSTRESS2, x);
       break;
     case OutputInfo::StressPR3:
       getPrincipalStress(geomState, allCorot, i, PSTRESS3, x);
       break;
     case OutputInfo::StrainPR1:
       getPrincipalStress(geomState, allCorot, i, PSTRAIN1, x);
       break;
     case OutputInfo::StrainPR2:
       getPrincipalStress(geomState, allCorot, i, PSTRAIN2, x);
       break;
     case OutputInfo::StrainPR3:
       getPrincipalStress(geomState, allCorot, i, PSTRAIN3, x);
       break;
     case OutputInfo::DispX:
       getPrimalScalar(i, xyz, numNodes, 0, x);
       break;
     case OutputInfo::DispY:
       getPrimalScalar(i, xyz, numNodes, 1, x);
       break;
     case OutputInfo::DispZ:
       getPrimalScalar(i, xyz, numNodes, 2, x);
       break;
     case OutputInfo::RotX:
       getPrimalScalar(i, xyz, numNodes, 3, x);
       break;
     case OutputInfo::RotY:
       getPrimalScalar(i, xyz, numNodes, 4, x);
       break;
     case OutputInfo::RotZ:
       getPrimalScalar(i, xyz, numNodes, 5, x);
       break;
     case OutputInfo::DispMod:
       if(!globVal) globVal = new Scalar[numNodes];
       for(inode=0; inode<numNodes; ++inode) {
         globVal[inode] = ScalarTypes::sqrt(xyz[inode][0]*xyz[inode][0] +
                                            xyz[inode][1]*xyz[inode][1] +
                                            xyz[inode][2]*xyz[inode][2]);
       }
       geoSource->outputNodeScalars(i, globVal, numNodes, x);
       break;
     case OutputInfo::RotMod:
       if(!globVal) globVal = new Scalar[numNodes];
       for(inode=0; inode<numNodes; ++inode) {
         globVal[inode] = ScalarTypes::sqrt(xyz[inode][3]*xyz[inode][3] +
                                            xyz[inode][4]*xyz[inode][4] +
                                            xyz[inode][5]*xyz[inode][5]);
       }
       geoSource->outputNodeScalars(i, globVal, numNodes, x);
       break;
     case OutputInfo::TotMod:
       if(!globVal) globVal = new Scalar[numNodes];
       for(inode=0; inode<numNodes; ++inode) {
         globVal[inode] = ScalarTypes::sqrt(xyz[inode][0]*xyz[inode][0] +
                                            xyz[inode][1]*xyz[inode][1] +
                                            xyz[inode][2]*xyz[inode][2] +
                                            xyz[inode][3]*xyz[inode][3] +
                                            xyz[inode][4]*xyz[inode][4] +
                                            xyz[inode][5]*xyz[inode][5]);
       }
       geoSource->outputNodeScalars(i, globVal, numNodes, x);
       break;
     default:
       filePrint(stderr," *** WARNING: Output case %d not implemented for non-linear FETI\n", i);
       break;
   }
  }
 }
 if(globVal) delete [] globVal;

 // --- Print Problem statistics -------------------------------------
/*
 if(x == 0.0 || x == 1.0) {
   int numnod = domain->numnodes;
   int numele = domain->numElements();
   int numsub = getNumSub();
   int numNeuman = domain->numNeuman;
   //int numOutput = numOutInfo;
   filePrint(stderr,"\n ---------- PROBLEM PARAMETERS --------");
   filePrint(stderr,"\n ... # subdomains         = %7d ...",numsub);
   filePrint(stderr,"\n ... # nodes              = %7d ...",numnod);
   filePrint(stderr,"\n ... # elements           = %7d ...",numele);
   filePrint(stderr,"\n ... # applied loads      = %7d ...",numNeuman);
   filePrint(stderr,"\n ... # output files       = %7d ...",numOutInfo);
   filePrint(stderr,"\n --------------------------------------\n");
 }
*/
}

// element vector distributed vector info
// each element has a certain length based on its element stiffness matrix
template<class Scalar>
DistrInfo*
GenDecDomain<Scalar>::elementVectorInfo()
{
  DistrInfo *eleVecInfo = new DistrInfo(numSub);
  int totLen = 0;
  for(int iSub = 0; iSub < numSub; ++iSub) {
    eleVecInfo->domLen[iSub] = subDomain[iSub]->maxNumDOF();
    totLen += eleVecInfo->domLen[iSub];
  }
  eleVecInfo->len = totLen;
  return eleVecInfo;
}

// prescribed boundary condition distributed vector info
template<class Scalar>
DistrInfo*
GenDecDomain<Scalar>::pbcVectorInfo()
{
 DistrInfo *bcVecInfo = new DistrInfo(numSub);
 int totLen = 0;
 int iSub;
 for(iSub = 0; iSub < numSub; ++iSub) {
    bcVecInfo->domLen[iSub] = subDomain[iSub]->nDirichlet();
    totLen += bcVecInfo->domLen[iSub];
 }
 bcVecInfo->len = totLen;
 return bcVecInfo;
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeInternalInfo()
{
 startTimerMemory(mt.makeInternalInfo, mt.memoryInternal);

 // Create internal Distributed information
 internalInfo.domLen = new int[numSub];
 internalInfo.numDom = numSub;
 int totLen = 0;
 for(int iSub = 0; iSub < numSub; ++iSub) {
   internalInfo.domLen[iSub] = subDomain[iSub]->numUncon();
   totLen += internalInfo.domLen[iSub];
 }
 internalInfo.len = totLen;
 if(domain->solInfo().inpc || domain->solInfo().timeIntegration == SolverInfo::Qstatic) {
   bool *internalMasterFlag = new bool[totLen];
   internalInfo.computeOffsets();
   for(int iSub = 0; iSub < numSub; ++iSub) {
     subDomain[iSub]->computeInternalMasterFlag();
     bool *subMasterFlag = subDomain[iSub]->getInternalMasterFlag();
     int subOffset = internalInfo.subOffset[iSub];
     for(int j=0; j<internalInfo.domLen[iSub]; ++j)
       internalMasterFlag[subOffset+j] = subMasterFlag[j];
   }
   internalInfo.setMasterFlag(internalMasterFlag);
 } else
 internalInfo.setMasterFlag();
 stopTimerMemory(mt.makeInternalInfo, mt.memoryInternal);
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeInternalInfo2()
{
 startTimerMemory(mt.makeInternalInfo, mt.memoryInternal);

 // Create internal Distributed information for all dofs, both constrained and unconstrained
 internalInfo2.domLen = new int[numSub];
 internalInfo2.numDom = numSub;
 int totLen = 0;
 for(int iSub = 0; iSub < numSub; ++iSub) {
   internalInfo2.domLen[iSub] = subDomain[iSub]->numdof();
   totLen += internalInfo2.domLen[iSub];
 }
 internalInfo2.len = totLen;
 internalInfo2.setMasterFlag();
 stopTimerMemory(mt.makeInternalInfo, mt.memoryInternal);
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeNodeInfo()
{
 startTimerMemory(mt.makeInternalInfo, mt.memoryInternal);

 // Create nodal Distributed information
 nodeInfo.domLen = new int[numSub];
 nodeInfo.numDom = numSub;
 int totLenNode = 0;
 for(int iSub = 0; iSub < numSub; ++iSub) {
   nodeInfo.domLen[iSub] = subDomain[iSub]->numNodes(); // this is used for nodal stress output
   totLenNode += nodeInfo.domLen[iSub];
 }
 nodeInfo.len = totLenNode;
#ifdef DISTRIBUTED
 nodeInfo.computeOffsets();
#else
 nodeInfo.setMasterFlag();
#endif
 stopTimerMemory(mt.makeInternalInfo, mt.memoryInternal);
}

template<class Scalar>
DistrInfo&
GenDecDomain<Scalar>::ndVecInfo()
{
 // Create nodal Distributed information
 if(!nodeVecInfo) {
   nodeVecInfo = new DistrInfo();
   nodeVecInfo->domLen = new int[numSub];
   nodeVecInfo->numDom = numSub;
   int totLenNode = 0;
   for(int iSub = 0; iSub < numSub; ++iSub) {
     nodeVecInfo->domLen[iSub] = subDomain[iSub]->numNodes(); // this is used for nodal stress output
     totLenNode += nodeVecInfo->domLen[iSub];
   }
   nodeVecInfo->len = totLenNode;
   nodeVecInfo->setMasterFlag();
 }
 return *nodeVecInfo;
}

template<class Scalar>
void
GenDecDomain<Scalar>::constructSubDomains(int iSub)
{
  subDomain[iSub] = GenSubDomainFactory<Scalar>::getFactory()->
    createSubDomain(*domain, iSub, *subToElem, *subToNode, localSubToGl[iSub]);
}


template<class Scalar>
void
GenDecDomain<Scalar>::renumberElements(int iSub)
{
  subDomain[iSub]->renumberElements();
}

template<class Scalar>
void
GenDecDomain<Scalar>::getSharedDOFs()
{
  startTimerMemory(mt.makeInterface, mt.memoryInterface);

  FSCommPattern<int> *nodeIntPat = new FSCommPattern<int>(communicator, cpuToSub, myCPU, FSCommPattern<int>::CopyOnSend);
  for(int i=0; i<numSub; ++i) subDomain[i]->setNodeCommSize(nodeIntPat);
  nodeIntPat->finalize();

  paralApplyToAll(numSub, subDomain, &GenSubDomain<Scalar>::sendDOFList, nodeIntPat);
  nodeIntPat->exchange();
  paralApply(numSub, subDomain, &GenSubDomain<Scalar>::gatherDOFList, nodeIntPat);
  //XXXXif(domain->solInfo().inpc || ((domain->solInfo().newmarkBeta==0.0)&&(domain->solInfo().isDynam()||domain->solInfo().acoustic)))
    paralApply(numSub, subDomain, &GenSubDomain<Scalar>::gatherDOFListPlus, nodeIntPat);
  delete nodeIntPat;
  
  stopTimerMemory(mt.makeInterface, mt.memoryInterface);
/*
  // PJSA: compute the total number of DOFs for the timing file
  int totalCornerDofs = 0;
  for(i=0; i<totalDofs; ++i) if(cornerWeight[i]) totalCornerDofs += (cornerWeight[i] - 1);
  int totalDofs = 0;
  double totalSharedDofs = 0.0;
  for(i=0; i<numSub; ++i) {
    totalDofs += subDomain[i]->getTotalDofCount();
    totalSharedDofs += subDomain[i]->getSharedDofCount();
  }
#ifdef DISTRIBUTED
  totalDofs = communicator->globalSum(totalDofs);
  totalSharedDofs = communicator->globalSum(totalSharedDofs);
#endif
  totalDofs -= (totalCornerDofs + int(floor(totalSharedDofs+0.5)));
  domain->setNumDofs(totalDofs);
*/
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeCorners()
{
  if(!(domain->solInfo().type == 2 && domain->solInfo().fetiInfo.version == FetiInfo::fetidp)) return;
  if(verboseFlag) filePrint(stderr, " ... Selecting the Corners          ...\n");

  FSCommPattern<int> cpat(communicator, cpuToSub, myCPU, FSCommPattern<int>::CopyOnSend);
  for(int i=0; i<numSub; ++i) subDomain[i]->setNodeCommSize(&cpat);
  cpat.finalize();
 
  SubCornerHandler **cornerHandler = new SubCornerHandler * [numSub]; // deleted by cornerMaker
  execParal(numSub, this, &GenDecDomain<Scalar>::makeCornerHandler, cornerHandler);
  CornerMaker cornerMaker(globalNumSub, numSub, cornerHandler, &cpat, communicator);
  cornerMaker.makeCorners();
  grToSub = cornerMaker.getGrToSub();
  execParal(numSub, this, &GenDecDomain<Scalar>::setLocalCorners, cornerHandler);
  
  paralApply(numSub, subDomain, &BaseSub::makeCCDSA);
}


template<class Scalar>
void
GenDecDomain<Scalar>::makeCornerHandler(int iSub, SubCornerHandler **cornerHandler)
{
  cornerHandler[iSub] = subDomain[iSub]->getCornerHandler();
}

template<class Scalar>
void
GenDecDomain<Scalar>::setLocalCorners(int iSub, SubCornerHandler **cornerHandler)
{
  subDomain[iSub]->setCorners(cornerHandler[iSub]->getNumCorners(), cornerHandler[iSub]->getCorners());
}

template<class Scalar>
void GenDecDomain<Scalar>::distributeBCs() 
{
  startTimerMemory(mt.distributeBCs,  mt.memoryDistBC);

  int i, iSub, subI;

  int *nDirichletPerSub = new int[numSub];
  int *nNeumannPerSub   = new int[numSub];
  int *nIDisPerSub      = new int[numSub];
  int *nIDis6PerSub     = new int[numSub];
  int *nIVelPerSub      = new int[numSub];

  // zero all bc counters
  for(iSub = 0; iSub < numSub; ++iSub)  {
    nIDisPerSub[iSub] = 0;
    nIDis6PerSub[iSub] = 0;
    nIVelPerSub[iSub] = 0;
    nDirichletPerSub[iSub] = 0;
    nNeumannPerSub[iSub] = 0;
  }
 
  // get bc's from geoSource
  BCond* dbc = 0;
  BCond* nbc = 0;
  BCond* cvbc = 0;
  BCond* iDis = 0;
  BCond* iDis6 = 0;
  BCond* iVel = 0;

  int numDirichlet = geoSource->getDirichletBC(dbc);
  int numNeuman    = geoSource->getNeumanBC(nbc);
  int numIDis      = geoSource->getIDis(iDis);
  int numIDis6     = geoSource->getIDis6(iDis6);
  int numIVel      = geoSource->getIVel(iVel);

  // Count the number of boundary conditions per subdomain
  int numDispDirichlet = 0; // number of displacement dirichlet BCs
  for(i = 0; i < numDirichlet; ++i) {
    int node = dbc[i].nnum;
    if(dbc[i].dofnum < 6) numDispDirichlet++; 
    for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
        nDirichletPerSub[ subI ]++;
    }
  }
  domain->setNumDispDirichlet(numDispDirichlet);

  for(i = 0; i < numNeuman; ++i) {
    int node = nbc[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub)
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
        nNeumannPerSub[subI]++;
  }

  for(i = 0; i < numIDis; ++i) {
    int node = iDis[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub)
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
        nIDisPerSub[subI]++;
  }

  for(i = 0; i < numIDis6; ++i) {
    int node = iDis6[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub)
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
        nIDis6PerSub[subI]++;
  }

  for(i = 0; i < numIVel; ++i) {
    int node = iVel[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub)
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
        nIVelPerSub[subI]++;
  }

  // Set BC's for all subdomains
  BCond **subBC = new BCond *[numSub];
  for(iSub = 0; iSub < numSub; ++iSub) {
    subBC[iSub] = new BCond[nDirichletPerSub[iSub]];
    nDirichletPerSub[iSub] = 0;
  }
  for(i = 0; i < numDirichlet; ++i) {
    int node = dbc[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
        subBC[subI][nDirichletPerSub[subI]] = dbc[i];
        nDirichletPerSub[subI]++;
      }
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub) {
    subDomain[iSub]->setDirichlet(nDirichletPerSub[iSub], subBC[iSub]);
    subBC[iSub] = new BCond[nNeumannPerSub[iSub]];
    nNeumannPerSub[iSub] = 0;
  }
  for(i = 0; i < numNeuman; ++i) {
    int node = nbc[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
        subBC[subI][nNeumannPerSub[subI]] = nbc[i];
        nNeumannPerSub[subI]++;
      }
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub) {
    subDomain[iSub]->setNeuman(nNeumannPerSub[iSub], subBC[iSub]);
    subBC[iSub] = new BCond[nIDisPerSub[iSub]];
    nIDisPerSub[iSub] = 0;
  }
  for(i = 0; i < numIDis; ++i) {
    int node = iDis[i].nnum;
    for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
        subBC[subI][nIDisPerSub[subI]] = iDis[i];
        nIDisPerSub[subI]++;
      }
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub) {
    subDomain[iSub]->setIDis(nIDisPerSub[iSub], subBC[iSub]);
    subBC[iSub] = new BCond[nIVelPerSub[iSub]];
    nIVelPerSub[iSub] = 0;
  }
  for(i = 0; i < numIVel; ++i) {
    int node = iVel[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
        subBC[subI][nIVelPerSub[subI]] = iVel[i];
        nIVelPerSub[subI]++;
      }
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub) {
    subDomain[iSub]->setIVel(nIVelPerSub[iSub], subBC[iSub]);
    subBC[iSub] = new BCond[nIDis6PerSub[iSub]];
    nIDis6PerSub[iSub] = 0;
  }
  for(i = 0; i < numIDis6; ++i) {
    int node = iDis6[i].nnum;
    for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
        subBC[subI][nIDis6PerSub[subI]] = iDis6[i];
        nIDis6PerSub[subI]++;
      }
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub) {
    subDomain[iSub]->setIDis6(nIDis6PerSub[iSub], subBC[iSub]);
  }

  delete [] nDirichletPerSub;
  delete [] nNeumannPerSub;
  delete [] nIDisPerSub;
  delete [] nIDis6PerSub;
  delete [] nIVelPerSub;
  delete [] subBC;

  if(domain->numSSN() > 0) {

    int totEle = (domain->getEset())->last();
    int *eleTouch = new int[totEle];
    int *eleCount = new int[totEle];

    // Sommerfeld
    for (i=0;i<totEle;i++) eleTouch[i] = -1;
    int *somToSub = new int[domain->numSommer];
    for (i=0;i<domain->numSommer;i++) somToSub[i] = -1;
    for(int iSommEle =0; iSommEle < domain->numSommer; ++iSommEle) {
      int iele = domain->sommer[iSommEle]->findEle(domain->nodeToElem,
                                                   eleTouch, eleCount, iSommEle,
                                                   &domain->getElementSet());
      if(iele < 0) {
        fprintf (stderr, "Error in the Sommerfeld b.c.'s - aborting\n");
        exit(0);
      }
      else {
        subI = glSubToLocal[*(*elemToSub)[iele]];
        if(subI >= 0) somToSub[iSommEle] = subI;
      }
    }

    // Scatterer
    for (i=0;i<totEle;i++) eleTouch[i] = -1;
    int *scaToSub = new int[domain->numScatter];
    for (i=0;i<domain->numScatter;i++) scaToSub[i] = -1;
    int *sBoundFlag = new int[domain->numNode()];
    for(i=0;i<domain->numNode();i++) sBoundFlag[i] = 0;
    for(int iScatter=0; iScatter < domain->numScatter; ++iScatter) {
      for(int iNode = 0;iNode<domain->scatter[iScatter]->numNodes();iNode++) {
        int ndNum = (domain->scatter[iScatter]->getNodes())[iNode];
        sBoundFlag[ndNum] = 1;
      }
      int iele = domain->scatter[iScatter]->findEle(domain->nodeToElem,
                                                    eleTouch, eleCount, iScatter,
                                                    &domain->getElementSet(),1);
      if(iele < 0) {
        fprintf (stderr, "Error in the scatterer b.c.'s - aborting %d \n",iScatter+1);
        exit(0);
      }
      else {
        subI = glSubToLocal[*(*elemToSub)[iele]];
        if(subI >= 0) scaToSub[iScatter] = subI;
      }
    }

   // Wet
   for (i=0;i<totEle;i++) eleTouch[i] = -1;
   int (*wetToSub)[2] = new int[domain->numWet][2];
   for (i=0;i<domain->numWet;i++) wetToSub[i][0] = wetToSub[i][1] = -1;
   int iWetEle;
   for(iWetEle =0; iWetEle < domain->numWet; ++iWetEle) {
     int iele[2];
     domain->wet[iWetEle]->findBothEle(domain->nodeToElem, eleTouch,
        eleCount, iWetEle, &domain->getElementSet(),iele);
     if(iele[0]<0 || iele[1]<0) {
       fprintf (stderr, "Error in the wet b.c.'s - aborting\n");
       exit(0);
     }
     else {
       int subI = glSubToLocal[*(*elemToSub)[iele[0]]];
       if(subI >= 0) wetToSub[iWetEle][0] = subI;
       int subI2 = glSubToLocal[*(*elemToSub)[iele[1]]];
       if(subI2 >= 0) wetToSub[iWetEle][1] = subI2;
     }
   }

    // Implicit Neumann
    for (i=0;i<totEle;i++) eleTouch[i] = -1;
    int *neumToSub = new int[domain->numNeum];
    for (i=0;i<domain->numNeum;i++) neumToSub[i] = -1;
    for(int iNeum=0; iNeum < domain->numNeum; ++iNeum) {
      int iele = domain->neum[iNeum]->findEle(domain->nodeToElem,
                                              eleTouch, eleCount, iNeum);
      if(iele < 0) {
        fprintf (stderr, "Error in the Neumann b.c.'s - aborting\n");
        exit(0);
      }
      else {
        subI = glSubToLocal[*(*elemToSub)[iele]];
        if(subI >= 0) neumToSub[iNeum] = subI;
      }
    }
 
    execParal(numSub, this, &GenDecDomain<Scalar>::distribBC, subDomain, domain,
              somToSub, scaToSub, neumToSub, wetToSub, sBoundFlag);

    delete[] sBoundFlag;
    delete[] somToSub;
    delete[] wetToSub;
    delete[] scaToSub;
    delete[] neumToSub;
    delete[] eleTouch;
    delete[] eleCount;
  }

  // complex nodal boundary conditions 
  if((domain->numComplexDirichlet > 0) || (domain->numComplexNeuman > 0)) {
    int *nComplexDirichletPerSub = new int[numSub]; 
    int *nComplexNeumannPerSub = new int[numSub];
    ComplexBCond **subCBC = new ComplexBCond *[numSub];
    for(iSub = 0; iSub < numSub; ++iSub)
      nComplexDirichletPerSub[iSub] = nComplexNeumannPerSub[iSub] = 0;

    // count Complex Dirichlet
    for(i =0; i < domain->numComplexDirichlet; ++i) {
      int node = domain->cdbc[i].nnum;
      for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
        if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
          nComplexDirichletPerSub[ subI ] += 1;
      }
    }

    // count Complex Neuman
    for(i = 0; i < domain->numComplexNeuman; ++i) {
      int node = domain->cnbc[i].nnum;
      for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
        if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1)
          nComplexNeumannPerSub[ subI ] += 1;
      }
    }

    // distribute Complex Dirichlet
    for(iSub = 0; iSub < numSub; ++iSub) {
      subCBC[iSub] = new ComplexBCond[nComplexDirichletPerSub[iSub]];
      nComplexDirichletPerSub[iSub] = 0;
    }
    for(i =0; i < domain->numComplexDirichlet; ++i) {
      int node = domain->cdbc[i].nnum;
      for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
        if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
          subCBC[subI][ nComplexDirichletPerSub[subI] ] = domain->cdbc[i];
          nComplexDirichletPerSub[subI] += 1;
        }
      }
    }
    for(iSub = 0; iSub < numSub; ++iSub) {
      subDomain[iSub]->setComplexDirichlet(nComplexDirichletPerSub[iSub], subCBC[iSub]);
    }

    // Distribute Complex Neuman
    for(iSub = 0; iSub < numSub; ++iSub) {
      subCBC[iSub] = new ComplexBCond[nComplexNeumannPerSub[iSub]];
      nComplexNeumannPerSub[iSub] = 0;
    }
    for(i =0; i < domain->numComplexNeuman; ++i) {
      int node = domain->cnbc[i].nnum;
      for(iSub =0; iSub < nodeToSub->num(node); ++iSub) {
        if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) {
          subCBC[subI][ nComplexNeumannPerSub[subI] ] = domain->cnbc[i];
          nComplexNeumannPerSub[ subI ] += 1;
        }
      }
    }
    for(iSub = 0; iSub < numSub; ++iSub) {
      subDomain[iSub]->setComplexNeuman(nComplexNeumannPerSub[iSub], subCBC[iSub]);
    }
  
    delete [] nComplexDirichletPerSub;
    delete [] nComplexNeumannPerSub; 
    delete [] subCBC;
  }

  stopTimerMemory(mt.distributeBCs, mt.memoryDistBC);
}

template<class Scalar>
void
GenDecDomain<Scalar>::distribBC(int iSub, GenSubDomain<Scalar> **sd, Domain *domain,
  int *somToSub, int *scaToSub, int *neumToSub, int (*wetToSub)[2], int *sBoundFlag)
{
 int iS;
 for(iS=0;iS<domain->numSommer;iS++) {
   if (somToSub[iS] == iSub) sd[iSub]->addSommer(domain->sommer[iS]);
 }
 for(iS=0;iS<domain->numScatter;iS++)
   if (scaToSub[iS] == iSub) sd[iSub]->addScatter(domain->scatter[iS]);
 for(iS=0;iS<domain->numNeum;iS++)
   if (neumToSub[iS] == iSub) sd[iSub]->addNeum(domain->neum[iS]);
 for(iS=0;iS<domain->numWet;iS++) {
   if (wetToSub[iS][0] == iSub && wetToSub[iS][1] == iSub) {
      sd[iSub]->addWet(domain->wet[iS]);
   }
   else if (wetToSub[iS][0] == iSub && (wetToSub[iS][0] != wetToSub[iS][1])) {
      SommerElement *se = domain->wet[iS]->clone();
      se->el2 = 0;
      sd[iSub]->addWet(se);
   }
   else if (wetToSub[iS][1] == iSub && (wetToSub[iS][0] != wetToSub[iS][1])) {
      SommerElement *se = domain->wet[iS]->clone();
      se->el = se->el2;
      se->el2 = 0;
      sd[iSub]->addWet(se);
   }
 }

 int i;
 for(i=0;i<domain->numNode();i++) {
   if(sBoundFlag[i]) 
     for(int jSub =0; jSub < nodeToSub->num(i); ++jSub) {
       int subJ = glSubToLocal[(*nodeToSub)[i][jSub]];
       if(subJ < 0) continue;
       if(subJ==iSub) sd[iSub]->addSBoundNode(i);
     }
 }

 // RT: Added for SO3
 int c=0;
 sd[iSub]->subScaToSca = new int[sd[iSub]->numScatter];
 for(iS=0;iS<domain->numScatter;iS++)
   if (scaToSub[iS] == iSub) sd[iSub]->subScaToSca[c++] = iS;
}


template<class Scalar>
void GenDecDomain<Scalar>::renumberBC() 
{
  // allocate memory for cluster to local node map
  int *clNodeToLoc = new int[nodeToSub->csize()];

  for(int iSub = 0; iSub < numSub; iSub++)  {
    // initialize node map
    for(int k = 0; k < subDomain[iSub]->numNode(); k++)
      clNodeToLoc[k] = -1;

    // create cluster to local node map
    for(int iNode = 0; iNode < subToNode->num(iSub); iNode++)
      clNodeToLoc[(*subToNode)[iSub][iNode]] = iNode;

    int error = subDomain[iSub]->renumberBC(clNodeToLoc);
    if (error < 0)  {
      fprintf(stderr,"no mapping for subdomain %d\n", iSub);
      exit (-1);
    }
  }
  delete [] clNodeToLoc;
}

template<class Scalar>
void GenDecDomain<Scalar>::setUserDefBC(double *usrDefDisp, double *usrDefVel) 
{
  int iSub;
  for (iSub = 0; iSub < numSub; iSub++)
    subDomain[iSub]->setUserDefBC(usrDefDisp, usrDefVel);
}

template<class Scalar>
void GenDecDomain<Scalar>::addUserForce(GenDistrVector<Scalar> &f, Scalar *userDefineForce) 
{
  int iSub;
  for (iSub = 0; iSub < numSub; iSub++)
    subDomain[iSub]->addUserForce(f.subData(iSub), userDefineForce);
}

template<class Scalar>
void GenDecDomain<Scalar>::addCtrl(GenDistrVector<Scalar> &force, Scalar *ctrfrc) 
{
  int iSub;
  for (iSub = 0; iSub < numSub; iSub++)
    subDomain[iSub]->addCtrl(force.subData(iSub), ctrfrc);
}

template<class Scalar>
void GenDecDomain<Scalar>::extractControlData(GenDistrVector<Scalar> &disp, GenDistrVector<Scalar> &vel,
                                              GenDistrVector<Scalar> &acc, Scalar *ctrdisp,
                                              Scalar *ctrvel, Scalar *ctracc) 
{
  int iSub;
  for (iSub = 0; iSub < numSub; iSub++)
    subDomain[iSub]->extractControlData(disp.subData(iSub),
        vel.subData(iSub), acc.subData(iSub), ctrdisp, ctrvel, ctracc);
}

template<class Scalar>
void GenDecDomain<Scalar>::distributeControlLawData() 
{
  // get Global Control Law
  ControlLawInfo *claw = geoSource->getControlLaw();
  if(!claw) return;

  int i, iSub, subI;

  // allocate Control Law Data
  int *nSensorsPerSub = new int[numSub];
  int *nActuatorsPerSub = new int[numSub];
  int *nUserDispPerSub = new int[numSub];
  int *nUserForcePerSub = new int[numSub];

  // zero all control law counters
  for (iSub = 0; iSub < numSub; ++iSub)  {
    nSensorsPerSub[iSub] = 0;
    nActuatorsPerSub[iSub] = 0;
    nUserDispPerSub[iSub] = 0;
    nUserForcePerSub[iSub] = 0;
  }

  // Count the number of control law data per subdomain
  for(i = 0; i < claw->numSensor; ++i) {
    int node = claw->sensor[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      nSensorsPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numActuator; ++i) {
    int node = claw->actuator[i].nnum;
    for (iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      nActuatorsPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numUserDisp; ++i) {
    int node = claw->userDisp[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      nUserDispPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numUserForce; ++i) {
    int node = claw->userForce[i].nnum;
    for (iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      nUserForcePerSub[subI]++;
    }
  }

  // Set Control Law for all subdomains
  BCond **subSensor = new BCond *[numSub];
  BCond **subActuator = new BCond *[numSub];
  BCond **subUserDisp = new BCond *[numSub];
  BCond **subUserForce = new BCond *[numSub];

  int **locToGlSensorData = new int*[numSub];
  int **locToGlActuatorData = new int*[numSub];
  int **locToGlUserDispData = new int*[numSub];
  int **locToGlUserForceData = new int*[numSub];

  for(iSub = 0; iSub < numSub; ++iSub) {
    subSensor[iSub] = new BCond[nSensorsPerSub[iSub]];
    subActuator[iSub] = new BCond[nActuatorsPerSub[iSub]];
    subUserDisp[iSub] = new BCond[nUserDispPerSub[iSub]];
    subUserForce[iSub] = new BCond[nUserForcePerSub[iSub]];
    locToGlSensorData[iSub] = new int[ nSensorsPerSub[iSub] ];
    locToGlActuatorData[iSub] = new int[ nActuatorsPerSub[iSub] ];
    locToGlUserDispData[iSub] = new int[ nUserDispPerSub[iSub] ];
    locToGlUserForceData[iSub] = new int[ nUserForcePerSub[iSub] ];

    // reinitialize the counters
    nSensorsPerSub[iSub] = 0;
    nActuatorsPerSub[iSub] = 0;
    nUserDispPerSub[iSub] = 0;
    nUserForcePerSub[iSub] = 0;
  }

  for(i = 0; i < claw->numSensor; ++i) {
    int node = claw->sensor[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      subSensor[subI][ nSensorsPerSub[subI] ] = claw->sensor[i];
      locToGlSensorData[subI][ nSensorsPerSub[subI] ] = i;
      nSensorsPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numActuator; ++i) {
    int node = claw->actuator[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      subActuator[subI][ nActuatorsPerSub[subI] ] = claw->actuator[i];
      locToGlActuatorData[subI][ nActuatorsPerSub[subI] ] = i;
      nActuatorsPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numUserDisp; ++i) {
    int node = claw->userDisp[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      subUserDisp[subI][ nUserDispPerSub[subI] ] = claw->userDisp[i];
      locToGlUserDispData[subI][ nUserDispPerSub[subI] ] = i;
      nUserDispPerSub[subI]++;
    }
  }

  for(i = 0; i < claw->numUserForce; ++i) {
    int node = claw->userForce[i].nnum;
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      subI = (*nodeToSub)[node][iSub];
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) == -1) continue;
      subUserForce[subI][ nUserForcePerSub[subI] ] = claw->userForce[i];
      locToGlUserForceData[subI][ nUserForcePerSub[subI] ] = i;
      nUserForcePerSub[subI]++;
    }
  }

  for(iSub = 0; iSub < numSub; ++iSub)  {
    ControlLawInfo *locCLaw = new ControlLawInfo();
    locCLaw->fileName = claw->fileName;
    locCLaw->routineName = claw->routineName;
    locCLaw->numSensor = nSensorsPerSub[iSub];
    locCLaw->sensor = subSensor[iSub];
    locCLaw->numActuator = nActuatorsPerSub[iSub];
    locCLaw->actuator = subActuator[iSub];
    locCLaw->numUserDisp = nUserDispPerSub[iSub];
    locCLaw->userDisp = subUserDisp[iSub];
    locCLaw->numUserForce = nUserForcePerSub[iSub];
    locCLaw->userForce = subUserForce[iSub];

    subDomain[iSub]->setControlData(locCLaw, locToGlSensorData[iSub],
                locToGlActuatorData[iSub], locToGlUserDispData[iSub],
                locToGlUserForceData[iSub]);
  }
  delete [] locToGlSensorData;
  delete [] locToGlActuatorData;
  delete [] locToGlUserDispData;
  delete [] locToGlUserForceData;
  delete [] subUserForce;
  delete [] subUserDisp;
  delete [] subActuator;
  delete [] subSensor;
  delete [] nUserForcePerSub;
  delete [] nUserDispPerSub;
  delete [] nActuatorsPerSub;
  delete [] nSensorsPerSub;
}

// NOTE: this creates global element connectivity in shared
//       and cluster element connectivity in distributed
template<class Scalar>
void GenDecDomain<Scalar>::createElemToNode() 
{
  mt.memoryElemToNode -= memoryUsed();
  int iSub;

  // get total number of elements
  int size = 0;
  for(iSub = 0; iSub < numSub; iSub++)
    size += subDomain[iSub]->numElements();

  // allocate connectivity pointers
  int *ptr = new int[size+1];

  // count number of targets
  int numTargetNodes = 0;
  for(iSub = 0; iSub < numSub; iSub++) {
    Elemset &elems = subDomain[iSub]->getElementSet();
    for(int iElem = 0; iElem < subDomain[iSub]->numElements(); iElem++) {
      numTargetNodes += elems[iElem]->numNodes();
      ptr[ elems[iElem]->getGlNum() ] = elems[iElem]->numNodes();
    }
  }

  // allocate targets
  int *nodeTargets = new int[numTargetNodes];

  // populate targets
  int lastPtr = 0;
  int count = 0;
  for(iSub = 0; iSub < numSub; iSub++) {
    Elemset &elems = subDomain[iSub]->getElementSet();
    for(int iElem = 0; iElem < subDomain[iSub]->numElements(); iElem++) {
      ptr[count++] = lastPtr;
      lastPtr += elems[iElem]->numNodes();  
      elems[iElem]->nodes(nodeTargets + elems[iElem]->getGlNum());
    }
  }
  ptr[size] = lastPtr;
  elemToNode = new Connectivity(size, ptr, nodeTargets);
      
  mt.memoryElemToNode += memoryUsed();
}

template<class Scalar>
void
GenDecDomain<Scalar>::getSharedMPCs()
{
  if(numDualMpc) {
    Connectivity *subToMpc = mpcToSub_dual->altReverse(); // reverse without reordering
    Connectivity *subToSub_mpc = subToMpc->altTranscon(mpcToSub_dual); // remove connection with self unless pure internal
    paralApply(numSub, subDomain, &BaseSub::makeMpcInterface, subToMpc, mpcToSub_dual, subToSub_mpc);
    delete subToMpc;
    delete subToSub_mpc;
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeMpcToMpc()
{
// This function constructs the mpcToMpc connectivity which will be used 
// to assemble the matrix CC^T for the generalized preconditioner (rixen/dual method)
 if(numDualMpc) {

  // build local mpcToMpc connectivities
  paralApplyToAll(numSub, subDomain, &GenSubDomain<Scalar>::makeLocalMpcToMpc);

  int i, j, iSub;
  int size = mpcToSub_dual->csize();
  int *flags = new int[size];
  for(i=0; i<size; ++i) flags[i] = -1;

  // identify the number of connections from MPC i
  int *np = new int[size+1];
  int cp = 0;
  for(i=0; i<size; ++i) {
    np[i] = cp;
    for(iSub=0; iSub < mpcToSub_dual->num(i); ++iSub) {
      int subi = glSubToLocal[(*mpcToSub_dual)[i][iSub]];
      if(subi == -1) continue;
      int ln = subDomain[subi]->globalToLocalMPC[i];
      Connectivity *locMpcToMpc = subDomain[subi]->localMpcToMpc;
      for(j=0; j < locMpcToMpc->num(ln); ++j) {
        int lj = (*locMpcToMpc)[ln][j];
        int gj = subDomain[subi]->localToGlobalMPC[lj];
        if(flags[gj] != i) {
          cp++;
          flags[gj] = i;
        }
      }
    }
  }
  np[size] = cp;

  // Now allocate and fill the new target
  for(i=0; i<size; ++i) flags[i] = -1;
  int *ntg = new int[cp];
  cp = 0;
  for(i=0; i<size; ++i) {
    for(iSub=0; iSub < mpcToSub_dual->num(i); ++iSub) {
      int subi = glSubToLocal[(*mpcToSub_dual)[i][iSub]];
      if(subi == -1) continue;
      int ln = subDomain[subi]->globalToLocalMPC[i];
      Connectivity *locMpcToMpc = subDomain[subi]->localMpcToMpc;
      for(j=0; j < locMpcToMpc->num(ln); ++j) {
        int lj = (*locMpcToMpc)[ln][j];
        int gj = subDomain[subi]->localToGlobalMPC[lj];
        if(flags[gj] != i) {
          ntg[cp] = gj;
          cp++;
          flags[gj] = i;
        }
      }
    }
  }
 
  delete [] flags;
  mpcToMpc = new Connectivity(size, np, ntg); // for all mpcs on this cpu
  if(domain->solInfo().fetiInfo.bmpc) {
    Connectivity *mpcToMpc_mod = mpcToMpc->modify(); 
    delete mpcToMpc;
    mpcToMpc = mpcToMpc_mod;
  }

#ifdef DISTRIBUTED
  makeGlobalMpcToMpc(mpcToMpc); // merge all cpus
#endif

  FetiInfo *finfo = &domain->solInfo().getFetiInfo();
  if(finfo->mpc_precno == FetiInfo::autoSelectCCt) {
    if(mpcToMpc->isDiagonal()) finfo->mpc_precno = FetiInfo::diagCCt;
    else if((numCPU > 1) || (numDualMpc < 100)) finfo->mpc_precno = FetiInfo::globalCCt;
    else finfo->mpc_precno = FetiInfo::superBlockDiagCCt;
  }
 }
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeGlobalMpcToMpc(Connectivity *_procMpcToMpc)
{
// This function constructs the global mpcToMpc connectivity which will be used
// to assemble the matrix CC^T for the generalized preconditioner (rixen/dual method)
  int i, j, k;
  int pid = myCPU;
  int *size = new int[numCPU];
  int *numtarget = new int[numCPU];
  for(i=0; i<numCPU; ++i) {
    if(i == pid) {
      size[i] = _procMpcToMpc->csize();
      numtarget[i] = _procMpcToMpc->numConnect();
    }
    else {
      size[i] = 0;
      numtarget[i] = 0;
    }
  }
#ifdef DISTRIBUTED
  communicator->globalSum(numCPU, size);
  communicator->globalSum(numCPU, numtarget);
#endif

  int totSize = 0;
  int totNumtarget = 0;
  for(i=0; i<numCPU; ++i) {
    totSize += size[i];
    totNumtarget += numtarget[i];
  }

  int *pointer = new int[totSize+numCPU];
  int *target = new int[totNumtarget];
  int startp = 0;
  int startt = 0;
  for(i=0; i<pid; ++i) {
    startp += size[i]+1;
    startt += numtarget[i];
  }
  for(i=0; i<totSize+numCPU; ++i) pointer[i] = 0;
  for(i=0; i<totNumtarget; ++i) target[i] = 0;
  for(i=0; i<size[pid]; ++i)
    pointer[startp + i] = _procMpcToMpc->offset(i);
  pointer[startp + size[pid]] = _procMpcToMpc->numConnect();
  for(i=0; i<numtarget[pid]; ++i)
    target[startt + i] = _procMpcToMpc->getTargetValue(i);
#ifdef DISTRIBUTED
  communicator->globalSum(totSize+numCPU, pointer);
  communicator->globalSum(totNumtarget, target);
#endif

  // now all the _procMpcToMpc connectivity data is stored in size, pointer and target
  Connectivity **tmpMpcToMpc = new Connectivity * [numCPU];
  for(i=0; i<numCPU; ++i) {
    startp = 0; startt = 0;
    for(j=0; j<i; ++j) {
      startp += size[j]+1;
      startt += numtarget[j];
    }
    tmpMpcToMpc[i] = new Connectivity(size[i], pointer+startp, target+startt);
  }
  // now each processor has the _procMpcToMpc connectivities for all other processors
  Connectivity *subToCpu = cpuToSub->reverse();
  mpcToCpu = mpcToSub_dual->transcon(subToCpu);
  delete subToCpu; subToCpu = 0;
  int csize = mpcToCpu->csize();
  int *flags = new int[csize];
  for(i=0; i<csize; ++i) flags[i] = -1;
  // identify the number of connections from MPC i
  int *np = new int[csize+1];
  int cp = 0;
  for(i=0; i<csize; ++i) {
    np[i] = cp;
    for(j=0; j < mpcToCpu->num(i); ++j) {
      int cpu = (*mpcToCpu)[i][j];
      for(k=0; k < tmpMpcToMpc[cpu]->num(i); ++k) {
        int mpck = (*tmpMpcToMpc[cpu])[i][k];
        if(flags[mpck] != i) {
          cp++;
          flags[mpck] = i;
        }
      }
    }
  }
  np[csize] = cp;

  // Now allocate and fill the new target
  for(i=0; i<csize; ++i) flags[i] = -1;
  int *ntg = new int[cp];
  cp = 0;
  for(i=0; i<csize; ++i) {
    np[i] = cp;
    for(j=0; j < mpcToCpu->num(i); ++j) {
      int cpu = (*mpcToCpu)[i][j];
      for(k=0; k < tmpMpcToMpc[cpu]->num(i); ++k) {
        int mpck = (*tmpMpcToMpc[cpu])[i][k];
        if(flags[mpck] != i) {
          ntg[cp] = mpck;
          cp++;
          flags[mpck] = i;
        }
      }
    }
  }

  delete [] size;
  delete [] numtarget;
  delete [] flags;
  delete [] pointer;
  delete [] target;
  delete [] tmpMpcToMpc;
#ifdef USE_MUMPS
  if(domain->solInfo().fetiInfo.cctSolver == FetiInfo::mumps && domain->solInfo().mumps_icntl[18] == 3) {
    procMpcToMpc = _procMpcToMpc;
  } else
#endif
  delete _procMpcToMpc;
  mpcToMpc = new Connectivity(csize, np, ntg);
}

template<class Scalar>
void GenDecDomain<Scalar>::addFsiElements()
{
 if ( (domain->solInfo().isCoupled) && (domain->solInfo().isMatching) &&
      (domain->solInfo().fetiInfo.fsi_corner != 0) ) {
// JLchange: replace addSubFsiElem() such that fsi elements are added only to structure elements. 
   for (int i=0; i< domain->getNumFSI(); ++i) {
     LMPCons *thisGlFSI = domain->getFsi(i);
     int glFluidNode = thisGlFSI->lmpcnum;
     for(int j=0; j< thisGlFSI->nterms; j++) {
       int glStrucNode = thisGlFSI->terms[j].nnum;
       // Find all subdomains which have a structure element covering these two nodes 

       int weight = 0;
       for (int fSub=0; fSub< nodeToSub->num(glFluidNode); ++fSub) {
         int thisSubNum = (*nodeToSub)[glFluidNode][fSub];
         if (nodeToSub->locate(glStrucNode,thisSubNum)) {
           bool foundOne = false;
           for (int iElem=0; iElem< (domain->nodeToElem)->num(glFluidNode); ++iElem) {
             if (foundOne == true) break;
             int thisElemNum = (*(domain->nodeToElem))[glFluidNode][iElem];
             if ((domain->nodeToElem)->locate(glStrucNode,thisElemNum) &&
                 domain->isStructureElement(thisElemNum) && subToElem->locate(thisSubNum,thisElemNum)) {
               weight++;
               foundOne = true;
             }
           }
         }
       }

       for (int fSub=0; fSub< nodeToSub->num(glFluidNode); ++fSub) {
         int thisSubNum = (*nodeToSub)[glFluidNode][fSub];
         int localSubNumber = glSubToLocal[thisSubNum];
         if ((localSubNumber >= 0) && nodeToSub->locate(glStrucNode,thisSubNum)) {
           bool foundOne = false;
           for (int iElem=0; iElem< (domain->nodeToElem)->num(glFluidNode); ++iElem) {
             if (foundOne == true) break;
             int thisElemNum = (*(domain->nodeToElem))[glFluidNode][iElem];
             if ((domain->nodeToElem)->locate(glStrucNode,thisElemNum) &&
                 domain->isStructureElement(thisElemNum) && subToElem->locate(thisSubNum,thisElemNum)) {
               LMPCons *localFsi = new LMPCons(glFluidNode, 0.0);
               LMPCTerm thisLmpcTerm(thisGlFSI->terms[j], 1.0/weight);
               localFsi->addterm(&(thisLmpcTerm));
               subDomain[localSubNumber]->addSingleFsi(localFsi);
               foundOne = true;
             }
           }
         }
       }
     }
   }
 }
}

template<class Scalar>
void GenDecDomain<Scalar>::preProcessFSIs()
{
  if(!domain->solInfo().isCoupled) return;
  if(verboseFlag) filePrint(stderr, " ... Applying the Fluid-Structure Interactions");
  domain->computeCoupledScaleFactors();
  domain->makeFsiToNode();
  wetInterfaceNodes = domain->getAllWetInterfaceNodes(numWetInterfaceNodes);
  distributeWetInterfaceNodes();
  if(verboseFlag) filePrint(stderr, " ...\n");
}

template<class Scalar>
void GenDecDomain<Scalar>::distributeWetInterfaceNodes()
{
  int i, iSub, subI;
  int *nWetInterfaceNodesPerSub = new int[numSub];
  for (iSub = 0; iSub < numSub; ++iSub)  {
    nWetInterfaceNodesPerSub[iSub] = 0;
  }
  for(i = 0; i < numWetInterfaceNodes; ++i) {
    int node = wetInterfaceNodes[i];
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub){
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) > -1) 
        nWetInterfaceNodesPerSub[subI]++;
    }
  }

  int **subWetInterfaceNodes = new int *[numSub];
  for(iSub = 0; iSub < numSub; ++iSub) {
    subWetInterfaceNodes[iSub] = new int[nWetInterfaceNodesPerSub[iSub]];
    nWetInterfaceNodesPerSub[iSub] = 0;
  }
  for(i = 0; i < numWetInterfaceNodes; ++i) {
    int node = wetInterfaceNodes[i];
    for(iSub = 0; iSub < nodeToSub->num(node); ++iSub) {
      if((subI = glSubToLocal[(*nodeToSub)[node][iSub]]) < 0) continue;
      subWetInterfaceNodes[subI][ nWetInterfaceNodesPerSub[subI] ] = node;  // global number
      nWetInterfaceNodesPerSub[subI]++;
    }
  }

  if ((domain->solInfo().isCoupled) && (domain->solInfo().fetiInfo.fsi_corner != 0)) { 
    execParal(numSub, this, &GenDecDomain<Scalar>::markSubWetInterface,
              nWetInterfaceNodesPerSub, subWetInterfaceNodes);
    //makeCorners();  // XXXX JLchange
  }  

  execParal(numSub, this, &GenDecDomain<Scalar>::setSubWetInterface, nWetInterfaceNodesPerSub, subWetInterfaceNodes);
  delete [] nWetInterfaceNodesPerSub;
  for(iSub = 0; iSub < numSub; ++iSub) delete [] subWetInterfaceNodes[iSub];
  delete [] subWetInterfaceNodes;
}

template<class Scalar>
void GenDecDomain<Scalar>::markSubWetInterface(int iSub, int *nWetInterfaceNodesPerSub, int **subWetInterfaceNodes)
{
  subDomain[iSub]->markWetInterface(nWetInterfaceNodesPerSub[iSub], subWetInterfaceNodes[iSub]);
}

template<class Scalar>
void GenDecDomain<Scalar>::setSubWetInterface(int iSub, int *nWetInterfaceNodesPerSub, int **subWetInterfaceNodes)
{
  subDomain[iSub]->setWetInterface(nWetInterfaceNodesPerSub[iSub], subWetInterfaceNodes[iSub]);
}

template<class Scalar>
void GenDecDomain<Scalar>::getSharedFSIs()
{
  if(domain->solInfo().isCoupled) {
    Connectivity *fsiToSub = domain->getFsiToNode()->transcon(nodeToSub);
    Connectivity *subToFsi = fsiToSub->altReverse(); // reverse without reordering
    Connectivity *subToSub_fsi = subToFsi->transcon(fsiToSub);
    paralApply(numSub, subDomain, &BaseSub::makeFsiInterface, subToFsi, fsiToSub, subToSub_fsi);
    delete fsiToSub;
    delete subToFsi;
    delete subToSub_fsi;
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::makeMpcToSub()
{
  if(numDualMpc) {
    int *pointer = new int[globalNumSub+1];
    for(int i =0; i < globalNumSub; ++i) pointer[i] = 0;
    paralApply(numSub, subDomain, &BaseSub::putNumMPC, pointer);
#ifdef DISTRIBUTED
    communicator->globalSum(globalNumSub, pointer);
#endif

    int total = 0;
    for(int i = 0; i < globalNumSub; ++i) {
      int tmp = pointer[i];
      pointer[i] = total;
      total += tmp;
    }
    pointer[globalNumSub] = total;

    int *target = new int[total];
    for(int i = 0; i < total; ++i) target[i] = 0;
    paralApply(numSub, subDomain, &BaseSub::putLocalToGlobalMPC, pointer, target);
#ifdef DISTRIBUTED
    communicator->globalSum(total, target);
#endif

    Connectivity *subToMpc = new Connectivity(globalNumSub, pointer, target);
    if(mpcToSub_dual) delete mpcToSub_dual;
    mpcToSub_dual = subToMpc->reverse();
    delete subToMpc;
  }
  if(numPrimalMpc) {
    int *pointer = new int[globalNumSub+1];
    for(int i =0; i < globalNumSub; ++i) pointer[i] = 0;
    paralApply(numSub, subDomain, &BaseSub::putNumMPC_primal, pointer);
#ifdef DISTRIBUTED
    communicator->globalSum(globalNumSub, pointer);
#endif

    int total = 0;
    for(int i = 0; i < globalNumSub; ++i) {
      int tmp = pointer[i];
      pointer[i] = total;
      total += tmp;
    }
    pointer[globalNumSub] = total;

    int *target = new int[total];
    for(int i = 0; i < total; ++i) target[i] = 0;
    paralApply(numSub, subDomain, &BaseSub::putLocalToGlobalMPC_primal, pointer, target);
#ifdef DISTRIBUTED
    communicator->globalSum(total, target);
#endif

    Connectivity *subToMpc = new Connectivity(globalNumSub, pointer, target);
    if(mpcToSub_primal) delete mpcToSub_primal;
    mpcToSub_primal = subToMpc->reverse();
    delete subToMpc;
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::buildOps(GenMDDynamMat<Scalar> &res, double coeM, double coeC, double coeK, Rbm **rbms, FullSquareMatrix **kelArray, bool make_feti)
{
 GenBasicAssembler<Scalar> *ba = 0;

 bool isFeti = domain->solInfo().type == 2;
 FetiInfo *finfo = &domain->solInfo().getFetiInfo();

 int solvertype = finfo->solvertype;

 int isFeti2 =  (isFeti && finfo->version == FetiInfo::feti2) ? 1 : 0;
 int isCtcOrDualMpc = (numDualMpc) ? 1 : 0;

 GenDomainGroupTask<Scalar> dgt(numSub, subDomain, coeM, coeC, coeK, rbms, kelArray,
                                domain->solInfo().alphaDamp, domain->solInfo().betaDamp, domain->numSommer, isFeti2,
                                solvertype, isCtcOrDualMpc);

 if(verboseFlag) filePrint(stderr," ... Assemble Subdomain Matrices    ... \n");
 execParal(numSub, &dgt, &GenDomainGroupTask<Scalar>::runFor, make_feti);

 if(domain->solInfo().inpc) {
   FSCommPattern<Scalar> *pat = new FSCommPattern<Scalar>(communicator, cpuToSub, myCPU, FSCommPattern<Scalar>::CopyOnSend);
   for(int i=0; i<numSub; ++i) subDomain[i]->setDofPlusCommSize(pat);
   pat->finalize();
   ba = new GenBasicAssembler<Scalar>(numSub, subDomain, pat);
 }

 res.K   = new GenSubDOp<Scalar>(numSub, dgt.K, ba);
 res.Kuc = new GenSubDOp<Scalar>(numSub, dgt.Kuc);

 if(dgt.C[0]) {
   res.C = new GenSubDOp<Scalar>(numSub, dgt.C);
   res.Cuc = new GenSubDOp<Scalar>(numSub, dgt.Cuc);
 }
 else {
   res.C   = 0; delete [] dgt.C;
   res.Cuc = 0; delete [] dgt.Cuc;
 }
 res.M   = new GenSubDOp<Scalar>(numSub, dgt.M);
 res.Muc = new GenSubDOp<Scalar>(numSub, dgt.Muc);

// RT
 if(dgt.C_deriv[0]) {
   res.C_deriv = new GenSubDOp<Scalar>*[1];
   (res.C_deriv)[0] = new GenSubDOp<Scalar>(numSub, dgt.C_deriv,0);
 } else {
   res.C_deriv   = 0; //delete [] dgt.C_deriv;
 }
 if(dgt.Cuc_deriv[0]) {
   res.Cuc_deriv = new GenSubDOp<Scalar>*[1];
   res.Cuc_deriv[0] = new GenSubDOp<Scalar>(numSub, dgt.Cuc_deriv,0);
 } else {
   res.Cuc_deriv = 0; //delete [] dgt.Cuc_deriv;
 }
// RT end

 if(isFeti) {
   if(make_feti) res.dynMat = getDynamicFetiSolver(dgt);
 } else
   res.dynMat = getDiagSolver(numSub, dgt.sd, dgt.dynMats);
}

template<class Scalar>
void
GenDecDomain<Scalar>::rebuildOps(GenMDDynamMat<Scalar> &res, double coeM, double coeC, double coeK, 
                                 FullSquareMatrix **kelArray, FullSquareMatrix **melArray)
{
 res.dynMat->reconstruct(); // do anything that needs to be done before zeroing and assembling the matrices

 execParal6R(numSub, this, &GenDecDomain<Scalar>::subRebuildOps, res, coeM, coeC, coeK, kelArray, melArray);

 res.dynMat->refactor(); // do anything that needs to be done after zeroing and assembling the matrices
}

template<class Scalar>
void
GenDecDomain<Scalar>::subRebuildOps(int iSub, GenMDDynamMat<Scalar> &res, double coeM, double coeC, double coeK, 
                                    FullSquareMatrix **kelArray, FullSquareMatrix **melArray)
{
  AllOps<Scalar> allOps;

  if(res.K)  allOps.K = (*res.K)[iSub];
  if(res.C)  allOps.C = (*res.C)[iSub];
  if(res.Cuc)  allOps.Cuc = (*res.Cuc)[iSub];
  if(res.M)  allOps.M = (*res.M)[iSub];
  if(res.Muc)  allOps.Muc = (*res.Muc)[iSub];
  if(res.Mcc)  allOps.Mcc = (*res.Mcc)[iSub];
  if(res.Kuc)  allOps.Kuc = (*res.Kuc)[iSub];

  allOps.zero();

  GenMultiSparse<Scalar> allMats(subDomain[iSub]->KrrSparse, subDomain[iSub]->KiiSparse, subDomain[iSub]->Kbb,
                                 subDomain[iSub]->Kib, subDomain[iSub]->Krc, subDomain[iSub]->Kcc);
  allMats.zeroAll();

  subDomain[iSub]->template makeSparseOps<Scalar>(allOps, coeK, coeM, coeC, &allMats, (kelArray) ? kelArray[iSub] : 0, 
                                                  (melArray) ? melArray[iSub] : 0);
}

template<class Scalar>
void
GenDecDomain<Scalar>::getWError(int iSub, GenDistrVector<Scalar> *u,
                                double *l2err, double *h1err, double *l2, double *h1)
{
 // subDomain[iSub]->wError( (BaseSub *) subDomain[iSub], l2err+iSub, h1err+iSub,
 //                         l2+iSub, h1+iSub, u->subData(iSub));
}

template<>
void
GenDecDomain<double>::buildFFP(GenDistrVector<double> &u, FILE *fffp);

template<>
void
GenDecDomain<DComplex>::buildFFP(GenDistrVector<DComplex> &u, FILE *fffp);

template<>
void
GenDecDomain<DComplex>::buildLocalFFP(int iSub, GenDistrVector<DComplex> *u,
                                      DComplex **ffp, int *numSample, double (*dir)[3]);

template<>
void
GenDecDomain<double>::buildLocalFFP(int iSub, GenDistrVector<double> *u,
                                    double **ffp, int *numSample, double (*dir)[3]);

template<class Scalar>
void
GenDecDomain<Scalar>::setNewProperties(int i)
{
  paralApply(numSub, subDomain, &Domain::setNewProperties, i);
}

template<class Scalar>
void
GenDecDomain<Scalar>::assignRandMat()
{
  paralApply(numSub, subDomain, &Domain::assignRandMat);
}

template<class Scalar>
void
GenDecDomain<Scalar>::retrieveElemset()
{
  paralApply(numSub, subDomain, &Domain::retrieveElemset);
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDecDomain<Scalar>::getElementAttr(int fileNumber,int iAttr, double time)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;
  int numNodes = geoSource->numNode();

  double *weights = (double *) dbg_alloca(numNodes*sizeof(double));
  double *props = (double *) dbg_alloca(numNodes*sizeof(double));
  for(int i=0; i<numNodes; ++i) weights[i] = props[i] = 0.0;

  // ... WRITE THE TIME VALUE
  //filePrint(oinfo[fileNumber].filptr,"%20.10e\n",domain->solInfo().getTimeStep());
  // ... OUTPUT precision
  //int p = oinfo[fileNumber].precision;

  for(int i=0; i<this->getNumSub(); ++i)
    { this->getSubDomain(i)->mergeElemProps(props, weights, iAttr); }

  if(avgnum == 1 ) 
    {
      for(int k=0; k<numNodes; ++k)  
	{
	  if(weights[k] > 0)
	    { props[k] /= weights[k]; }
	}      
      geoSource->outputNodeScalars(fileNumber, props, numNodes, time);
    }
  else
    {
      // not implemented
      assert(0);
    }
  return;
}

template<class Scalar>
void GenDecDomain<Scalar>::setContactGap(DistrGeomState *geomState, GenFetiSolver<Scalar> *fetiSolver)
{
  // XXXX this should only be done for linear mpcs
  if(numDualMpc) {
    GenDistrVector<Scalar> *cx = new GenDistrVector<Scalar>(fetiSolver->interfInfo());
    GenDistrVector<Scalar> *x = new GenDistrVector<Scalar>(fetiSolver->localInfo());
    execParal2R(numSub, this, &GenDecDomain<Scalar>::extractPosition, *geomState, *x);
    ((GenFetiDPSolver<Scalar> *)fetiSolver)->multC(*x, *cx); // cx = C*x
    (*cx) *= -1.0;
    execParal(this->numSub, this, &GenDecDomain<Scalar>::setMpcRhs, *cx);
    delete cx; delete x;
  }
}

template<class Scalar>
void
GenDecDomain<Scalar>::extractPosition(int iSub, DistrGeomState &geomState, GenDistrVector<Scalar> &x)
{
  geomState[iSub]->extract(x.subData(subDomain[iSub]->localSubNum()));
}

template<class Scalar>
void
GenDecDomain<Scalar>::setMpcRhs(int iSub, GenDistrVector<Scalar> &cu)
{
  subDomain[iSub]->setMpcRhs(cu.subData(subDomain[iSub]->localSubNum()));
}

template<class Scalar>
void
GenDecDomain<Scalar>::printLMPC()
{
  for(int i=0; i<numSub; ++i) subDomain[i]->printLMPC();
}
