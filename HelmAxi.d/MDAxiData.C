#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/SubDomain.h>
#include <Driver.d/PolygonSet.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/SparseSet.h>
#include <Math.d/matrix.h>
#include <Threads.d/Paral.h>
#include <Timers.d/Timing.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <HelmAxi.d/AxiHElem.h>
#include <HelmAxi.d/FourierHelmBCs.h>
#include <HelmAxi.d/MPCData.h>
#include <HelmAxi.d/MDAxiData.h>

#if defined(sgi) && ! defined(_OPENMP)
#include <ulocks.h>
extern ulock_t allocLock;
#endif

extern long totMemSky;
extern long totMemSparse;
extern ThreadManager *threadManager;


/////////////////////////////////////
// This is included for debugging
#include <Timers.d/GetTime.h>
/////////////////////////////////////

MDAxiData::MDAxiData(Domain *d, int n, Connectivity *cn, Connectivity *sToN) :
 SubDomain(*d, n, *cn, *sToN, n) {

  locBCs = new FourierHelmBCs();

  hParams = new GenHelmParam;
  hParams->kappaHelm = geoSource->kappa();

  KC   = 0;
  KucC = 0;
  KbbC = 0;

  KibC = 0;
  KiiC = 0;
  KccC = 0;

  numInterface   = 0;
  InterfElem     = 0;
  InterfaceIndex = 0;

  locToglDir = 0;

  InterfSign = 0 ;

  coarseEqNums = 0;

  numInQSet = 1;               // default value for a simple FETI-H-AXI

  locMPCs = 0;
  CKQ  = 0;
  CKCt = 0;
  localToGlobalMPC = 0;
  globalToLocalMPC = 0;

}


void
MDAxiData::extractBCs(FourierHelmBCs *glBCs, Connectivity *nodeToSub) {

// UH - 9/9/99
//
// For a node on the interface with Neumann BC, the force will be
// distributed as full in THE subdomain containing the element and
// zero for all the other subdomains.
// The same idea is applied for Sommerfeld.
//

 int i;

 locBCs->setModes(glBCs->numModes);
 locBCs->setSlices(glBCs->numSlices);
 locBCs->setDir(glBCs->dirX, glBCs->dirY, glBCs->dirZ);
 locBCs->setConst(glBCs->cstant);

 for(i =0; i < glBCs->numDir; ++i) {
   int node = glBCs->dirBC[i].nnum;
   if (nodeToSub->offset(node,subNumber) != -1)
      locBCs->addDirichlet(glBCs->dirBC[i]);
 }

 if (locBCs->numDir) {
   locToglDir = new int[locBCs->numDir];
   int offSet = 0;
   for (i=0; i<glBCs->numDir; ++i) {
     int node = glBCs->dirBC[i].nnum;
     if (nodeToSub->offset(node,subNumber) != -1)
       locToglDir[offSet++] = i;
   }
 }

 for(i = 0; i < glBCs->numNeu; ++i) {
   int node1 = glBCs->neuBC[i].nnum1;
   int node2 = glBCs->neuBC[i].nnum2;

   if ((nodeToSub->offset(node1,subNumber) == -1) ||
       (nodeToSub->offset(node2,subNumber) == -1))
       continue;

   locBCs->addNeuman(glBCs->neuBC[i]);
 }

 if (glBCs->numSomm) {

   locBCs->setSomType(glBCs->somType);
   locBCs->setSurf(glBCs->surfR, glBCs->surfZ);

   int *node = 0;
   for (i = 0; i < glBCs->numSomm; ++i) {

     node = glBCs->somBC[i]->nodes(node);

     if ((nodeToSub->offset(node[0],subNumber) == -1) ||
        (nodeToSub->offset(node[1],subNumber) == -1))
         continue;

     locBCs->addSommer(glBCs->somBC[i]);

   }
 }

}


void
MDAxiData::extractMPCs(MPCData *glMPCs, Connectivity *mpcToNode,
           Connectivity *nodeToSub) {

 numInQSet += 1;

 int i, j;
 int numMPC = 0;
 int glNumMPC = (mpcToNode) ? mpcToNode->csize() : 0;
 int *activeMPC = new int[glNumMPC];

 // PASS 1 : count number of local mpcs
 // Note : We distribute the MPC equally to each subdomain for
 //        a shared node when we build the matrix

 for (i=0; i<glNumMPC; ++i) {
   // Note : One MPC involves one node but several "Fourier-dofs" of
   //        this particular node
   if (nodeToSub->offset((*mpcToNode)[i][0],subNumber) != -1) {
     activeMPC[i] = 1;
     numMPC++;
   }
   else
     activeMPC[i] = -1;
 }

 if (numMPC == 0)
    return;

 locMPCs = new MPCData();
 locMPCs->totalNumMPC = numMPC;

 localToGlobalMPC = new int[numMPC];
 globalToLocalMPC = new int[glNumMPC];
 for (i=0; i<glNumMPC; ++i)
    globalToLocalMPC[i] = -1;

 int offSetLoc = 0;
 int offSetGl = 0;

 // PASS 2 : Get the mpc values

 for (i=0; i<glMPCs->numMPCSet; ++i) {
   if (activeMPC[offSetGl]==1) {
      locMPCs->addMPC(glMPCs->term[i]);
      for (j=0; j<(glMPCs->term[i].slices + 1); ++j) {
          localToGlobalMPC[offSetLoc]= offSetGl;
          globalToLocalMPC[offSetGl] = offSetLoc;
          offSetGl++;
          offSetLoc++;
      }
   }
   else
      offSetGl += glMPCs->term[i].slices + 1;
 }

 delete[] activeMPC;
 activeMPC = 0;

}


int*
MDAxiData::getInterfaceIndex() {

  return InterfaceIndex;

}


void
MDAxiData::renumberData() {

 int globalNMax = 0;
 int i;

 for (i = 0; i < numnodes; ++i) {
   if(glNums[i] > globalNMax) globalNMax = glNums[i];
 }

 globalNMax += 1;
 int *glToLocalNode = new int[globalNMax];

 for (i=0; i < numnodes; ++i)
   glToLocalNode[glNums[i]] = i;

 int nc = scomm->sharedNodes->numConnect();
 int *allC = (*scomm->sharedNodes)[0];

 for(i = 0; i < nc; ++i)
   allC[i] = glToLocalNode[allC[i]];

 for(i=0; i < numele; ++i)
   packedEset[i]->renum(glToLocalNode);

 for(i=0; i<locBCs->numDir; ++i)
   locBCs->dirBC[i].nnum = glToLocalNode[locBCs->dirBC[i].nnum];

 for(i=0; i<locBCs->numNeu; ++i) {
   locBCs->neuBC[i].nnum1 = glToLocalNode[locBCs->neuBC[i].nnum1];
   locBCs->neuBC[i].nnum2 = glToLocalNode[locBCs->neuBC[i].nnum2];
   if (locBCs->neuBC[i].nnum3>-1)
     locBCs->neuBC[i].nnum3 = glToLocalNode[locBCs->neuBC[i].nnum3];
 }

 for(i=0; i < locBCs->numSomm; ++i)
   locBCs->somBC[i]->renum(glToLocalNode);

 if (locMPCs)
    for (i=0; i<locMPCs->numMPCSet; ++i) {
      locMPCs->term[i].nnum = glToLocalNode[locMPCs->term[i].nnum];
    }

 delete[] glToLocalNode;
 glToLocalNode = 0;

}


void
MDAxiData::getInterface(PolygonSet ***allPs) {

 // Set up interface elements

 int *map = new int[numnodes];

 PolygonSet **ps = new PolygonSet *[scomm->sharedNodes->csize()];

 allPs[subNum()] = ps;

 int i;
 for (i=0; i<numnodes; i++)
   map[i] = -1;

 if (elemToNode == 0)
      elemToNode = new Connectivity(&packedEset);

 if (nodeToElem == 0)
      nodeToElem = elemToNode->reverse();

 Connectivity *edgeToElem = scomm->sharedNodes->transcon(nodeToElem);

 int iInter;
 for (iInter=0;iInter<scomm->sharedNodes->csize(); iInter++) {

   for(i = 0; i < scomm->sharedNodes->num(iInter); ++i)
      map[(*scomm->sharedNodes)[iInter][i]] = i;

   ps[iInter] = new PolygonSet(map);

   for (int iele=0; iele<edgeToElem->num(iInter); iele++)
     packedEset[(*edgeToElem)[iInter][iele]]->addFaces(ps[iInter]);

   for (i = 0; i < scomm->sharedNodes->num(iInter); ++i)
      map[(*scomm->sharedNodes)[iInter][i]] = -1;

   setSendData(iInter, ps[iInter]);
 }

 delete[] map;
 map = 0;

}


void
MDAxiData::finishInterface(PolygonSet ***allPs) {

  PolygonSet **ps = allPs[subNum()];

  int iInter;

  int numNeighb = scomm->sharedNodes->csize();

  InterfaceIndex = new int[numNeighb+1];

  int *iSNum = new int[scomm->sharedNodes->csize()];
  for (iInter=0;iInter<scomm->sharedNodes->csize(); iInter++)
    iSNum[iInter] = ps[iInter]->getSize();

  numInterface = 0;
  InterfaceIndex[0] = 0;

  for (iInter=0;iInter<scomm->sharedNodes->csize(); iInter++) {
    PolygonSet *neighbPS = (PolygonSet *) scomm->getExchangeData(iInter);
    numInterface += (iSNum[iInter] = ps[iInter]->selfIntersect(neighbPS));
    InterfaceIndex[iInter+1] = numInterface;
  }

  InterfElem = new SommerElement *[numInterface];
  numInterface = 0;

  for (iInter=0;iInter<scomm->sharedNodes->csize(); iInter++) {
    ps[iInter]->getAxiSommerElem((*scomm->sharedNodes)[iInter],
                                InterfElem+numInterface);
    numInterface += iSNum[iInter];
  }

 delete[] iSNum;
 iSNum = 0;
}

void
MDAxiData::sendDOFList()
{
 Connectivity &sharedNodes = *(scomm->sharedNodes);

 // First make our own dsa and c_dsa
 if (elemToNode==0)
   elemToNode = new Connectivity(&packedEset);

 // This routines creates the connectivities nodeToElem, nodeToNode
 Renumber rnum = getRenumbering();

 dsa = new DofSetArray(numnodes, packedEset, rnum.renumb);

 //int numdofs = dsa->size();

 BCond *buffer;
 int i;

 buffer = new BCond[locBCs->numDir];
 for (i=0; i<locBCs->numDir; ++i) {
    buffer[i].nnum = locBCs->dirBC[i].nnum;
    buffer[i].dofnum = 7;
    buffer[i].val = 0.0;
 }

 c_dsa = new ConstrainedDSA( *dsa, locBCs->numDir, buffer);

 int iSub, iNode;

 for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    DofSet *interDOFs = new DofSet[sharedNodes.num(iSub)];
    for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode)
      interDOFs[iNode] = (*c_dsa)[ sharedNodes[iSub][iNode] ];

    setSendData(iSub , (void *) interDOFs);
 }

 delete[] buffer;
 buffer = 0;
}

void
MDAxiData::gatherDOFList(DofSet ***allBoundary)
{
  int iSub, iNode;
  Connectivity &sharedNodes = *(scomm->sharedNodes);
  DofSet **boundaryDOFs = new DofSet *[scomm->numNeighb];
  allBoundary[subNum()] = boundaryDOFs;
  int nbdofs = 0;

  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    boundaryDOFs[iSub] = (DofSet *) scomm->getExchangeData(iSub);
    for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
      boundaryDOFs[iSub][iNode] = boundaryDOFs[iSub][iNode] &
                              (*c_dsa)[ sharedNodes[iSub][iNode] ];
      nbdofs += boundaryDOFs[iSub][iNode].count();
    }
  }

  int *allBoundDofs = new int[nbdofs];
  int *boundDofPointer = new int[scomm->numNeighb+1];
  nbdofs = 0;

  for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
     boundDofPointer[iSub] = nbdofs;
     for(iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
        int nnn = (*c_dsa).number(sharedNodes[iSub][iNode],
                    boundaryDOFs[iSub][iNode],  allBoundDofs+nbdofs);
        nbdofs += boundaryDOFs[iSub][iNode].count();
        if(boundaryDOFs[iSub][iNode].count() != nnn)
          fprintf(stderr,"Strange: %d %d\n",
                  boundaryDOFs[iSub][iNode].count(),nnn);
     }
  }

  boundDofPointer[scomm->numNeighb] = nbdofs;
  scomm->sharedDOFs =
      new Connectivity(scomm->numNeighb, boundDofPointer, allBoundDofs);

// Now count the number of subdomains touching one DOF
  int iDof;
  int ndof = c_dsa->size();
  weight = new int[ndof];
  for (iDof = 0; iDof < ndof; ++iDof)
     weight[iDof] = 1;
  for (iDof = 0; iDof < nbdofs; ++iDof)
     weight[allBoundDofs[iDof]] += 1;

// nbdofs = number of boundary dofs
  scaling = new double[nbdofs];

  for(iDof = 0; iDof < nbdofs; ++iDof)
       scaling[iDof] = 1.0/weight[allBoundDofs[iDof]];

}


void
MDAxiData::defineKs(int **glBoundMap, int **glInternalMap) {

 KC = new ComplexSolver *[locBCs->numModes+1];
 KbbC = new DBComplexSparseMatrix *[locBCs->numModes+1];

 if (solInfo().getFetiInfo().precno == FetiInfo::dirichlet) {
   if (internalLen>0)
     KiiC = new ComplexSolver *[locBCs->numModes+1];
 }
 else
   KiiC = 0;

 if (locBCs->numDir>0) {
   KucC = new CuCComplexSparse *[locBCs->numModes+1];
   KccC = new DBComplexSparseMatrix *[locBCs->numModes+1];
 }

 makeAllDOFs();

 glBoundMap[subNumber] = makeBMaps(c_dsa);
 if (KiiC)
   glInternalMap[subNumber] = makeIMaps(c_dsa);
 else
   glInternalMap[subNumber] = 0;

}


void
MDAxiData::makeKs(int *glBoundMap, int *glInternalMap, int m1, int m2) {

// Memory allocation of the different matrices
// The size of the arrays is locBCs->numModes+1 for all
// the basis functions 1,cos,cos2,..,cosN
// Indeed, cos and sin give the same matrices

 int mode;
 int mstart = m1;
 int mstop = (m2<0) ? locBCs->numModes+1 : m2;

 long matSize;

 int storage = solInfo().getFetiInfo().local_cntl->subtype;

 switch (storage) {
   default:
   case 0:
     for (mode=mstart; mode<mstop; ++mode) {
         KC[mode] = new SkyMatrixC(nodeToNode, dsa, c_dsa, sinfo.solvercntl->trbm);
         matSize = (KC[mode]) ? KC[mode]->size() : 0;
#if defined(sgi) && ! defined(_OPENMP)
         __add_and_fetch(&totMemSky, matSize);
         __add_and_fetch(&memK, matSize);
#else
         totMemSky += matSize;
         memK += matSize;
#endif

         if (KiiC) {
            KiiC[mode] = new SkyMatrixC(nodeToNode, dsa, sinfo.solvercntl->trbm,
                                        glInternalMap,0);
            matSize = KiiC[mode]->size();
#if defined(sgi) && ! defined(_OPENMP)
            __add_and_fetch(&memPrec, matSize);
#else
            memPrec += matSize;
#endif
         }
      }
      break;
   case 1:
      for (mode=mstart; mode<mstop; ++mode) {
         KC[mode] = new BLKSparseMatrixC(nodeToNode, dsa, c_dsa,
                                         sinfo.solvercntl->trbm, *sinfo.solvercntl);
         matSize = (KC[mode]) ? KC[mode]->size() : 0;
#if defined(sgi) && ! defined(_OPENMP)
         __add_and_fetch(&totMemSparse, matSize);
         __add_and_fetch(&memK, matSize);
#else
         totMemSparse += matSize;
         memK += matSize;
#endif

         if (KiiC) {
            KiiC[mode] = new BLKSparseMatrixC(nodeToNode,dsa,glInternalMap,
                             sinfo.solvercntl->trbm, *sinfo.solvercntl);
            matSize = KiiC[mode]->size();
#if defined(sgi) && ! defined(_OPENMP)
            __add_and_fetch(&memPrec, matSize);
#else
            memPrec += matSize;
#endif
         }
      }
      break;
 }

 for (mode=mstart; mode<mstop; ++mode) {

    KbbC[mode] = new DBComplexSparseMatrix(nodeToNode, dsa, glBoundMap);
    matSize = KbbC[mode]->size();
#if defined(sgi) && ! defined(_OPENMP)
    __add_and_fetch(&memPrec, matSize);
#else
    memPrec += matSize;
#endif

 }

 if (locBCs->numDir>0) {

   int cLen = dsa->size();
   int* cBoundMap = (int*) dbg_alloca(sizeof(int)*cLen);
   int cBoundLen = 0;

   for (int iDof = 0; iDof < cLen; ++iDof) {
     int dofI = c_dsa->getRCN(iDof);
     cBoundMap[iDof] = (dofI < 0) ? cBoundLen++:-1;
   }

   for (mode=mstart; mode<mstop; ++mode) {
     KucC[mode] = new CuCComplexSparse(nodeToNode, dsa, c_dsa);
     KccC[mode] = new DBComplexSparseMatrix(nodeToNode, dsa, cBoundMap);
   }

 }

}


void
MDAxiData::Assemble(int m1, int m2) {

 int i;
 int iele;
 int maxele = numElements();
 int maxDof = maxNumDOF();
 //int MaxNodes = numNode();
 int mode;
 int mstart = m1;
 int mstop = (m2<0) ? locBCs->numModes+1 : m2;

 int storage = solInfo().getFetiInfo().local_cntl->subtype;

 double *firstpointer = (double *) dbg_alloca(sizeof(double)*maxDof*maxDof);
 double *secondpointer = (double *) dbg_alloca(sizeof(double)*maxDof*maxDof);
 double kappa=hParams->kappaHelm;

 ComplexSparseMatrix *spm;
 FullSquareMatrix keltmp(maxDof, secondpointer);

//
// All matrices constructors put the matrix at zero
// So no initialisation of the Ks is performed
//

 for (iele=0; iele<maxele; ++iele) {

    AxiHElement *elem = dynamic_cast<AxiHElement *>(packedEset[iele]);
    if (elem == 0)  {
       fprintf(stderr,"Element chosen non axisymmetric. Aborting\n");
       exit(1);
    }

    FullSquareMatrix kel = elem->stiffness(nodes,firstpointer);

    kel *= 2*M_PI;
    int done = 0;

    for (mode=mstart; mode<mstop; ++mode) {

       if ((done==0) && (mode>0)) {
         kel *= 0.5;
         done = 1;
       }

       switch (storage) {
         default:
         case 0:
           spm = (SkyMatrixC*) KC[mode];
           spm->add(kel,(*allDOFs)[iele]);
           if (KiiC) {
              spm = (SkyMatrixC*) KiiC[mode];
              spm->add(kel,(*allDOFs)[iele]);
           }
           break;
         case 1:
           spm = (BLKSparseMatrixC*) KC[mode];
           spm->add(kel,(*allDOFs)[iele]);
           if (KiiC) {
              spm = (BLKSparseMatrixC*) KiiC[mode];
              spm->add(kel,(*allDOFs)[iele]);
           }
           break;
       }

       KbbC[mode]->add(kel,(*allDOFs)[iele]);

       if (KucC) KucC[mode]->add(kel,(*allDOFs)[iele]);
       if (KibC) KibC[mode]->add(kel,(*allDOFs)[iele]);
       if (KccC) KccC[mode]->add(kel,(*allDOFs)[iele]);

    }

    FullSquareMatrix ktheta = elem->stiffteta(nodes,firstpointer);

    ktheta *= M_PI;

    for (mode=mstart; mode<mstop; ++mode) {

       if (mode==0)
         continue;

       keltmp.copy(ktheta.data());

       keltmp *= mode*mode;

       switch (storage) {
         default:
         case 0:
           spm = (SkyMatrixC*) KC[mode];
           spm->add(keltmp,(*allDOFs)[iele]);
           if (KiiC) {
              spm = (SkyMatrixC*) KiiC[mode];
              spm->add(keltmp,(*allDOFs)[iele]);
           }
           break;
         case 1:
           spm = (BLKSparseMatrixC*) KC[mode];
           spm->add(keltmp,(*allDOFs)[iele]);
           if (KiiC) {
              spm = (BLKSparseMatrixC*) KiiC[mode];
              spm->add(keltmp,(*allDOFs)[iele]);
           }
           break;
       }

       KbbC[mode]->add(keltmp,(*allDOFs)[iele]);

       if (KucC) KucC[mode]->add(keltmp,(*allDOFs)[iele]);
       if (KibC) KibC[mode]->add(keltmp,(*allDOFs)[iele]);
       if (KccC) KccC[mode]->add(keltmp,(*allDOFs)[iele]);

    }

 }

 if (locBCs->numSomm>0) {

   int somDofs = locBCs->somBC[0]->numDofs();
   int *somP = (int *) dbg_alloca(sizeof(int)*somDofs);
   DComplex *buffSom = (DComplex *) dbg_alloca(sizeof(DComplex)*somDofs*somDofs);

   for (i=0; i<locBCs->numSomm; ++i) {

      for (mode=mstart; mode<mstop; ++mode) {

         locBCs->somBC[i]->dofs(*dsa, somP);

         FullSquareMatrixC ks = locBCs->somBC[i]->turkelMatrix(nodes,
                                kappa, mode, buffSom);

         if (mode==0)
            ks *= DComplex(-2*M_PI,0.0);
         else
            ks *= DComplex(-M_PI,0.0);

         switch (storage) {
           default:
           case 0:
             spm = (SkyMatrixC*) KC[mode];
             spm->add(ks, somP);
             if (KiiC) {
                spm = (SkyMatrixC*) KiiC[mode];
                spm->add(ks, somP);
             }
             break;
           case 1:
             spm = (BLKSparseMatrixC*) KC[mode];
             spm->add(ks, somP);
             if (KiiC) {
                spm = (BLKSparseMatrixC*) KiiC[mode];
                spm->add(ks, somP);
             }
             break;
         }

         KbbC[mode]->add(ks,somP);

         if (KucC) KucC[mode]->add(ks,somP);
         if (KibC) KibC[mode]->add(ks,somP);
         if (KccC) KccC[mode]->add(ks,somP);

      }

   }

 }

}


void
MDAxiData::addInterface(int *SubSign, int m1, int m2) {

 int iInter;
 int i;
 int mode;
 int mstart = m1;
 int mstop = (m2<0) ? locBCs->numModes+1 : m2;
 int numNeighb = scomm->numNeighb;
 double kappa=hParams->kappaHelm;

 ComplexSparseMatrix *spm;
 FullSquareMatrix ksi;

 int storage = solInfo().getFetiInfo().local_cntl->subtype;
 int typeInterf = solInfo().getFetiInfo().lumpedinterface;

 int numDofs;

 // Assumption : all elements on the interface have the same
 //              number of dofs

 if (numNeighb==0)
    numDofs = 0;
 else
    numDofs = InterfElem[0]->numDofs();

 int *dofs = (int *) dbg_alloca(sizeof(int)*numDofs);
 double *pointer = (double *) dbg_alloca(sizeof(double)*numDofs*numDofs);

 for (iInter = 0; iInter < numNeighb; ++iInter) {

   if (SubSign[scomm->getNeighb<DComplex>(iInter)->subNum()] != InterfSign) {

      for (i=InterfaceIndex[iInter]; i<InterfaceIndex[iInter+1]; i++) {

        InterfElem[i]->dofs(*dsa, dofs);

        if (typeInterf==1)
           ksi = InterfElem[i]->interfMatrixLumped(nodes,pointer);
        else
           ksi = InterfElem[i]->interfMatrixConsistent(nodes,pointer);

        ksi *= kappa*InterfSign;

        for (mode=mstart; mode<mstop; ++mode) {

           switch (storage) {
             default:
             case 0:
               spm = (SkyMatrixC*) KC[mode];
               spm->addImaginary(ksi,dofs);
               if (KiiC) {
                  spm = (SkyMatrixC*) KiiC[mode];
                  spm->addImaginary(ksi,dofs);
               }
               break;
             case 1:
               spm = (BLKSparseMatrixC*) KC[mode];
               spm->addImaginary(ksi,dofs);
               if (KiiC) {
                  spm = (BLKSparseMatrixC*) KiiC[mode];
                  spm->addImaginary(ksi,dofs);
               }
               break;
           }

           KbbC[mode]->addImaginary(ksi,dofs);

           if (KibC) KibC[mode]->addImaginary(ksi,dofs);
           if (KucC) KucC[mode]->addImaginary(ksi,dofs);

        }

      }

    }

 }

}


void
MDAxiData::deleteGlMap(int **glBoundMap, int **glInternalMap) {

 delete[] glBoundMap[subNumber];
 delete[] glInternalMap[subNumber];

 glBoundMap[subNumber] = 0;
 glInternalMap[subNumber] = 0;

}


void
MDAxiData::factorKC(int m1, int m2) {

 int mode;
 int mstart = m1;
 int mstop = (m2<0) ? locBCs->numModes+1 : m2;

 for (mode=mstart; mode<mstop; ++mode)
     KC[mode]->factor();

}


void
MDAxiData::factorKiiC(int m1, int m2) {

 if (KiiC) {
   int mode;
   int mstart = m1;
   int mstop = (m2<0) ? locBCs->numModes+1 : m2;
   for (mode=mstart; mode<mstop; ++mode)
      KiiC[mode]->factor();
 }

}


void
MDAxiData::prepareCoarseData(DofSet ***allBoundary, int *counterInterf,
           int **FList, int **nonZE, int **subCount, DComplex ***subPosition) {

// Memory allocation for QSet
// Computation of connectivites for matrices Q

 DofSet **boundaryDOFs = allBoundary[subNum()];

 QSet = new ComplexSparseSet *[2*locBCs->numModes+1];

 localCVec = new ComplexVector *[2*locBCs->numModes+1];
 BKQ  = new DComplex *[2*locBCs->numModes+1];
 if (locMPCs) {
   CKQ = new DComplex *[2*locBCs->numModes+1];
 }

 int numWaveDir = solInfo().solvercntl->fetiInfo.numcgm;
 int numbEdges = 0;
 int *nonZeroEdges = new int[scomm->numNeighb];

 Connectivity &sharedNodes = *(scomm->sharedNodes);
 edgeDofSize = new int[scomm->numNeighb];

 // 1. Count number of edge dofs
 int iSub, iNode;
 int InterfCounter = 0;
 for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    edgeDofSize[iSub]=0;
    nonZeroEdges[iSub]=0;
    for (iNode = 0; iNode < sharedNodes.num(iSub); ++iNode)
      if (boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) {
         edgeDofSize[iSub] = numWaveDir;
         InterfCounter += 1;
         nonZeroEdges[iSub] += 1;
      }
    if (edgeDofSize[iSub]>0)
       numbEdges += 1;
 }

 int *Count = new int[numWaveDir*numbEdges];
 int **List = (int**) dbg_alloca(sizeof(int*)*scomm->numNeighb);
 DComplex **Position = new DComplex *[scomm->numNeighb];

 int i;
 int counter = 0;

 for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
     if (edgeDofSize[iSub]==0)
        continue;
     for (i=0; i<numWaveDir; ++i)
         Count[counter++] = nonZeroEdges[iSub];
 }

 for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
   if (edgeDofSize[iSub]==0)
      continue;
   List[iSub] = (int*) dbg_alloca(sizeof(int)*nonZeroEdges[iSub]);
   Position[iSub] = new DComplex[nonZeroEdges[iSub]];
   counter = 0;
   for (iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
     if (boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) {
       int dof = c_dsa->locate(sharedNodes[iSub][iNode],DofSet::Helm);
       if (dof>=0) {
         List[iSub][counter] = dof;
         Node nd = *nodes[sharedNodes[iSub][iNode]];
         Position[iSub][counter] = DComplex(nd.x,nd.y);
         counter++;
       }
     }
   }
 }

 int *FinalList = new int[numWaveDir*InterfCounter];

 counter = 0;
 for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
   if (edgeDofSize[iSub]==0)
      continue;
   for (i=0; i<numWaveDir; ++i) {
      for (iNode=0; iNode<nonZeroEdges[iSub]; ++iNode)
          FinalList[counter++] = List[iSub][iNode];
   }
 }

 setNumCoarseDofs();

 counterInterf[subNumber] = InterfCounter;
 FList[subNumber] = FinalList;
 nonZE[subNumber] = nonZeroEdges;
 subCount[subNumber] = Count;
 subPosition[subNumber] = Position;

}


void
MDAxiData::makeCoarseData(int *counterInterf, int **FList, int **nonZE,
           int **subCount, DComplex ***subPosition, int f1, int f2) {

// Computation of the matrices Q for enforcing Q^T.r = 0
// The size of the arrays is 2*locBCs->numModes+1 for all
// the basis functions 1,cos,sin,cos2,sin2,..,cosN,sinN

 int counter;
 int i, iSub, iNode;
 int numbEdges;
 int numWaveDir = solInfo().solvercntl->fetiInfo.numcgm;

 int InterfCounter = counterInterf[subNumber];
 int *FinalList = FList[subNumber];
 int *nonZeroEdges = nonZE[subNumber];
 DComplex **Position = subPosition[subNumber];
 int *Count = subCount[subNumber];

 int Fourier;
 int fstart = f1;
 int fstop = (f2==-1) ? 2*locBCs->numModes+1 : f2;

 CuCComplexSparse *Ktemp;
 DComplex *Coefs;

 for (Fourier=fstart; Fourier<fstop; ++Fourier) {

    Coefs = new DComplex[numWaveDir*InterfCounter];

    int *tmp = new int[numWaveDir*InterfCounter];
    for (i=0; i<numWaveDir*InterfCounter; ++i)
       tmp[i] = FinalList[i];

    counter = 0;
    numbEdges = 0;

    for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
       if (edgeDofSize[iSub]==0)
         continue;
       numbEdges += 1;
       double sign = (scomm->subNums[iSub] < subNumber) ? 1.0 : -1.0;
       for (i=0; i<numWaveDir; ++i) {
          // Setup the planar wave direction
          double kdx;
          double kdy;
          if (i == 0) {
            kdx = 0.0;
            kdy = 0.0;
            }
          else {
            kdx = hParams->kappaHelm*cos((2.0*M_PI*(i-1))/(numWaveDir-1));
            kdy = hParams->kappaHelm*sin((2.0*M_PI*(i-1))/(numWaveDir-1));
            }

          for (iNode=0; iNode<nonZeroEdges[iSub]; ++iNode)
             Coefs[counter++] = sign*CoarseGridValue(Fourier, kdx, kdy,
                                Position[iSub][iNode]);

       }
    }

    Ktemp = new CuCComplexSparse(numWaveDir*numbEdges, c_dsa->size(),
                                 Count, tmp, Coefs);

    QSet[Fourier] = new ComplexSparseSet(numInQSet);

    // NOTE : UH - 10/16/99
    // May need to generalize the position number in the SparseSet
    // QSet[Fourier]->setSparseMatrices(0, Ktemp);
    QSet[Fourier]->addSparseMatrix(Ktemp);

 }

}


void
MDAxiData::deleteCoarseInfo(int **FList, int **nonZE, int **subCount,
           DComplex ***subPosition) {

 int iSub;

 delete[] FList[subNumber];
 FList[subNumber] = 0;
 delete[] nonZE[subNumber];
 nonZE[subNumber] = 0;
 delete[] subCount[subNumber];
 subCount[subNumber] = 0;

 for (iSub=0; iSub<scomm->numNeighb; ++iSub) {
   delete[] subPosition[subNumber][iSub];
   subPosition[subNumber][iSub] = 0;
 }

 delete[] subPosition[subNumber];
 subPosition[subNumber] = 0;

}


DComplex
MDAxiData::CoarseGridValue(int Fourier, double kdx, double kdy,
           DComplex z) {

 int mode;

 mode = (Fourier%2==0) ? mode = Fourier/2 :
                         mode = (Fourier+1)/2;

 DComplex value;

 switch (solInfo().getFetiInfo().nonLocalQ) {
   default:
   case 0:
     value = exp(DComplex(0.0,kdx*real(z)));
     break;
   case 1:
     if (Fourier%2==0)
        value = DComplex(jn(-mode, kdx*real(z)),0.0);
     else
        value = DComplex(jn(mode, kdx*real(z)),0.0);
     break;
 }

 value *= exp(DComplex(0.0,kdy*imag(z)));

 return value;

}


void
MDAxiData::setNumCoarseDofs() {

 numCGridDofs = 0;
 int iSub;
 for (iSub = 0; iSub < scomm->numNeighb; ++iSub)
    numCGridDofs += edgeDofSize[iSub];

}


int
MDAxiData::numCoarseDofs() {

 return numCGridDofs;

}


void
MDAxiData::makeCKCt() {

 if (locMPCs) {
   CKCt = new SymFullMatrixC(locMPCs->totalNumMPC);
 }

}


void
MDAxiData::prepareMPCSet(int f1, int f2) {

// Computation of the matrices C for enforcing MPCs
// The size of the arrays is 2*locBCs->numModes+1 for all
// the basis functions 1,cos,sin,cos2,sin2,..,cosN,sinN

 int Fourier;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;
 int totalFourier = fstop-fstart;

 CuCComplexSparse *Ktemp;

 if (locMPCs == 0) {

   for (Fourier=fstart; Fourier<fstop; ++Fourier) {
      Ktemp = new CuCComplexSparse(0, c_dsa->size(), 0, 0, 0);
      // QSet[Fourier]->setSparseMatrices(1, Ktemp);
      QSet[Fourier]->addSparseMatrix(Ktemp);
   }

 }
 else {

   // UH - 1/4/00 -
   // We are in the particular case of one non-zero coefficient per
   // column for each Fourier mode

   int *Count = (int*) dbg_alloca(sizeof(int)*locMPCs->totalNumMPC);
   int *List = (int*) dbg_alloca(sizeof(int)*locMPCs->totalNumMPC);
   DComplex **Coefs = new DComplex*[totalFourier];

   int i,j;
   int counter = 0;

   for (i=0; i<locMPCs->totalNumMPC; ++i)
       Count[i] = 1;

   for (Fourier=fstart; Fourier<fstop; ++Fourier)
       Coefs[Fourier-fstart] = new DComplex[locMPCs->totalNumMPC];

   counter = 0;

   for (i=0; i<locMPCs->numMPCSet; ++i) {

      int dof = c_dsa->locate(locMPCs->term[i].nnum, (DofSet::Helm));

      if (dof<0) {
         fprintf(stderr,"\n\n ... MPC applied to constrained node ... \n");
         continue;
      }

      double theta;

      for (j=0; j<locMPCs->term[i].slices+1; ++j) {

        List[counter] = dof;

        if (locMPCs->term[i].angle1==locMPCs->term[i].angle2)
          theta=locMPCs->term[i].angle1;
        else
          theta=locMPCs->term[i].angle1 + (locMPCs->term[i].angle2-
                locMPCs->term[i].angle1)*j/locMPCs->term[i].slices;

        for (Fourier=fstart; Fourier<fstop; ++Fourier) {
           if (Fourier==0)
              Coefs[Fourier-fstart][counter] = DComplex(1.0/weight[dof], 0.0);
           else
                if (Fourier%2==0) {
                   int mode = Fourier/2;
                   Coefs[Fourier-fstart][counter] = DComplex(sin(mode*theta) /
                                             weight[dof], 0.0);
                }
                else {
                   int mode = (Fourier+1)/2;
                   Coefs[Fourier-fstart][counter] = DComplex(cos(mode*theta) /
                                             weight[dof], 0.0);
                }
        }

        counter++;

      }

   }

   for (Fourier=fstart; Fourier<fstop; ++Fourier) {

     int *tmp = new int[locMPCs->totalNumMPC];
     for (i=0; i<locMPCs->totalNumMPC; ++i)
         tmp[i] = List[i];

     Ktemp = new CuCComplexSparse(locMPCs->totalNumMPC, c_dsa->size(), Count,
                                  tmp, Coefs[Fourier-fstart]);

     // NOTE : UH - 10/16/99
     // May need to generalize the position number in the SparseSet

     // QSet[Fourier]->setSparseMatrices(1, Ktemp);
     QSet[Fourier]->addSparseMatrix(Ktemp);

   }

 }

}


void
MDAxiData::assembleCoarse(int Fourier, ComplexSparseMatrix *spm, FullMC **CoarseMPC) {

// Computation of Q^T K^-1 Q
// Storage of B K^-1 Q and C K^-1 Q and C K^-1 C^t

 int iRHS, iDof;
 int jk, jstart, jstop;
 //int iP = Fourier%threadManager->numThr();

 // Assumption : Sizes do not depend on the current mode
 int nRHS = QSet[Fourier]->numCol();
 int numbEqs = KC[0]->neqs();
 int InterfLength = interfLen();
 int *allBoundDofs = (*scomm->sharedDOFs)[0];
 int locMPCSize = getMPCSize();

 int stepRHS = 16;
 stepRHS = (stepRHS>nRHS) ? nRHS : stepRHS;

/*
 DComplex **KQ = new DComplex*[stepRHS];
 for (iRHS=0; iRHS < stepRHS; ++iRHS)
   KQ[iRHS] = new DComplex[numbEqs];
*/

 DComplex **KQ = new DComplex*[stepRHS];
 DComplex *firstPointer = new DComplex[stepRHS*numbEqs];
 for (iRHS=0; iRHS < stepRHS; ++iRHS)
   KQ[iRHS] = firstPointer + iRHS*numbEqs;
// Changer le reSolve(nRHS, KQ) en reSolve(nRHS, KQ[0])

 DComplex *QtKQ = new DComplex[stepRHS*nRHS];

 //DComplex *pointer = (DComplex *) dbg_alloca(sizeof(DComplex)*4);
 //int *dofs = (int *) dbg_alloca(sizeof(int)*2);

 switch (solInfo().getFetiInfo().nonLocalQ) {
   default:
   case 0:
     localCVec[2*Fourier] = new ComplexVector(nRHS, DComplex(0.0, 0.0));
     BKQ[2*Fourier] = new DComplex[nRHS*InterfLength];
     if (locMPCSize)
       CKQ[2*Fourier] = new DComplex[numCGridDofs*locMPCSize];
     if (Fourier>0) {
       localCVec[2*Fourier-1] = new ComplexVector(nRHS, DComplex(0.0, 0.0));
       BKQ[2*Fourier-1] = new DComplex[nRHS*InterfLength];
       if (locMPCSize)
         CKQ[2*Fourier-1] = new DComplex[numCGridDofs*locMPCSize];
     }

     for (jk=0; jk<nRHS; jk+=stepRHS) {
       jstart= jk;
       jstop = stepRHS+jstart;
       jstop = (jstop>nRHS) ? nRHS : jstop;
       // Initialize QtKQ at 0.0
       for (iRHS=0; iRHS<stepRHS; ++iRHS)
         for (iDof=0; iDof<nRHS; ++iDof)
           QtKQ[iRHS*nRHS+iDof] = DComplex(0.0, 0.0);
       // Compute Qt K^-1 Q
       computeQtKQ(2*Fourier, jstart, jstop, KQ, QtKQ);
       // Compute B x ( K^-1 Q )
       for (iDof=0; iDof<InterfLength; iDof++) {
         for (iRHS=jstart; iRHS < jstop; ++iRHS) {
           BKQ[2*Fourier][iRHS*InterfLength+iDof]
                         =KQ[iRHS-jstart][allBoundDofs[iDof]];
         }
       }
       // Assemble QtKQ
       assembleQtKQ(2*Fourier, jstart, jstop, QtKQ, spm);
       // Assemble parts on MPC
       if (locMPCSize)
         assembleCKQSet(2*Fourier, jstart, jstop, QtKQ, CoarseMPC);
       if (Fourier>0) {
         // Initialize QtKQ at 0.0
         for (iRHS=0; iRHS<stepRHS; ++iRHS)
           for (iDof=0; iDof<nRHS; ++iDof)
             QtKQ[iRHS*nRHS+iDof] = DComplex(0.0, 0.0);
         // Compute Qt K^-1 Q
         computeQtKQ(2*Fourier-1, jstart, jstop, KQ, QtKQ);
         // Compute B x ( K^-1 Q )
         for (iDof=0; iDof<InterfLength; iDof++) {
           for (iRHS=jstart; iRHS < jstop; ++iRHS) {
             BKQ[2*Fourier-1][iRHS*InterfLength+iDof]
                           =KQ[iRHS-jstart][allBoundDofs[iDof]];
           }
         }
         // Assemble parts on MPC
         if (locMPCSize)
           assembleCKQSet(2*Fourier-1, jstart, jstop, QtKQ, CoarseMPC);
       }
     }
     break;
   case 1:
     localCVec[Fourier] = new ComplexVector(nRHS, DComplex(0.0, 0.0));
     BKQ[Fourier] = new DComplex[nRHS*InterfLength];
     if (locMPCSize)
       CKQ[Fourier] = new DComplex[numCGridDofs*locMPCSize];

     for (jk=0; jk<nRHS; jk+=stepRHS) {
       jstart= jk;
       jstop = jstart+stepRHS;
       jstop = (jstop>nRHS) ? nRHS : jstop;
       // Initialize QtKQ at 0.0
       for (iRHS=0; iRHS<stepRHS; ++iRHS)
         for (iDof=0; iDof<nRHS; ++iDof)
           QtKQ[iRHS*nRHS+iDof] = DComplex(0.0, 0.0);
       // Compute Qt K^-1 Q
       computeQtKQ(Fourier, jstart, jstop, KQ, QtKQ);
       // Compute B x ( K^-1 Q )
       for (iDof=0; iDof<InterfLength; iDof++) {
         for (iRHS=jstart; iRHS < jstop; ++iRHS) {
           BKQ[Fourier][iRHS*InterfLength+iDof]=
                                         KQ[iRHS-jstart][allBoundDofs[iDof]];
         }
       }
       // Assemble QtKQ
       assembleQtKQ(Fourier, jstart, jstop, QtKQ, spm);
       // Assemble parts on MPC
       if (locMPCSize)
         assembleCKQSet(Fourier, jstart, jstop, QtKQ, CoarseMPC);
     }
     break;
 }

 delete[] QtKQ;
 QtKQ = 0;

/*
 for (iRHS=0; iRHS<stepRHS; ++iRHS)
   delete[] KQ[iRHS];
 delete[] KQ;
*/

 delete[] firstPointer;
 firstPointer = 0;
 delete[] KQ;
 KQ = 0;

}


void
MDAxiData::computeQtKQ(int Fourier, int jstart, int jstop, DComplex **KQ,
           DComplex *QtKQ) {

 int mode = (Fourier%2==0) ? Fourier/2 : (Fourier+1)/2;

 // QSet stores actually B^T Q and C^T
 QSet[Fourier]->multIdentity(KQ, jstart, jstop);

 // Compute K^-1 x Q
 KC[mode]->reSolve(jstop-jstart, KQ[0]);

 // Compute  - Q^T K^-1 Q
 QSet[Fourier]->transposeMultSubtract(jstop, KQ, QtKQ, jstart);

}


void
MDAxiData::assembleQtKQ(int Fourier, int jstart, int jstop, DComplex *QtKQ,
           ComplexSparseMatrix *spm) {

 int iRHS, iDof;
 int nRHS = QSet[Fourier]->numCol();

 if (jstart<numCGridDofs) {
   DComplex *pointer = (DComplex *) dbg_alloca(sizeof(DComplex)*4);
   pointer[0] = DComplex(0.0,0.0);
   pointer[2] = DComplex(0.0,0.0);
   pointer[3] = DComplex(0.0,0.0);
   int *dofs = (int *) dbg_alloca(sizeof(int)*2);
   int colStop = (jstop>numCGridDofs) ? numCGridDofs : jstop;
   for (iRHS=jstart; iRHS<colStop; ++iRHS) {

     dofs[1] = coarseEqNums[iRHS];

     for (iDof=0; iDof<iRHS; ++iDof) {
       FullSquareMatrixC k2(2, pointer);
       k2[0][1] = QtKQ[iDof+(iRHS-jstart)*nRHS];
       dofs[0] = coarseEqNums[iDof];
       spm->add(k2, dofs);
     }

     FullSquareMatrixC k1(1, pointer+1);
     k1[0][0] = QtKQ[iDof+(iRHS-jstart)*nRHS];
     spm->add(k1,dofs+1);

     for (iDof=iRHS+1; iDof<numCGridDofs; ++iDof) {
       FullSquareMatrixC k2(2, pointer);
       k2[0][1] = QtKQ[iDof+(iRHS-jstart)*nRHS];
       dofs[0] = coarseEqNums[iDof];
       spm->add(k2, dofs);
     }

   }
 }

}


void
MDAxiData::assembleCKQSet(int Fourier, int jstart, int jstop, DComplex *QtKQ,
           FullMC **CoarseMPC) {

 int iRHS, iDof;
 int nRHS = QSet[Fourier]->numCol();
 int locMPCSize = getMPCSize();

 for (iRHS=jstart; iRHS<jstop; ++iRHS) {
   if (iRHS<numCGridDofs) {
     // Store - C K^-1 Q & Update CoarseMPC[Fourier]
     int offSet = iRHS*locMPCSize;
     for (iDof=0; iDof<locMPCSize; ++iDof)
       CKQ[Fourier][offSet+iDof]=QtKQ[iDof+numCGridDofs+(iRHS-jstart)*nRHS];
     for (iDof=0; iDof<locMPCSize; ++iDof) {
       (*CoarseMPC[Fourier])[coarseEqNums[iRHS]][localToGlobalMPC[iDof]]
            += QtKQ[iDof+numCGridDofs+(iRHS-jstart)*nRHS];
     }
   }
   else {
     // Update CKCt
     int colRight=iRHS-numCGridDofs;
#if defined(sgi) && ! defined(_OPENMP)
     ussetlock(allocLock);
     for (iDof=0; iDof<=colRight; ++iDof) {
       (*CKCt)[iDof][colRight-iDof]+=QtKQ[iDof+numCGridDofs+(iRHS-jstart)*nRHS];
     }
     usunsetlock(allocLock);
#else
     for (iDof=0; iDof<=colRight; ++iDof) {
       (*CKCt)[iDof][colRight-iDof]+=QtKQ[iDof+numCGridDofs+(iRHS-jstart)*nRHS];     }
#endif
   }
 }

}



void
MDAxiData::deleteConn() {

 delete elemToNode;
 delete nodeToElem;
 elemToNode = 0;
 nodeToElem = 0;

}


int*
MDAxiData::getCoarseGridDofs(DofSetArray *coarseEqs, Connectivity &subToEdge) {

 if (coarseEqNums)
    return coarseEqNums;
 else
    coarseEqNums = new int[numCGridDofs];

 int offset = 0;
 int myNum = subNum();
 int iNeighb, i;

 for (iNeighb=0; iNeighb<scomm->numNeighb; ++iNeighb) {
    int fDof = coarseEqs->firstdof(subToEdge[myNum][iNeighb]);
    for (i=0; i<edgeDofSize[iNeighb]; ++i)
       coarseEqNums[offset++] = fDof+i;
 }

 return coarseEqNums;

}


void
MDAxiData::makeCoarseGridDofs(DofSetArray *coarseEqs, Connectivity *subToEdge) {

 if (coarseEqNums)
    return;
 else
    coarseEqNums = new int[numCGridDofs];

 int offset = 0;
 int myNum = subNum();
 int iNeighb, i;

 for (iNeighb=0; iNeighb<scomm->numNeighb; ++iNeighb) {
    int fDof = coarseEqs->firstdof((*subToEdge)[myNum][iNeighb]);
    for (i=0; i<edgeDofSize[iNeighb]; ++i)
       coarseEqNums[offset++] = fDof+i;
 }

}


void
MDAxiData::buildRHS(DComplex *force, int f1, int f2) {

// Build the local RHS forces from Neumann and Dirichlet
// boundary conditions.

 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;
 int localSize = localLen();
 int Fourier;
 int i;
 int offSet;

 // ... COMPUTE EXTERNAL FORCE FROM COMPLEX NEUMAN BC
 for (i=0; i < locBCs->numNeu; ++i) {

  double *x = (double *) dbg_alloca(3*sizeof(double));
  double *y = (double *) dbg_alloca(3*sizeof(double));

  int dof1 = c_dsa->locate(locBCs->neuBC[i].nnum1, (1 << 7));
  if (dof1 < 0)
     continue;
  Node nd = nodes.getNode(locBCs->neuBC[i].nnum1);
  x[0] = nd.x;
  y[0] = nd.y;

  int dof2 = c_dsa->locate(locBCs->neuBC[i].nnum2, (1 << 7));
  if (dof2 < 0)
     continue;
  nd = nodes.getNode(locBCs->neuBC[i].nnum2);
  x[1] = nd.x;
  y[1] = nd.y;

  int dof3 = (locBCs->neuBC[i].nnum3==-1) ? -1 :
              c_dsa->locate(locBCs->neuBC[i].nnum3, (1 << 7));
  if (dof3>=0) {
    nd = nodes.getNode(locBCs->neuBC[i].nnum3);
    x[2] = nd.x;
    y[2] = nd.y;
  }

  double nx = 0;
  double ny = 0;

  getNormal2D(locBCs->neuBC[i].nnum1,locBCs->neuBC[i].nnum2,nx,ny);

  offSet = fstart*localSize;

  for (Fourier=fstart; Fourier<fstop; ++Fourier) {
     DComplex *f = (DComplex*) dbg_alloca(sizeof(DComplex)*3);
     f[0] = DComplex(0.0, 0.0);
     f[1] = DComplex(0.0, 0.0);
     f[2] = DComplex(0.0, 0.0);
     if (dof3==-1) {
       locBCs->IntNeu2(Fourier,i,hParams->kappaHelm,x,y,nx,ny,f);
       if (dof1 >= 0) force[dof1+offSet] += f[0];
       if (dof2 >= 0) force[dof2+offSet] += f[1];
     }
     else {
       locBCs->IntNeu3(Fourier,i,hParams->kappaHelm,x,y,nx,ny,f);
       if (dof1 >= 0) force[dof1+offSet] += f[0];
       if (dof2 >= 0) force[dof2+offSet] += f[1];
       if (dof3 >= 0) force[dof3+offSet] += f[2];
     }
     offSet += localSize;
  }

 }

 if (locBCs->numDir>0) {

    offSet = fstart*localSize;
    int mode;

    ComplexVector Vc(locBCs->numDir, DComplex(0.0,0.0));

    for (Fourier=fstart; Fourier<fstop; ++Fourier) {

      mode = (Fourier%2==0) ? mode = Fourier/2 :
                              mode = (Fourier+1)/2;

      // CONSTRUCT THE NON-HOMOGENEOUS COMPLEX DIRICHLET BC VECTOR
      for (i=0; i< locBCs->numDir; ++i) {
         int dof = dsa->locate(locBCs->dirBC[i].nnum, (1 << 7));
         if (dof < 0) continue;
         int dof1 = c_dsa->invRCN(dof);
         if (dof1 >= 0) {
            Node nd = nodes.getNode(locBCs->dirBC[i].nnum);
            Vc[dof1] = locBCs->Dir2Four(Fourier,i,hParams->kappaHelm,
                                        nd.x,nd.y);
         }
      }
      // ... PERFORM MULTIPLICATION TO GET NON-HOMOGENEOUS FORCE:
      //                    Fnh -= [Kuc]*[Vc]
      if (KucC) KucC[mode]->multSubtract(Vc.data(),force+offSet);
      offSet += localSize;
      Vc.copy(DComplex(0.0,0.0));

    }

 }

}


void
MDAxiData::assembleDUDN(DComplex *u, DComplex **dudn, int f1, int f2) {

// Compute on the boundary of the mesh with Dirichlet BC
// the integral of du/dn and the shape functions

 if (locBCs->numDir==0)
   return;

 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;

 DComplex *buff1 = (DComplex *) dbg_alloca(sizeof(DComplex)*locBCs->numDir);
 DComplex *buff2 = (DComplex *) dbg_alloca(sizeof(DComplex)*locBCs->numDir);

 int i;
 int Fourier, mode;
 int localSize = localLen();
 int offSet = fstart*localSize;

 for (Fourier=fstart; Fourier<fstop; ++Fourier) {

    for (i=0; i<locBCs->numDir; ++i) {
      int dof = dsa->locate(locBCs->dirBC[i].nnum, (1 << 7));
      if (dof < 0) continue;
      int dof1 = c_dsa->invRCN(dof);
      if (dof1 >= 0) {
          Node nd = nodes.getNode(locBCs->dirBC[i].nnum);
          buff1[dof1] = locBCs->Dir2Four(Fourier,i,hParams->kappaHelm,
                                      nd.x,nd.y);
      }
    }

    mode = (Fourier%2==0) ? mode = Fourier/2 :
                            mode = (Fourier+1)/2;

    if (KccC) KccC[mode]->mult(buff1, buff2);
    if (KucC) KucC[mode]->transposeMultAdd(u+offSet, buff2);

    for (i=0; i<locBCs->numDir; ++i) {
      int dof2 = dsa->locate(locBCs->dirBC[i].nnum,(1 << 7));
      dof2 = c_dsa->invRCN(dof2);
      if(dof2 >= 0) {
        dudn[Fourier][locToglDir[i]] += buff2[dof2];
      }
    }

    offSet += localSize;

 }

}


void
MDAxiData::localDUDN(DComplex **u, DComplex **dudn, int *nodePos,
                     int *glToSurf, int f1, int f2) {

// Compute on an inner surface of the mesh defined by nodePos
// the integral of du/dn and the shape functions

 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;

 int i, j, iEle;
 int Fourier, mode;
 int maxele = numElements();
 int maxDof = maxNumDOF();
 //int MaxNodes = numNode();

 double *pointer = (double *) dbg_alloca(sizeof(double)*maxDof*maxDof);
 double *pointert = (double *) dbg_alloca(sizeof(double)*maxDof*maxDof);
 int *poinTmp = (int *) dbg_alloca(sizeof(int)*maxDof);

 for (iEle=0; iEle<maxele; ++iEle) {

   AxiHElement *elem = dynamic_cast<AxiHElement *>(packedEset[iEle]);
   int *dofs = elem->nodes(poinTmp);
   int buffer = 0;

   for (i=0; i<maxDof; ++i)
     if (nodePos[glNums[dofs[i]]]==-1) {
       buffer=-1;
       break;
     }
     else
       buffer += nodePos[glNums[dofs[i]]];

   if (buffer<1)
      continue;

   FullSquareMatrix kel = elem->stiffness(nodes,pointer);
   kel *= M_PI;
   FullSquareMatrix kelt = elem->stiffteta(nodes,pointert);
   kelt *= M_PI;

   for (i=0; i<maxDof; ++i) {

     if (nodePos[glNums[dofs[i]]]==0)
       continue;

     for (j=0; j<maxDof; ++j) {

       for (Fourier=fstart; Fourier<fstop; ++Fourier) {
         mode = (Fourier%2==0) ? mode = Fourier/2 :
                                 mode = (Fourier+1)/2;
         dudn[Fourier][glToSurf[glNums[dofs[i]]]] += kel[i][j]*
                      u[Fourier][glNums[dofs[j]]];
         if (Fourier==0) {
           dudn[Fourier][glToSurf[glNums[dofs[i]]]] += kel[i][j]*
                        u[Fourier][glNums[dofs[j]]];
         }
         else {
           dudn[Fourier][glToSurf[glNums[dofs[i]]]] +=
                        mode*mode*kelt[i][j]*u[Fourier][glNums[dofs[j]]];
         }
       }

     }

   }

 }

}


DComplex *
MDAxiData::getLocalCVec(int Fourier) {

 return localCVec[Fourier]->data();

}


void
MDAxiData::zeroLocalCVec(int f1, int f2) {

 int Fourier;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;

 for (Fourier=fstart; Fourier<fstop; ++Fourier)
   localCVec[Fourier]->zero();

}


void
MDAxiData::multLocalQtBK(DComplex *f, int f1, int f2) {

// Compute [ K^-1 B^T Q | K^-1 C^T ]^T ( -f )
// Store the result locally in localCVec

 int mode;
 int Fourier;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;
 int localSize = localLen();

 if (localSize == 0)
    return;

 DComplex *force = (DComplex *) dbg_alloca(sizeof(DComplex)*localSize);

 int iDof;

 for (Fourier=fstart; Fourier<fstop; ++Fourier) {

    localCVec[Fourier]->zero();

    for (iDof = 0; iDof<localSize; ++iDof)
       force[iDof] = -f[iDof+Fourier*localSize];

    mode = (Fourier%2==0) ? mode = Fourier/2 :
                            mode = (Fourier+1)/2;

    if (KC[mode])
       KC[mode]->reSolve(force);

    if (QSet[Fourier])
       QSet[Fourier]->transposeMultAdd(force,localCVec[Fourier]->data());

 }

}


void
MDAxiData::multKf(DComplex *u, DComplex *f, int shiftu, int shiftf, int f1,
           int f2) {

// For each mode, compute u = K^-1 f

 int Fourier;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;
 int localSize = localLen();
 int mode;

 for (Fourier=fstart; Fourier<fstop; ++Fourier) {
    mode = (Fourier%2==0) ? mode = Fourier/2 :
                            mode = (Fourier+1)/2;
    KC[mode]->solve(f+(Fourier-shiftf)*localSize,u+(Fourier-shiftu)*localSize);
 }

}


void
MDAxiData::computeBw(DComplex *Bw, DComplex *w, int f1, int f2) {

 int iDof;
 int totalInterfSize = scomm->sharedDOFs->numConnect();
 int *allBoundDofs = (*scomm->sharedDOFs)[0];
 int Fourier;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;
 int localSize = localLen();

 // Distribute the vector to the interface
 for(iDof = 0; iDof < totalInterfSize; ++iDof) {
   for (Fourier=fstart; Fourier<fstop; ++Fourier)
      Bw[iDof+Fourier*totalInterfSize] =
                           w[allBoundDofs[iDof]+(Fourier-fstart)*localSize];
 }

}


void
MDAxiData::sendInterf(DComplex *Vec) {

 int interfSize = interfLen();
 int iDof, iSub;
 int offset;

 for (iDof = 0; iDof < interfSize; ++iDof)
     interfBuff[iDof] = Vec[iDof];

 offset = 0;

 for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
     setSendData(iSub, interfBuff+offset);
     offset += scomm->sharedDOFs->num(iSub);
 }
}

void
MDAxiData::interfaceJump(DComplex *Vec)
{
 int iSub, iDof;
 int offset = 0;

 for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    DComplex *buff = (DComplex *) scomm->getExchangeData(iSub);
    for(iDof = 0; iDof < scomm->sharedDOFs->num(iSub); ++iDof)
       Vec[offset+iDof] -= buff[iDof];
    offset += scomm->sharedDOFs->num(iSub);
 }
}

void
MDAxiData::multKbbC(DComplex *rVec, DComplex *zVec, DComplex *deltaUVec,
                    DComplex *deltaFVec, int f1, int f2) {

 int *allBoundDofs = (*scomm->sharedDOFs)[0];
 int interfSize = interfLen();
 int localSize = localLen();

 DComplex *v   = (DComplex *) dbg_alloca(boundLen*sizeof(DComplex));
 DComplex *res = (DComplex *) dbg_alloca(boundLen*sizeof(DComplex));

 int iDof;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;
 int Fourier, mode;

 for (Fourier=fstart; Fourier<fstop; ++Fourier) {

    for(iDof = 0; iDof < boundLen; ++iDof) {
         v[iDof] = res[iDof] = DComplex(0.0, 0.0);
    }

    for (iDof = 0; iDof < interfSize; ++iDof) {
       v[dualToBoundary[iDof]] += rVec[Fourier*interfSize+iDof] *
                                  DComplex(scaling[iDof],0.0);
       if (deltaUVec)
          deltaUVec[Fourier*localSize+allBoundDofs[iDof]] =
                                              -v[dualToBoundary[iDof]];
    }

    mode = (Fourier%2==0) ? mode = Fourier/2 :
                            mode = (Fourier+1)/2;

    KbbC[mode]->mult(v,res);

    // NOTE - UH - 10/8/99
    // For the moment, the Dirichlet preconditioner is not available

    if (deltaFVec)
      for (iDof = 0; iDof < boundLen; ++iDof) {
         deltaFVec[Fourier*localSize+boundMap[iDof]] = res[iDof];
      }

    for (iDof = 0; iDof < interfSize; ++iDof)
       zVec[Fourier*interfSize+iDof] = res[dualToBoundary[iDof]] *
                                       DComplex(scaling[iDof],0.0);

 }

}


void
MDAxiData::sendDeltaF(DComplex *deltaF) {

 int offset = 0;
 int iSub, jDof;

 for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    int subLen =  scomm->sharedDOFs->num(iSub);
    for(jDof = 0; jDof < subLen; ++jDof)
       interfBuff[offset+jDof] = deltaF[(*scomm->sharedDOFs)[iSub][jDof]];
    setSendData(iSub, interfBuff+offset);
    offset += subLen;
 }
}

DComplex
MDAxiData::collectAndDotDeltaF(DComplex *deltaF)
{
 DComplex dot(0.0, 0.0);
 int iSub,iDof,jDof;

 for(iDof = 0; iDof < localLen(); ++iDof )
   dot += conj(deltaF[iDof])*deltaF[iDof];

 for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    DComplex *buff = (DComplex *) scomm->getExchangeData(iSub);
    int subLen = scomm->sharedDOFs->num(iSub);
    for(jDof = 0; jDof < subLen; ++jDof)
       dot += conj(deltaF[(*scomm->sharedDOFs)[iSub][jDof]])*buff[jDof];
 }
 return dot;
}


void
MDAxiData::multTransposeBKQ(DComplex *zl, int f1, int f2) {

// Compute [ (B K^-1 B^T Q)^T | (B K^-1 C^T)^T ] [ -zl ]

 int nRHS;
 int Fourier;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;
 int InterfLength = interfLen();
 int i,j;

 for (Fourier=fstart; Fourier<fstop; ++Fourier)
     localCVec[Fourier]->zero();

 if ((InterfLength == 0) || (localLen() == 0))
    return;

 for (Fourier=fstart; Fourier<fstop; ++Fourier) {
    nRHS = QSet[Fourier]->numCol();
    int offSet = 0;
    for (i=0; i<nRHS; ++i) {
       for (j=0; j<InterfLength; ++j) {
         (*localCVec[Fourier])[i] -= BKQ[Fourier][offSet+j]*
                                     zl[j+Fourier*InterfLength];
       }
       offSet += InterfLength;
    }
 }

}


void
MDAxiData::finishFPz(ComplexVector **w , DComplex *z, DComplex *FPz, int f1,
           int f2) {

// Assumption : same number of column in QSet per mode

 int Fourier;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;
 int totalFourier = fstop-fstart;
 int iDof;
 int InterfLength = interfLen();
 int totalInterfSize = scomm->sharedDOFs->numConnect();
 int localSize = localLen();
 int *allBoundDofs = (*scomm->sharedDOFs)[0];

 int numCDofs = QSet[0]->numCol();
 DComplex *wLocal = (DComplex *) dbg_alloca(sizeof(DComplex)*numCDofs);

 DComplex *f = new DComplex[totalFourier*localSize];

 for (iDof=0; iDof<totalFourier*localSize; ++iDof)
   f[iDof] = DComplex(0.0,0.0);

 for (Fourier=fstart; Fourier<fstop; ++Fourier) {

    for (iDof=0; iDof<numCDofs; ++iDof)
      wLocal[iDof] = (*w[Fourier])[coarseEqNums[iDof]];

    // Apply B^T z
    for (iDof=0; iDof<totalInterfSize; ++iDof)
      f[allBoundDofs[iDof]+(Fourier-fstart)*localSize] =
                                                  z[iDof+Fourier*InterfLength];

    // Perform f -= B^T Q w
    if (QSet[Fourier]) QSet[Fourier]->multSub(wLocal,
                                              f+(Fourier-fstart)*localSize);

 }

 DComplex *buff = new DComplex[totalFourier*localSize];
 multKf(buff, f, fstart, fstart, fstart, fstop);
 computeBw(FPz, buff, fstart, fstop);

 delete[] f;
 delete[] buff;

}


void
MDAxiData::getHalfInterf(DComplex *s, DComplex *t) {

 int iTg  = 0;
 int iOff = 0;
 int iSub, i, istart, istop;
 int numShared;

 for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    numShared = scomm->sharedDOFs->num(iSub);
    if(scomm->subNums[iSub] < subNumber) {
      istart = 0;
      istop  = numShared/2;
    } else {
      istart = numShared/2;
      istop  = numShared;
    }
    for(i = istart; i < istop; ++i)
      t[iTg++] = s[iOff+i];
    iOff += numShared;
 }

}


void
MDAxiData::getHalfInterf(DComplex *s,DComplex *t,DComplex *ss,
           DComplex *tt) {

 int iTg  = 0;
 int iOff = 0;
 int iSub, i, istart, istop;
 int numShared;

 for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    numShared = scomm->sharedDOFs->num(iSub);
    if (scomm->subNums[iSub] < subNumber) {
      istart = 0;
      istop = numShared/2;
    } else {
      istart = numShared/2;
      istop  = numShared;
    }
    for(i = istart; i < istop; ++i) {
      t[iTg]    = s[iOff+i];
      tt[iTg++] = ss[iOff+i];
    }
    iOff += numShared;
 }

}


void
MDAxiData::scatterHalfInterf(DComplex *s, DComplex *buffer) {

 int iTg = 0;
 int iOff = 0;
 int iSub, i, istop, imid;
 for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
    istop = scomm->sharedDOFs->num(iSub);
    imid = istop/2;
    if(scomm->subNums[iSub] < subNumber) {
      for(i = 0; i < imid; ++i)
        buffer[iOff+i] = s[iTg++];
    } else {
      for(i = imid; i < istop; ++i)
        buffer[iOff+i] = s[iTg++];
    }
   iOff += istop;
 }
}

void
MDAxiData::rebuildInterf(DComplex *v)
{
 int iOff = 0;
 int iSub, i, istop, imid;
 for(iSub = 0; iSub < scomm->numNeighb; ++iSub) {
   DComplex *buffer = (DComplex *) scomm->getExchangeData(iSub);
   istop = scomm->sharedDOFs->num(iSub);
   imid = istop/2;
   if(scomm->subNums[iSub] < subNumber)
     for(i = imid; i < istop; ++i)
       v[iOff+i] = -buffer[i];
   else
     for(i = 0; i < imid; ++i)
       v[iOff+i] = -buffer[i];
   iOff += istop;
 }
}

void
MDAxiData::subtractBt(DComplex *z, DComplex *f, int f1, int f2)
{
// Compute the vector f <- f - B^T z

 int Fourier;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;
 int iDof;
 int *allBoundDofs = (*scomm->sharedDOFs)[0];
 int localSize = localLen();
 int interfSize = interfLen();

 for (Fourier=fstart; Fourier<fstop; ++Fourier) {
   for (iDof=0; iDof<interfSize; ++iDof)
      f[Fourier*localSize+allBoundDofs[iDof]] -= z[iDof+Fourier*interfSize];
 }

}


void
MDAxiData::finishBtPz(ComplexVector **w, DComplex *z, DComplex *f) {

// Compute the vector f <- f - B^T z + B^T Q w

 if (localLen()==0)
    return;

 int Fourier;
 int totalFourier = 2*locBCs->numModes+1;
 int iDof;
 int *allBoundDofs = (*scomm->sharedDOFs)[0];
 int localSize = localLen();
 int interfSize = interfLen();

 for (Fourier=0; Fourier<totalFourier; ++Fourier) {
   for (iDof=0; iDof<interfSize; ++iDof)
      f[Fourier*localSize+allBoundDofs[iDof]] -= z[iDof+Fourier*interfSize];
   QSet[Fourier]->multAdd(w[Fourier]->data(),f+Fourier*localSize);
 }

}


int
MDAxiData::getMPCSize() {

 int size;

 if (locMPCs == 0)
    size = 0;
 else
    size = locMPCs->totalNumMPC;

 return size;

}


void
MDAxiData::mergeSolution(DComplex **globalSol, DComplex *u) {

 int dof, dof1;
 int inode;
 int Fourier;
 int totalFourier = 2*locBCs->numModes+1;
 int offSet = 0;
 int localSize = localLen();
 int i;

 DComplex** Vc = new DComplex *[totalFourier];

 if (locBCs->numDir>0) {
   for (Fourier=0; Fourier<totalFourier; ++Fourier) {
     Vc[Fourier] = new DComplex[locBCs->numDir];
     // CONSTRUCT THE NON-HOMOGENEOUS COMPLEX DIRICHLET BC VECTOR
     for (i=0; i< locBCs->numDir; ++i) {
         dof = dsa->locate(locBCs->dirBC[i].nnum, DofSet::Helm);
         if (dof < 0) continue;
         dof1 = c_dsa->invRCN(dof);
         if (dof1 >= 0) {
           Node nd = nodes.getNode(locBCs->dirBC[i].nnum);
           Vc[Fourier][dof1] = locBCs->Dir2Four(Fourier,i,
                               hParams->kappaHelm,nd.x,nd.y);
         }
     }
   }
 }

 for (inode = 0; inode < numnodes; ++inode) {

   if ((dof = c_dsa->locate(inode, (DofSet::Helm))) >= 0) {
     // Free dof
     offSet = 0;
     for (Fourier = 0; Fourier<totalFourier; ++Fourier) {
        globalSol[Fourier][glNums[inode]] = u[dof+offSet];
        offSet += localSize;
     }
   }
   else if ((dof = dsa->locate(inode, (DofSet::Helm))) >= 0) {
     // Constrained dof
     for (Fourier = 0; Fourier<totalFourier; ++Fourier) {
        dof1 = c_dsa->invRCN(dof);
        globalSol[Fourier][glNums[inode]] = Vc[Fourier][dof1];
     }
   }
   else {
     // Dof does not exist
     for (Fourier = 0; Fourier<totalFourier; ++Fourier)
        globalSol[Fourier][glNums[inode]] = DComplex(0.0,0.0);
   }

 }

 delete[] Vc;

}


void
MDAxiData::subtractQw(DComplex *f, ComplexVector **w, ComplexVector *mu,
           int f1, int f2) {

// Compute f -= Q [ w | mu ]
// Assumption : same number of column in QSet per mode

 int localSize = localLen();
 int locMPCSize = getMPCSize();
 DComplex *wLocal = (DComplex *) dbg_alloca(sizeof(DComplex)*
                    (locMPCSize+numCGridDofs));

 int i;
 int Fourier;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes + 1 : f2;

 for (i=0; i<locMPCSize; ++i)
    wLocal[i+numCGridDofs] = (*mu)[localToGlobalMPC[i]];

 for (Fourier=fstart; Fourier<fstop; ++Fourier) {

    for (i=0; i<numCGridDofs; ++i)
      wLocal[i] = (*w[Fourier])[coarseEqNums[i]];

    // Perform f = f - Q w
    if (QSet[Fourier]) QSet[Fourier]->multSub(wLocal, f+Fourier*localSize);

 }

}


void
MDAxiData::addCKQw(ComplexVector **lambdaVec, ComplexVector *mu,
                        int f1, int f2) {

// Compute - [ - ( C K^-1 B^T Q) | - (C K^-1 C^T) ] [lambdaVec | mu ]
// Put the resulting vector in localCVec[Fourier] on the "mu" part

 int locMPCSize = getMPCSize();
 if (locMPCSize==0)
    return;

 int i,k;
 int Fourier;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;

 for (Fourier=fstart; Fourier<fstop; Fourier++) {
   for (i=0; i<locMPCSize; ++i) {
     for (k=0; k<numCGridDofs; ++k) {
       (*localCVec[Fourier])[i+numCGridDofs] -= CKQ[Fourier][k*locMPCSize+i] *
                               (*lambdaVec[Fourier])[coarseEqNums[k]];
     }
   }
 }

 if (fstart==0) {
   for (i=0; i<locMPCSize; ++i) {
     (*localCVec[0])[i+numCGridDofs] -=
                       (*CKCt)[i][0]*(*mu)[localToGlobalMPC[i]];
     for (k=i+1; k<locMPCSize; ++k)
         (*localCVec[0])[i+numCGridDofs] -=
                       (*CKCt)[i][k-i]*(*mu)[localToGlobalMPC[k]];
     for (k=0; k<i; ++k)
         (*localCVec[0])[i+numCGridDofs] -=
                       (*CKCt)[k][i-k]*(*mu)[localToGlobalMPC[k]];
   }
 }

}


void
MDAxiData::finishFzl(DComplex *z, DComplex *mu, ComplexVector **lambda,
           DComplex *FPz, int f1, int f2) {

// Assumption : same number of column in QSet per mode

 int Fourier;
 int fstart = f1;
 int fstop = (f2<0) ? 2*locBCs->numModes+1 : f2;
 int totalFourier = fstop-fstart;
 int iDof;
 int InterfLength = interfLen();
 int totalInterfSize = scomm->sharedDOFs->numConnect();
 int localSize = localLen();
 int *allBoundDofs = (*scomm->sharedDOFs)[0];

 DComplex *f = new DComplex[totalFourier*localSize];

 for (iDof=0; iDof<totalFourier*localSize; ++iDof)
   f[iDof] = DComplex(0.0,0.0);

 for (Fourier=fstart; Fourier<fstop; ++Fourier) {
    // Apply B^T z
    for (iDof=0; iDof<totalInterfSize; ++iDof)
      f[allBoundDofs[iDof]+(Fourier-fstart)*localSize] =
                                               z[iDof+Fourier*InterfLength];
 }

 DComplex *buff = new DComplex[totalFourier*localSize];
 multKf(buff, f, fstart, fstart, fstart, fstop);
 computeBw(FPz, buff, fstart, fstop);

 delete[] f;
 delete[] buff;

 int i,j;
 int locMPCSize = getMPCSize();

 for (Fourier=fstart; Fourier<fstop; ++Fourier) {
    for (i=0; i<InterfLength; ++i) {
       for (j=0; j<numCGridDofs; ++j) {
          FPz[i+Fourier*InterfLength] -= BKQ[Fourier][j*InterfLength+i]*
                                     (*lambda[Fourier])[coarseEqNums[j]];
       }
       for (j=0; j<locMPCSize; ++j) {
          FPz[i+Fourier*InterfLength] -=
                                 BKQ[Fourier][(j+numCGridDofs)*InterfLength+i]
                                 * mu[localToGlobalMPC[j]];
       }
    }
 }

}



