#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>

#include <Driver.d/SubDomain.h>
#include <Feti.d/Feti.h>
#include <Feti.d/GMRESOrthoSet.h>
#include <Feti.d/GCROrthoSet.h>
#include <Feti.d/CGOrthoSet.h>
#include <Threads.d/PHelper.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/Skyline.d/DistSky.h>
#include <Math.d/Skyline.d/BlockSky.h>
#include <Math.d/Skyline.d/DistBlockSky.h>
#include <Math.d/matrix.h>
#include <Math.d/DistBLKSparse.h>
#include <Math.d/mathUtility.h>
#include <Utils.d/linkfc.h>
#include <Feti.d/CCtSolver.d/GlobalCCt.h>
#include <Feti.d/CCtSolver.d/BlockCCt.h>
#include <Feti.d/CCtSolver.d/SuperBlockCCt.h>
#include <Feti.d/CCtSolver.d/SubBlockCCt.h>
#include <Driver.d/Mpc.h>

#ifdef TFLOP 
#include <Math.d/mathUtility.h>
#endif

#include <Math.d/BLKSparseMatrix.h>
#include <Solvers.d/Rbm.h>
#include <Timers.d/GetTime.h>
#include <Math.d/VecSet.h>
#include <Math.d/FullMatrix.h>
#include <Utils.d/print_debug.h>

#ifndef MIN
#define MIN(X,Y) ((X < Y) ? X : Y)
#endif
#ifndef MAX
#define MAX(X,Y) (ScalarTypes::greaterThan(X,Y) ? X : Y)
#endif

#include <limits>

extern int verboseFlag;

extern "C" {
   void _FORTRAN(dsvdc)(double *, int &, int &, int&, double *,
                        double *, double *, int &, double *, int &,
                        double *, const int &, int &);

   void _FORTRAN(zsvdc)(complex<double> *, int &, int &, int&, complex<double> *,
                        complex<double> *, complex<double> *, int &, complex<double> *, int &,
                        complex<double> *, const int &, int &);
}

inline void Tsvdc(double *a, int &b, int &c, int&d, double *e,
                  double *f, double *g, int &h, double *i, int &j,
                  double *k, const int &l, int &m)
{
 _FORTRAN(dsvdc)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}

inline void Tsvdc(complex<double> *a, int &b, int &c, int&d, complex<double> *e,
                  complex<double> *f, complex<double> *g, int &h, complex<double> *i, int &j,
                  complex<double> *k, const int &l, int &m)
{
 _FORTRAN(zsvdc)(a,b,c,d,e,f,g,h,i,j,k,l,m);
}

extern double t0;
extern double t1;
extern double t2;
extern double t4;
extern double t5;
extern double t6;

// number of main loop iterations + number of friction orthogonalisation loop iterations
extern int iterTotal;

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
#endif

#include <Utils.d/Memory.h>

inline double DABS(double x) { return (x>0.0) ? x : -x; }

// New constructor for both shared and distributed memory
template<class Scalar>
GenFetiDPSolver<Scalar>::GenFetiDPSolver(int _nsub, GenSubDomain<Scalar> **_sd,
                  Connectivity *_subToSub, FetiInfo *_fetiInfo, FSCommunicator *_fetiCom,
                  int *_glSubToLoc, Connectivity *_mpcToSub, Connectivity *_mpcToSub_primal, Connectivity *_mpcToMpc,
                  Connectivity *_mpcToCpu, Connectivity *_cpuToSub, Connectivity *_bodyToSub,
                  GenSolver<Scalar> **sysMatrices, GenSparseMatrix<Scalar> **sysSparse,
                  Rbm **, int sandiaFlag, bool _computeRbms)
 : GenFetiSolver<Scalar>(_nsub, threadManager->numThr()), internalR(_nsub), internalC(_nsub), internalWI(_nsub) 
{
 filePrint(stderr," ... Create FETI-DP(H) solver");
 initialize();

 // Compute memory used by FETI Solver
 t6 -= getTime();
#ifdef DISTRIBUTED
 this->times.memoryFETI = -threadManager->memoryUsed();
#else
 this->times.memoryFETI -= memoryUsed();
#endif

 this->nsub       = _nsub;        // Number of subdomains
 this->sd         = _sd;          // pointer to Array of all Subdomains
 this->subToSub   = _subToSub;    // subdomain to subdomain connectivity
 this->mpcToSub   = _mpcToSub;    // MPC to subdomain connectivity
 this->glNumMpc = (this->mpcToSub) ? this->mpcToSub->csize() : 0;
 this->mpcToSub_primal = _mpcToSub_primal;
 this->glNumMpc_primal = (this->mpcToSub_primal) ? this->mpcToSub_primal->csize() : 0;
 mpcToMpc   = _mpcToMpc;    // MPC to MPC connectivity (PJSA: used for CC^t preconditioner
 mpcToCpu   = _mpcToCpu;
 this->cpuToSub   = _cpuToSub;
 this->fetiInfo   = _fetiInfo;    // Feti solver information
 this->fetiCom    = _fetiCom;       // PJSA
 this->glSubToLoc = _glSubToLoc;

  globalFlagCtc = domain->getNumCTC();
#ifdef DISTRIBUTED
  globalFlagCtc = this->fetiCom->globalMax((int) globalFlagCtc);
#endif

 //if(sandiaFlag) this->isDynamic = (sandiaFlag == 3) ? 1 : 0;
 //else this->isDynamic = (sysMatrices) ? 1: 0;

 //if((this->isDynamic && domain->probType() != SolverInfo::Modal) || this->fetiInfo->dph_flag) //CBM
 //  computeRbms = false;
 //else
 //  computeRbms = (this->glNumMpc || _computeRbms);
 if(geoSource->isShifted()) this->fetiInfo->project_g = false;

 computeRbms = _computeRbms;

 this->myCPU = this->fetiCom->cpuNum();
 this->numCPUs = this->fetiCom->size();
 bodyToSub = _bodyToSub;
 subToBody = bodyToSub->reverse();

 // Define FETI tolerance and maximum number of iterations
 double fetiTolerance = this->fetiInfo->tol;
 this->epsilon2 = fetiTolerance*fetiTolerance;
 this->maxiter  = this->fetiInfo->maxit;

 int iSub;
 // create this->vPat FSCommPattern object, used to send/receive a scalar vector (interfaceDOFs)
 this->vPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<Scalar>::CopyOnSend);
 for(iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setDofCommSize(this->vPat);
 this->vPat->finalize();

 // create this->sPat FSCommPattern objects, used to send/receive a single integer
 this->sPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
 for(iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setCommSize(this->sPat, 1);
 this->sPat->finalize();

 if(globalFlagCtc) {
   mpcPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
   for(iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setMpcCommSize(mpcPat);
   mpcPat->finalize();
 }

// DEBUG NONLINEAR
 if(this->fetiInfo->type == FetiInfo::nonlinear) {
   mpcSPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<Scalar>::CopyOnSend);
   for(iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setMpcCommSize(mpcSPat);
   mpcSPat->finalize();
 }

 // Classes to organize parallel execution of tasks
 this->fetiOps   = new GenFetiOp<Scalar> *[this->nsub];

 // Compute total this->interface length, total internal length
 // and total half this->interface length
 int tInterfLen    = 0;
 int tLocalLen     = 0;
 int tLocalRLen    = 0;
 this->halfSize = 0;
 for(iSub = 0; iSub < this->nsub; ++iSub) {
   this->fetiOps[iSub] = new GenFetiOp<Scalar>;
   this->interface.domLen[iSub] = this->sd[iSub]->interfLen();
   this->internalDI.domLen[iSub]  = this->sd[iSub]->numUncon();
   internalR.domLen[iSub] = this->sd[iSub]->localRLen();

   this->sd[iSub]->computeMasterFlag(this->mpcToSub);
   this->fetiOps[iSub]->setHalfOffset(this->halfSize);
   this->halfSize      += this->sd[iSub]->halfInterfLen();
   tInterfLen    += this->interface.domLen[iSub];
   tLocalLen     += this->internalDI.domLen[iSub];
   tLocalRLen    += internalR.domLen[iSub];
 }
 this->interface.len = tInterfLen;
 this->internalDI.len  = tLocalLen;
 internalR.len = tLocalRLen;

 // PJSA: compute the masterFlags
 bool *interfaceMasterFlag = new bool[tInterfLen];
 this->interface.computeOffsets();
 for(iSub = 0; iSub < this->nsub; ++iSub) {
   bool *subMasterFlag = this->sd[iSub]->getMasterFlag();
   int subOffset = this->interface.subOffset[iSub];
   int j;
   for(j=0; j<this->interface.domLen[iSub]; ++j)
     interfaceMasterFlag[subOffset+j] = subMasterFlag[j];
 }
 this->interface.setMasterFlag(interfaceMasterFlag);
 this->internalDI.setMasterFlag();
 internalR.setMasterFlag();
 // don't delete interfaceMasterFlag

 // Allocate space for reorthogonalization set
 this->times.memoryOSet -= memoryUsed();
 if((this->fetiInfo->outerloop == 0) || (this->fetiInfo->outerloop == 3))
   this->oSetCG = (this->fetiInfo->maxortho > 0) ? new GenCGOrthoSet<Scalar>(this->halfSize, this->fetiInfo->maxortho, this->fetiCom) : 0;
 else if(this->fetiInfo->outerloop == 1)
   this->oSetGMRES = new GenGMRESOrthoSet<Scalar>(this->halfSize, this->fetiInfo->maxortho, this->fetiCom);
 else
   this->oSetGCR = new GenGCROrthoSet<Scalar>(this->halfSize, this->fetiInfo->maxortho, this->fetiCom);
 this->times.memoryOSet += memoryUsed();
 filePrint(stderr,"       ... \n");

 if(sysMatrices == 0) {
   filePrint(stderr," ... Build Edge Augmentation (Q)");
   computeLocalWaveNumbers();
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::makeQ);  // build augmentation matrix
   if(this->fetiInfo->augment == FetiInfo::Gs) {
     // exchange number of each neighbors rbms
     paralApply(this->nsub, this->sd, &BaseSub::sendNumNeighbGrbm, this->sPat);
     this->sPat->exchange();
     paralApply(this->nsub, this->sd, &BaseSub::recvNumNeighbGrbm, this->sPat);
   }
   filePrint(stderr,"    ... \n");

   filePrint(stderr," ... Construct Subdomain Matrices");
   startTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);
   timedParal(this->times.consMatrix, this->nsub, this, &GenFetiSolver<Scalar>::constructMatrices);
   stopTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);
   filePrint(stderr,"   ... \n");

   filePrint(stderr," ... Assemble Subdomain Matrices");
   startTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);
   timedParal(this->times.consMatrix, this->nsub, this, &GenFetiSolver<Scalar>::assembleMatrices);
   stopTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);
   filePrint(stderr,"    ... \n");
 }
 else {
   for(iSub = 0; iSub < this->nsub; ++iSub) {
     this->sd[iSub]->Krr = sysMatrices[iSub];
     this->sd[iSub]->KrrSparse = sysSparse[iSub]; // XXXX
     this->fetiOps[iSub]->setSysMatrix(this->sd[iSub]->Krr, sysSparse[iSub]);
   }
 }

 // Compute stiffness scaling if required
 // HB: this also perform the LMPCs stiffness scaling/splitting for the primal method
 paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::initScaling);
 if((this->fetiInfo->scaling == FetiInfo::kscaling) || ((this->fetiInfo->mpc_scaling == FetiInfo::kscaling) && (this->glNumMpc_primal > 0)) 
    || (this->fetiInfo->augment == FetiInfo::WeightedEdges)) {
   execParal(this->nsub, this, &GenFetiSolver<Scalar>::sendScale);
   this->vPat->exchange();
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::collectScaling, this->vPat);
 }
 if(domain->solInfo().isCoupled) 
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::scaleAndSplitKww); // HB

 if(this->fetiInfo->augment == FetiInfo::WeightedEdges)
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::weightEdgeGs); // PJSA: W*Q

 // MPCs 
 mpcPrecon = false; 
 if(this->glNumMpc > 0) { 
   if(this->fetiInfo->mpc_scaling == FetiInfo::kscaling) { // MPC stiffness scaling
     FSCommPattern<Scalar> *mpcDiagPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU, 
                                                                   FSCommPattern<Scalar>::CopyOnSend);
     for(iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setMpcDiagCommSize(mpcDiagPat);
     mpcDiagPat->finalize();
     paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::sendMpcDiag, mpcDiagPat);
     mpcDiagPat->exchange();
     paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::collectMpcDiag, mpcDiagPat);
     delete mpcDiagPat;
   }

   if(this->fetiInfo->c_normalize) normalizeC();
 
   if(this->fetiInfo->mpc_precno == FetiInfo::diagCCt) {
     // use W scaling for preconditioning mpcs, don't need to build & invert CC^t
     paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::sendMpcScaling, this->vPat);
     this->vPat->exchange();  // neighboring subs mpc weights
     paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::collectMpcScaling, this->vPat);
   }
   else if(this->fetiInfo->mpc_precno != FetiInfo::noMpcPrec) {
     // used generalized proconditioner for mpcs, need to build and invert CC^t
     mpcPrecon = true;
     paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::initMpcScaling);
     buildCCt();
   }
 }

 // Factor matrices: K and Kii (if dirichlet preconditioner))
 filePrint(stderr," ... Factor Subdomain Matrices      ... \n");
 startTimerMemory(this->times.factor, this->times.memoryFactor);
 if(this->fetiInfo->solvertype != FetiInfo::spooles)
   timedParal(this->times.factorMat, this->nsub, this, &GenFetiDPSolver<Scalar>::factorLocalMatrices);
 else // PJSA: unresolved problem when factoring local spooles in parallel
   for(iSub=0; iSub<this->nsub; ++iSub) factorLocalMatrices(iSub);
 
 stopTimerMemory(this->times.factor, this->times.memoryFactor);
 if(this->fetiInfo->augment == FetiInfo::Gs) {
   this->makeRbmPat();
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::precondGrbm);
   // Get all the numbers of rigid body modes and dispatch RBMs to neighbors
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::sendInterfaceGrbm, this->rbmPat);
   this->rbmPat->exchange();
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::receiveInterfaceGrbm, this->rbmPat);
 }

 // Make coarse problems (Kcc^* and GtG)
 makeKcc();
 if(ngrbms /*&& this->glNumMpc*/) makeGtG();  // currently G = C^T*R (ie restriction of R to mpc interface)

 // build CC^t for preconditioning mpc residual if necessary
 // done earlier if(mpcPrecon) buildCCt();

 if(domain->solInfo().isCoupled) wetInterfaceComms();  // PJSA

 int tLocalCLen = 0;
 for(iSub = 0; iSub < this->nsub; ++iSub) {
   internalC.domLen[iSub] = this->sd[iSub]->numCoarseDofs();
   tLocalCLen += internalC.domLen[iSub];
 }
 internalC.len = tLocalCLen;

 if(domain->solInfo().isCoupled) {
   int tLocalWILen = 0;
   for(iSub = 0; iSub < this->nsub; ++iSub) {
     internalWI.domLen[iSub] = this->sd[iSub]->numWetInterfaceDofs();
     tLocalWILen += internalWI.domLen[iSub];
   }
   internalWI.len = tLocalWILen;
   internalWI.setMasterFlag();
 }

 // Allocate Distributed Vectors necessary for FETI solve loop
 this->times.memoryDV -= memoryUsed();
 int numC = (KccSolver) ? KccSolver->neqs() : 0;
 this->wksp = new GenFetiWorkSpace<Scalar>(this->interface, internalR, internalWI, ngrbms, numC, globalFlagCtc);
 this->times.memoryDV += memoryUsed();

 t6 += getTime();
}


template<class Scalar>
void
GenFetiDPSolver<Scalar>::computeLocalWaveNumbers()
{
  int iSub;
  if(this->fetiInfo->numdir > 0) { // PJSA 1-15-08, support for multiple fluids 
    if(this->fetiInfo->waveMethod == FetiInfo::uniform) {
      paralApply(this->nsub, this->sd, &BaseSub::computeWaveNumbers);
    }
    else if(this->fetiInfo->waveMethod == FetiInfo::averageK) {
      paralApply(this->nsub, this->sd, &BaseSub::computeWaveNumbers);
      // send/receive wave numbers for FETI-DPH EdgeWs augmentation
      FSCommPattern<double> *kPat = new FSCommPattern<double>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<double>::CopyOnSend);
      for(iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setCommSize(kPat, 4);
      kPat->finalize();
      paralApply(this->nsub, this->sd, &BaseSub::sendWaveNumbers, kPat);
      kPat->exchange();
      paralApply(this->nsub, this->sd, &BaseSub::collectWaveNumbers,kPat);
      delete kPat;
    }
    else {
      paralApply(this->nsub, this->sd, &BaseSub::averageMatProps);
      // send/receive neighb material properties for FETI-DPH EdgeWs augmentation
      FSCommPattern<double> *matPat = new FSCommPattern<double>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<double>::CopyOnSend);
      for(iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setCommSize(matPat, 5);
      matPat->finalize();
      paralApply(this->nsub, this->sd, &BaseSub::sendMatProps, matPat);
      matPat->exchange();
      paralApply(this->nsub, this->sd, &BaseSub::collectMatProps, matPat);
      delete matPat;
    }
  }
}


template<class Scalar> 
void
GenFetiDPSolver<Scalar>::makeKcc()
{
 int glNumSub = this->subToSub->csize();

 startTimerMemory(this->times.coarse1, this->times.memoryGtG);

#ifdef TFLOP
 map<int, int, less<int> > glCornerMap;
#else
 map<int, int> glCornerMap;
#endif

 int uniqueCorner=0;

 int *pointer = new int[glNumSub+1];

 int i, iSub;
#ifdef DISTRIBUTED
 for(i=0; i<glNumSub+1; ++i)
   pointer[i] = 0;

 for(iSub=0; iSub<this->nsub; ++iSub)
   pointer[this->sd[iSub]->subNum()] = this->sd[iSub]->numCorners();

 this->fetiCom->globalSum(glNumSub, pointer);

 int total = 0;
 for(iSub=0; iSub < glNumSub; ++iSub) {
   int tmp = pointer[iSub];
   pointer[iSub] = total;
   total += tmp;
 }
 pointer[glNumSub] = total;

 int *glCornerNodes = new int[total];
 for(i=0; i<total; ++i)
   glCornerNodes[i] = 0;

 for(iSub=0; iSub<this->nsub; ++iSub) {
   int numCorner = this->sd[iSub]->numCorners();
// RT - made this change so glCornerMap is constructed correctly when
// this function is called the second time
//   int *cornerNodes = this->sd[iSub]->getCornerNodes();
   int *localCornerNodes = this->sd[iSub]->getLocalCornerNodes();
   int *glN = this->sd[iSub]->getGlNodes();
   int iCorner;
   for(iCorner=0; iCorner<numCorner; ++iCorner)
// RT
//     glCornerNodes[pointer[this->sd[iSub]->subNum()]+iCorner] = cornerNodes[iCorner];
     glCornerNodes[pointer[this->sd[iSub]->subNum()]+iCorner] = glN[localCornerNodes[iCorner]];
 }

 this->fetiCom->globalSum(total, glCornerNodes);

 int iCorner;
 for(iCorner=0; iCorner<total; ++iCorner)
   if(glCornerMap.find(glCornerNodes[iCorner]) == glCornerMap.end() )
      glCornerMap[ glCornerNodes[iCorner] ] = uniqueCorner++;

 delete [] glCornerNodes;

#else
 int total = 0;
 for(iSub=0; iSub<this->nsub; ++iSub) {
   pointer[iSub] = total;
   int numCorner = this->sd[iSub]->numCorners();
   total += numCorner;

   // make sure this function returns global corner node numbers
// RT - made this change so glCornerMap is constructed correctly when
// this function is called the second time
//   int *cornerNodes = this->sd[iSub]->getCornerNodes();
   int *localCornerNodes = this->sd[iSub]->getLocalCornerNodes();
   int *glN = this->sd[iSub]->getGlNodes();

   int iCorner;
   for(iCorner=0; iCorner<numCorner; ++iCorner)
     if(glCornerMap.find(glN[localCornerNodes[iCorner]]) == glCornerMap.end()){
        glCornerMap[ glN[localCornerNodes[iCorner]] ] = uniqueCorner++;
     }
 }
 pointer[this->nsub] = total;
#endif

 glNumCorners = uniqueCorner;

 if(verboseFlag)
   filePrint(stderr," ... Total Number of Corners %5d  ...\n", glNumCorners);

 int *target = new int[total];
#ifdef DISTRIBUTED
 for(i=0; i<total; ++i) target[i] = 0;
#endif

 // set these numbers back 
 for(iSub=0; iSub<this->nsub; ++iSub) {
   int numCorner    = this->sd[iSub]->numCorners();
   int *cornerNodes = this->sd[iSub]->getCornerNodes();
   int *localCornerNodes = this->sd[iSub]->getLocalCornerNodes();
   int *glN = this->sd[iSub]->getGlNodes();
   int iCorner;
//    fprintf(stderr, "Sub %d found %d corners\n", iSub, numCorner);
   for(iCorner=0; iCorner<numCorner; ++iCorner) {
// RT - made this change so glCornerMap is constructed correctly when
// this function is called the second time
//      cornerNodes[iCorner] = glCornerMap[cornerNodes[iCorner] ];
      cornerNodes[iCorner] = glCornerMap[glN[localCornerNodes[iCorner] ] ];
      target[iCorner+pointer[this->sd[iSub]->subNum()]] = cornerNodes[iCorner];
   }
 }
#ifdef DISTRIBUTED
 this->fetiCom->globalSum(total, target);
#endif

 // PJSA: check for unsafe faces causing singularities in augmented part of coarse problem
 if(domain->solInfo().debug_icntl[4] == 1) {
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::initializeFaceSafety);
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::locateUnsafeFaces);
   paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::sendFaceSafetyInfo, this->sPat);
   this->sPat->exchange();
   paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::receiveFaceSafetyInfo, this->sPat);
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::locateUnsafeFaces);
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::printUnsafeFaces);
 }

 paralApply(this->nsub, this->sd, &BaseSub::findEdgeNeighbors);  // PJSA

 Connectivity *subToCorner = new Connectivity(glNumSub, pointer, target);

 if(cornerToSub) delete cornerToSub;  // PJSA
 cornerToSub = subToCorner->reverse();

 Connectivity *coarseConnectivity=0;
 Connectivity *coarseToSub=cornerToSub;
 Connectivity *subToCoarse=subToCorner;
 switch(this->fetiInfo->augment) {
   default:
   case FetiInfo::none:
     if(this->glNumMpc_primal > 0) {  // XXPA
       if(glNumCorners) { // standard case with corners and mpc equations
         coarseToSub = cornerToSub->merge(this->mpcToSub_primal);
         subToCoarse = coarseToSub->reverse();
         coarseConnectivity = coarseToSub->transcon(subToCoarse);
         //if(subToCoarse) { delete subToCoarse; subToCoarse=0; }
         //if(coarseToSub) { delete coarseToSub; coarseToSub=0; }
       } 
       else { // one processor case with mpc equations
         subToCoarse = this->mpcToSub_primal->reverse();
         coarseConnectivity = this->mpcToSub_primal->transcon(subToCoarse);
         //if(subToCoarse) { delete subToCoarse; subToCoarse=0; }
       }
     } 
     else { // standard case with corners only
       coarseConnectivity = cornerToSub->transcon(subToCorner);
     }
     break;
   case FetiInfo::Gs:
     if(this->glNumMpc_primal > 0) {  // XXPA
       Connectivity *temp = this->subToSub->merge(this->mpcToSub_primal);
       coarseToSub = cornerToSub->merge(temp);
       delete temp;
     } 
     else {
       coarseToSub = cornerToSub->merge(this->subToSub);
     }
     subToCoarse = coarseToSub->reverse();
     coarseConnectivity = coarseToSub->transcon(subToCoarse);
     //delete subToCoarse; subToCoarse=0;
     //delete coarseToSub; coarseToSub=0;
     break;
   case FetiInfo::WeightedEdges:
   case FetiInfo::Edges:
     if(!this->edgeToSub) makeEdgeConnectivity();
     if(this->glNumMpc_primal > 0) {  // XXPA
       Connectivity *temp = this->edgeToSub->merge(this->mpcToSub_primal);
       coarseToSub = cornerToSub->merge(temp);
       delete temp;
     } 
     else {
       coarseToSub = cornerToSub->merge(this->edgeToSub);
     }
     subToCoarse = coarseToSub->reverse();
     coarseConnectivity = coarseToSub->transcon(subToCoarse);
     //delete subToCoarse; subToCoarse=0;
     //delete coarseToSub; coarseToSub=0;
     break;
 }

 // Renumber the coarse problem equations
 int renumFlag = (this->fetiInfo->gtgSolver == FetiInfo::sparse) ? 0 : domain->solInfo().renum;
 compStruct renumber = coarseConnectivity->renumByComponent(renumFlag);
 if(cornerEqs) delete cornerEqs; // PJSA
 cornerEqs = new DofSetArray(coarseConnectivity->csize(), renumber.renum, 1);
 delete [] renumber.xcomp;
 // --------------------------------------------------------
 // HB: see if there are any independant components in Kcc
 //coarseConnectivity->print();
 if(verboseFlag) filePrint(stderr, " ... Kcc made of %3d blocks         ... \n",renumber.numComp);
 //fprintf(stderr, " ... cornerEqs.size() = %d\n",cornerEqs->size());
 //fprintf(stderr, " ... coarseConnectivity->csize()= %d\n",coarseConnectivity->csize());
 // --------------------------------------------------------

#if defined(USE_MUMPS) && defined(DISTRIBUTED)
 if(this->fetiInfo->gtgSolver == FetiInfo::mumps && domain->solInfo().mumps_icntl[18] == 3) { // matrix is distributed, use local graph for matrix structure
   delete coarseConnectivity;
   Connectivity *subToCPU = this->cpuToSub->reverse();
   coarseConnectivity = coarseToSub->transcon(subToCoarse, subToCPU->getTarget(), this->myCPU);
   delete subToCPU;
 }
#endif
 if(coarseToSub != cornerToSub) delete coarseToSub;
 if(subToCoarse != subToCorner) delete subToCoarse;

#ifdef DISTRIBUTED
 int mknd = DofSet::max_known_nonL_dof;
 int *glCornerDofs = new int[mknd*glNumCorners];
 for(i=0; i<mknd*glNumCorners; ++i)
   glCornerDofs[i] = 0;

 for(iSub=0; iSub<this->nsub; ++iSub) {
   int *glCornerNums = this->sd[iSub]->getCornerNodes();
   int numCorner     = this->sd[iSub]->numCorners();
   DofSet *cornerDofs= this->sd[iSub]->getCornerDofs();
   int iCorner;
   for(iCorner=0; iCorner<numCorner; ++iCorner) {
     if(cornerDofs[iCorner].contains(DofSet::Xdisp))
       glCornerDofs[mknd*glCornerNums[iCorner]+0] = 1;
     if(cornerDofs[iCorner].contains(DofSet::Ydisp))
       glCornerDofs[mknd*glCornerNums[iCorner]+1] = 1;
     if(cornerDofs[iCorner].contains(DofSet::Zdisp))
       glCornerDofs[mknd*glCornerNums[iCorner]+2] = 1;
     if(cornerDofs[iCorner].contains(DofSet::Xrot))
       glCornerDofs[mknd*glCornerNums[iCorner]+3] = 1;
     if(cornerDofs[iCorner].contains(DofSet::Yrot))
       glCornerDofs[mknd*glCornerNums[iCorner]+4] = 1;
     if(cornerDofs[iCorner].contains(DofSet::Zrot))
       glCornerDofs[mknd*glCornerNums[iCorner]+5] = 1;
     if(cornerDofs[iCorner].contains(DofSet::Temp))
       glCornerDofs[mknd*glCornerNums[iCorner]+6] = 1;
     if(cornerDofs[iCorner].contains(DofSet::Helm))
       glCornerDofs[mknd*glCornerNums[iCorner]+7] = 1; 
     if(cornerDofs[iCorner].contains(DofSet::IntPress))
       glCornerDofs[mknd*glCornerNums[iCorner]+8] = 1; 
   }
 }
 this->fetiCom->globalSum(mknd*glNumCorners, glCornerDofs);//mknd was 8 on 10-21-05 JF
 for(iSub=0; iSub<glNumSub; ++iSub) { 
   int numCorner = subToCorner->num(iSub);
   int iC;
   for(iC=0; iC<numCorner; ++iC) {
     int iCorner = (*subToCorner)[iSub][iC];
     if(glCornerDofs[mknd*iCorner+0] >= 1)
       cornerEqs->mark(iCorner, DofSet::Xdisp);
     if(glCornerDofs[mknd*iCorner+1] >= 1)
       cornerEqs->mark(iCorner, DofSet::Ydisp);
     if(glCornerDofs[mknd*iCorner+2] >= 1)
       cornerEqs->mark(iCorner, DofSet::Zdisp);
     if(glCornerDofs[mknd*iCorner+3] >= 1)
       cornerEqs->mark(iCorner, DofSet::Xrot);
     if(glCornerDofs[mknd*iCorner+4] >= 1)
       cornerEqs->mark(iCorner, DofSet::Yrot);
     if(glCornerDofs[mknd*iCorner+5] >= 1)
       cornerEqs->mark(iCorner, DofSet::Zrot);
     if(glCornerDofs[mknd*iCorner+6] >= 1)
       cornerEqs->mark(iCorner, DofSet::Temp);
     if(glCornerDofs[mknd*iCorner+7] >= 1)
       cornerEqs->mark(iCorner, DofSet::Helm);
     if(glCornerDofs[mknd*iCorner+8] >= 1)
       cornerEqs->mark(iCorner, DofSet::IntPress);
   }
 }
 delete [] glCornerDofs;

 // KHP: MPC MODIFICATION
 // number mpc equations after the corner equations
 if(this->glNumMpc_primal > 0) {  // XXPA
   // compute mpc equation offset
   mpcOffset = glNumCorners;
   if(this->fetiInfo->augment == FetiInfo::Gs) mpcOffset += glNumSub;
   if(this->fetiInfo->isEdgeAugmentationOn()) mpcOffset += this->edgeToSub->csize();
   int iMPC;
   for(iMPC=0; iMPC<this->mpcToSub_primal->csize(); ++iMPC)
     cornerEqs->setWeight(mpcOffset+iMPC, 1);
 }
 else mpcOffset = 0;

 int glNumRbm = 0;
 if(this->fetiInfo->augment == FetiInfo::Gs) {
   int *numRBMPerSub = new int[glNumSub];
   for(iSub=0; iSub<glNumSub; ++iSub)
     numRBMPerSub[iSub] = 0;  
   for(iSub=0; iSub<this->nsub; ++iSub) {
     numRBMPerSub[this->sd[iSub]->subNum()] = this->sd[iSub]->numRBM();
   }
   this->fetiCom->globalSum(glNumSub,numRBMPerSub);

   for(iSub=0; iSub<glNumSub; ++iSub) {
     cornerEqs->setWeight(glNumCorners + iSub, 
                          numRBMPerSub[iSub]);
     glNumRbm += numRBMPerSub[iSub];
   }
   delete [] numRBMPerSub;
 }

 if(this->fetiInfo->isEdgeAugmentationOn()) {
   int *edgeWeights = new int[this->edgeToSub->csize()];
   for(i=0; i<this->edgeToSub->csize(); ++i)
     edgeWeights[i] = 0;
   
   this->times.numEdges = 0;
   int jSub;
   for(iSub=0; iSub<this->nsub; ++iSub) {
     int numNeighbor = this->sd[iSub]->numNeighbors();
     int myNum = this->sd[iSub]->subNum();
     int jEdgeN = 0;
     for(jSub=0; jSub<numNeighbor; ++jSub) {
       int subJ = this->sd[iSub]->getSComm()->subNums[jSub];
       if(this->sd[iSub]->isEdgeNeighbor(jSub)) {
         if(myNum < subJ) { // PJSA
           int nEdge = this->sd[iSub]->numEdgeDofs(jSub);
           edgeWeights[(*this->subToEdge)[this->sd[iSub]->subNum()][jEdgeN]] = nEdge;
           //cornerEqs->setWeight(glNumCorners+(*this->subToEdge)[iSub][jSub], nEdge);
         }
         jEdgeN++;
       }
     }
   }
   this->fetiCom->globalSum(this->edgeToSub->csize(), edgeWeights);
   
   int iEdge;
   for(iEdge=0; iEdge<this->edgeToSub->csize(); ++iEdge) { 
     cornerEqs->setWeight(glNumCorners+iEdge, edgeWeights[iEdge]);
     this->times.numEdges += edgeWeights[iEdge];
   }
   // make sure I can delete this
   delete [] edgeWeights;
 }

#else
 for(iSub=0; iSub<this->nsub; ++iSub) {
   int *glCornerNums = this->sd[iSub]->getCornerNodes();
   int numCorner     = this->sd[iSub]->numCorners();
   DofSet *cdofs     = this->sd[iSub]->getCornerDofs();
   int iCorner;
   for(iCorner=0; iCorner<numCorner; ++iCorner) {
     cornerEqs->mark(glCornerNums[iCorner], cdofs[iCorner].list());
   }
 }

 if(this->glNumMpc_primal > 0) {  // XXPA
   // compute mpc equation offset
   mpcOffset = glNumCorners;
   if(this->fetiInfo->augment == FetiInfo::Gs) mpcOffset += glNumSub;
   if(this->fetiInfo->isEdgeAugmentationOn()) mpcOffset += this->edgeToSub->csize();
   int iMPC;
   for(iMPC=0; iMPC<this->mpcToSub_primal->csize(); ++iMPC)
     cornerEqs->setWeight(mpcOffset+iMPC, 1);
 }
 else mpcOffset = 0;

 int glNumRbm = 0;
 if(this->fetiInfo->augment == FetiInfo::Gs) {
   for(i=0; i<this->nsub; ++i) {
     int numRBM = this->sd[i]->numRBM();
     cornerEqs->setWeight(glNumCorners+this->sd[i]->subNum(), numRBM);
     glNumRbm += numRBM;
   }
   this->times.numRBMs = glNumRbm;
 } 
 // KHP: EDGE #dof = number of grbm for that edge's subdomain
 if(this->fetiInfo->isEdgeAugmentationOn()) {
   int j;
   this->times.numEdges = 0;
   for(i=0; i<this->nsub; ++i) {
     int numNeighbor = this->sd[i]->numNeighbors();
     int myNum = this->sd[i]->subNum();
     int jEdgeN = 0;
     for(j=0; j<numNeighbor; ++j) {
       int subJ = this->sd[i]->getSComm()->subNums[j];
       if(this->sd[i]->isEdgeNeighbor(j)) {
         if(myNum < subJ) {
           int nEdge = this->sd[i]->numEdgeDofs(j);
           cornerEqs->setWeight(glNumCorners+(*this->subToEdge)[i][jEdgeN], nEdge);
           this->times.numEdges += nEdge; 
         }
         jEdgeN++;
       }
     }
   }
 }
#endif
 delete subToCorner; subToCorner = 0;

 cornerEqs->finish();
 this->times.numCRNs = cornerEqs->size();

 if(verboseFlag) {
   filePrint( stderr, " ... Size of Interface  %7d     ...\n", 
              this->interface.len);
   filePrint( stderr, " ... Size of Global Kcc %7d     ...\n", 
              cornerEqs->size());
   filePrint( stderr, " ... global_rbm_tol     = %9.3e ...\n", 
              this->fetiInfo->grbm_tol);
   filePrint( stderr, " ... global_cor_rbm_tol = %9.3e ...\n", 
              this->fetiInfo->crbm_tol);
 }

 // PJSA: start new code **************************************************

 int nBodies = (bodyToSub) ? bodyToSub->csize() : 0;
 bool mbgflag = false;
 Connectivity *groupToBody = 0;
 if(this->glNumMpc_primal > 0) {

   // Find out which bodies are connected together by mpcs 
   // a collection of inter-connected bodies is referred to as a group
   Connectivity *subToMpc = this->mpcToSub_primal->reverse();
   Connectivity *bodyToMpc = bodyToSub->transcon(subToMpc);
   Connectivity *mpcToBody = bodyToMpc->reverse();
   Connectivity *bodyToBodyTmp = bodyToMpc->transcon(mpcToBody);
   Connectivity *bodyToBody = bodyToBodyTmp->modify();
   compStruct renumber = bodyToBody->renumByComponent(1);  // 1 = sloan renumbering
   nGroups = renumber.numComp;
   if(nGroups < nBodies) { // at least one multi-body group exists
     mbgflag = true;
     // make groupToBody connectivity
     renumber.order = new int[nBodies];
     for(i = 0; i < nBodies; ++i)
       renumber.order[renumber.renum[i]] = i;
     int *pointer = new int[nGroups + 1];
     pointer[0] = 0;
     for(i = 0; i < nGroups; ++i) {
       int nbod = renumber.xcomp[i + 1] - renumber.xcomp[i];
       pointer[i + 1] = pointer[i] + nbod;
     }
     groupToBody = new Connectivity(nGroups, pointer, renumber.order);
     groupToSub = groupToBody->transcon(bodyToSub);
     subToGroup = groupToSub->reverse();
   }
   if(renumber.xcomp) delete [] renumber.xcomp;
   if(renumber.renum) delete [] renumber.renum;
   delete subToMpc;
   delete bodyToMpc;
   delete mpcToBody;
   delete bodyToBodyTmp;
   delete bodyToBody;
 }
 if(!mbgflag) { // one body per group
   groupToSub = bodyToSub;
   subToGroup = subToBody;
   nGroups = nBodies;
 }

 // tell each subDomain what group it is in and find groups on each processor
 paralApply(this->nsub, this->sd, &BaseSub::setGroup, subToGroup);
 if(groups) delete [] groups;
 groups = new int[nGroups];  // groups represented on this processor
#ifdef DISTRIBUTED
 if(this->sd) {
   groups[0] = (*subToGroup)[this->sd[0]->subNum()][0];
   int n=1;
   for(i=1; i<this->nsub; ++i) {
     int group = (*subToGroup)[this->sd[i]->subNum()][0];
     int j;
     for(j=0; j<n; ++j) if(group == groups[j]) j=n+1;
     if(j==n) groups[n++] = group;
   }
   nGroups1 = n;  // number of groups represented on this processor
 }
 else nGroups1 = 0;
#else
 for(i=0; i<nGroups; ++i) groups[i] = i; 
 nGroups1 = nGroups;  
#endif
 if(verboseFlag) filePrint(stderr, " ... Number of bodies = %3d         ...\n", nGroups);
 
 int *groupProc = 0;
 if(computeRbms) {
  if((this->fetiInfo->corners == FetiInfo::noCorners) && (this->glNumMpc_primal > 0)) {
     // subdomain ZEMs need to be projected or eliminated
     // I think adding the dofs of primal mpcs to the "c" dofs would resolve this issue
     filePrint(stderr, " *** ERROR: mpc_type 2 is not supported with no corners *** \n"); 
     exit(-1);
  }
  if(verboseFlag) filePrint(stderr, " ... Computing multi-body GRBMs     ...\n");
  // calculate the centroid of each group
  double *centroid = new double[nGroups*3];   // pseudo centroid of each group
  double *nNodes = new double[nGroups];       // number of nodes in each group;
  for(i=0; i<nGroups; ++i) {
    nNodes[i] = 0.0;
    for(int j=0; j<3; ++j) centroid[3*i+j] = 0.0;
  }
  for(i=0; i<this->nsub; ++i) this->sd[i]->addNodeXYZ(centroid, nNodes);  // groups could be done in parallel
#ifdef DISTRIBUTED
  this->fetiCom->globalSum(nGroups, nNodes);
  this->fetiCom->globalSum(nGroups*3, centroid);
#endif
  for(i=0; i<nGroups; ++i) {
    if(nNodes[i] > 0) 
      for(int j=0; j<3; ++j) 
        centroid[3*i+j] = centroid[3*i+j]/nNodes[i];
  }
  // note: this centroid calculation assumes nodes are equally spaced, and also
  // interface nodes from a group split over more than one interface are used more than once.
  // however, the objective is only to find a point inside the group to be used as 
  // a reference for the geometric rbm calculation.  it is not necessary to use the exact
  // geometric centroid.
  delete [] nNodes;

  // make Zstar and R matrices for each subdomain
  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::makeZstarAndR, centroid);
  delete [] centroid;

  Connectivity *groupToMpc = 0;
  if(this->glNumMpc_primal > 0) {
    Connectivity *subToMpc = this->mpcToSub_primal->reverse();
    groupToMpc = groupToSub->transcon(subToMpc);
    paralApply(this->nsub, this->sd, &BaseSub::makeLocalToGroupMPC, groupToMpc);
    delete subToMpc;
  }

  // assemble global Zstar matrix for each body
  FullM **globalZstar = new FullM * [nGroups];
  int *zRow = new int[nGroups];
  int *zRowDim = new int[nGroups];
  int *zColDim = new int[nGroups];
  int *zColOffset = new int[nBodies];
  int zColDim1 = (this->sd) ? this->sd[0]->zColDim() : 0;  // (6 for 3D, 3 for 2D)
#ifdef DISTRIBUTED
  zColDim1 = this->fetiCom->globalMax(zColDim1);  // enforce it to be the same
#endif
  for(i=0; i<nGroups; ++i) {
    zRowDim[i] = 0; 
    if(mbgflag) {
      int cOff = 0;
      for(int j=0; j < groupToBody->num(i); ++j) {
        int body = (*groupToBody)[i][j];
        zColOffset[body] = cOff;
        cOff += zColDim1;
      }
      zColDim[i] = cOff;
    }
    else {
      zColOffset[i] = 0;
      zColDim[i] = zColDim1;
    }
  }
  for(iSub=0; iSub<this->nsub; ++iSub) {
    int subGroup = (*subToGroup)[this->sd[iSub]->subNum()][0];
    zRowDim[subGroup] += this->sd[iSub]->zRowDim();
  }
  if(this->glNumMpc_primal > 0) {
    for(i=0; i<nGroups; ++i) zRowDim[i] += groupToMpc->num(i);
  }
  if(groupToBody) delete groupToBody;

#ifdef DISTRIBUTED
  int *zRowOffset = new int[this->numCPUs*nGroups];
  for(i=0; i<this->numCPUs*nGroups; ++i) zRowOffset[i] = 0;
  for(i=0; i<nGroups1; ++i) {
    int iGroup = groups[i];
    for(int j=this->myCPU+1; j<this->numCPUs; ++j) zRowOffset[iGroup*this->numCPUs +j] = zRowDim[iGroup];
  }
  this->fetiCom->globalSum(nGroups, zRowDim);
  this->fetiCom->globalSum(this->numCPUs*nGroups, zRowOffset);
  for(i=0; i<nGroups; ++i) zRow[i] = zRowOffset[i*this->numCPUs + this->myCPU];
  delete [] zRowOffset;
#else
  for(i=0; i<nGroups; ++i) zRow[i] = 0;
#endif
  for(i=0; i<nGroups; ++i) {
    globalZstar[i] = new FullM(zRowDim[i], zColDim[i]);
    globalZstar[i]->zero();
  }
  // could do this in parallel (by groups)
  for(iSub=0; iSub<this->nsub; ++iSub) {
    int subBody = (*subToBody)[this->sd[iSub]->subNum()][0];
    int subGroup = (*subToGroup)[this->sd[iSub]->subNum()][0];
    if(this->sd[iSub]->zRowDim() > 0) 
      this->sd[iSub]->addSPCsToGlobalZstar(globalZstar[subGroup], zRow[subGroup], zColOffset[subBody]);
    if(this->sd[iSub]->numMPCs_primal() > 0) {
      int startRow = zRowDim[subGroup] - groupToMpc->num(subGroup);
      this->sd[iSub]->addMPCsToGlobalZstar(globalZstar[subGroup], startRow, zColOffset[subBody], zColDim1);
    }
  }
  if(this->glNumMpc_primal > 0) execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::setBodyRBMoffset, zColOffset);
  delete [] zColOffset;
  if(groupToMpc) delete groupToMpc;

  groupProc = new int[nGroups];
#ifdef DISTRIBUTED
  for(i=0; i<nGroups; ++i) {
    this->fetiCom->globalSum(zRowDim[i]*zColDim[i], globalZstar[i]->data());
    groupProc[i] = -1;
  }
  for(i=0; i<nGroups1; ++i) groupProc[groups[i]] = this->myCPU;
  for(i=0; i<nGroups; ++i) groupProc[i] = this->fetiCom->globalMax(groupProc[i]);
#else
  for(i=0; i<nGroups; ++i) groupProc[i] = this->myCPU;
#endif

  // now do svd on globalZstar for each group to get globalQ for each group
  ngrbmGr = new int[nGroups];
  for(i=0; i<nGroups; ++i) ngrbmGr[i] = 0;
  ngrbm = 0;  // total of all groups
  FullM  **Qtranspose;
  Qtranspose = new FullM * [nGroups];
  for(i=0; i<nGroups1; ++i) {
    int iGroup = groups[i];
    int ncol = zColDim[iGroup];  
    int nrow = zRowDim[iGroup];
    FullM U(ncol,ncol); U.zero(); 
    int rank = 0;
    singularValueDecomposition(*globalZstar[iGroup], U, ncol, nrow, rank, domain->solInfo().tolsvd);
    int ngrbmGrTmp = ncol - rank;
    globalZstar[iGroup]->clean_up();
    if(groupProc[iGroup] == this->myCPU) {
      ngrbmGr[iGroup] = ngrbmGrTmp;
      ngrbm += ngrbmGr[iGroup];
      fprintf(stderr, " ... Number of GRBMs for body %d: %d ...\n", iGroup, ngrbmGrTmp);
    }
    Qtranspose[iGroup] = new FullM(U, ngrbmGrTmp, rank, ncol, 0);
  }
#ifdef DISTRIBUTED
  ngrbms = this->fetiCom->globalSum(ngrbm);  // total number of rigid body modes for all processes
#else
  ngrbms = ngrbm;
#endif
  if(verboseFlag) filePrint(stderr, " ... total number of GRBMs = %5d  ...\n", ngrbms);

  delete [] zRow;
  delete [] zRowDim;
  delete [] zColDim;
  for(i=0; i<nGroups; ++i) delete globalZstar[i];
  delete [] globalZstar;

  // make local Rstar (and Rcstar if necessary)
  paralApply(this->nsub, this->sd, &BaseSub::makeLocalRstar, Qtranspose, this->fetiInfo->constrain_kcc);
  for(i=0; i<nGroups1; ++i) delete Qtranspose[groups[i]];
  delete [] Qtranspose;
 }
 else ngrbm = ngrbms = 0;

 if(c_cornerEqs) delete c_cornerEqs; // PJSA
 if(ngrbms && this->fetiInfo->constrain_kcc) {
  // assemble global Rc (restriction of block RBMs to corners) using new parallel algorithm &
  // create condensed cornerEqs DofSetArray so Kcc will always be non-singular
  int numCRNdof = cornerEqs->size();
  FullM globalRcstar(numCRNdof, ngrbm);
  globalRcstar.zero();
  paralApply(this->nsub, this->sd, &BaseSub::assembleGlobalRcstar, cornerEqs, &globalRcstar, ngrbmGr);

  // locate singularities (not necessary for dynamic analysis)
  int *sing = new int[ngrbm];
  for(i=0; i<ngrbm; ++i) sing[i] = -1;
  int svdcount = 0;
  int totRank = 0;
  double trbm2 = 0.0; // previously was domain->solInfo().trbm
  if(ngrbm > 0) {
    // assemble glCrnDofGroup array
    int *glCrnDofGroup = new int[numCRNdof];
    for(i=0;i<numCRNdof;++i) glCrnDofGroup[i] = -2;
    paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::makeGlCrnDofGroup, cornerEqs, glCrnDofGroup);
    if(verboseFlag) filePrint(stderr, " ... Locate singularities in Kcc    ...\n");
    int *rank = new int[nGroups];
    int *rows = new int[numCRNdof*nGroups]; // temp fix
    int *rowsptr = new int[nGroups];
    int *nrow = new int[nGroups];
    int *ncol = new int[nGroups];
    int *startcol = new int[nGroups];
    for(i=0;i<nGroups;++i) {
      rank[i] = 0;
      nrow[i] = 0;
      ncol[i] = ngrbmGr[i];
      startcol[i] = 0;
      for(int j=0;j<i;++j) startcol[i] += ngrbmGr[j];
      rowsptr[i] = i*numCRNdof;  // temp fix
    }
    for(i=(numCRNdof-1); i>=0; --i) {
      int group = glCrnDofGroup[i];  // = -2 for augmented dof, -1 for one body only,
                                     // >= 0 if one or more groups is a floating body
      if(group >= -1) {  // dof is not from augmentation
#ifdef DISTRIBUTED
       if(groupProc[group] == this->myCPU) {
#endif
        if(group == -1) group = 0;     // this means there is only one body
        if(rank[group] < ngrbmGr[group]) {  // if true, then haven't found all singularities for this group yet
          double maxval = 0;
          for(int j=0;j<ngrbmGr[group];++j) 
            maxval = myMax(maxval,fabs(globalRcstar[i][startcol[group]+j])); 
          if(maxval > trbm2) {  // ignore zero rows
            rows[rowsptr[group]+nrow[group]] = i;
            nrow[group]++;
            FullM rcsub(globalRcstar, nrow[group], rows+rowsptr[group], ncol[group], startcol[group]);
            FullM Utmp(ncol[group], ncol[group]);
            Utmp.zero();
            int subrank = 0;
            singularValueDecomposition(rcsub,Utmp,ncol[group],nrow[group],subrank,this->fetiInfo->constrain_kcc_tol);
            svdcount++;
            if(subrank == (rank[group]+1)) { // check for singularity
              rank[group]++;
              sing[totRank++] = i;
              if(totRank>=ngrbm) i = -1;  // found all singularities, terminate search
            }
          }
        }
#ifdef DISTRIBUTED
       }
#endif
     }
    }
    delete [] rank;
    delete [] rows;
    delete [] rowsptr;
    delete [] nrow;
    delete [] ncol;
    delete [] startcol;
    delete [] glCrnDofGroup;
  }

#ifdef DISTRIBUTED
  int *sings = new int[ngrbms];
  int *counts = new int[this->numCPUs];
  int *displacements = new int[this->numCPUs];
  this->fetiCom->allGather(&ngrbm, 1, counts, 1);
  if(displacements) displacements[0] = 0;
  for(i=1; i<this->numCPUs; ++i) 
    displacements[i] = counts[i-1] + displacements[i-1];
  this->fetiCom->allGatherv(sing, ngrbm, sings, counts, displacements);
  svdcount = this->fetiCom->globalSum(svdcount);
  totRank = this->fetiCom->globalSum(totRank);
  delete [] sing; sing = sings;
  delete [] counts;
  delete [] displacements;
#endif
  if(verboseFlag) {
    filePrint(stderr, " ... Found %4d singularities       ...\n", totRank);
    filePrint(stderr, " ... Number of SVDs required = %4d ...\n", svdcount);
    filePrint(stderr, " ... Eliminate singularities        ...\n");
  }
  if((totRank < ngrbms) && (cornerEqs->size() > 0)) { cerr << "error: didn't find enough singularities, try reducing FetiInfo::constrain_kcc_tol\n"; exit(-1); } // XXXX
  // condense singular dofs from cornerEqs dsa
  c_cornerEqs = new ConstrainedDSA(*cornerEqs, ngrbms, sing);
  delete [] sing; 
 }
 else {
  c_cornerEqs = new ConstrainedDSA(*cornerEqs, 0);
 }

 if(ngrbms) { // need to do this global sum after constrain_kcc
#ifdef DISTRIBUTED
   delete [] groupProc;
   this->fetiCom->globalSum(nGroups, ngrbmGr);
#endif
   paralApply(this->nsub, this->sd, &BaseSub::setNumGroupRBM, ngrbmGr);
   paralApply(this->nsub, this->sd, &BaseSub::deleteLocalRBMs);
 }

 // end new code ********************************************************

 KccSparse = 0;

 if(c_cornerEqs->size() > 0) {  // replaced cornerEqs with c_cornerEqs

   if(verboseFlag) filePrint(stderr, " ... Assemble Kcc solver            ");

   double tolerance = this->fetiInfo->crbm_tol;  // this is set in input file using global_cor_rbm_tol
   GenBLKSparseMatrix<Scalar> *BLKMatrix = 0;
   GenBlockSky<Scalar> *blockSky = 0;
   GenSkyMatrix<Scalar> *sky = 0;
#ifdef DISTRIBUTED
   int firstAlpha,nRowAlpha;  
   int neq       = c_cornerEqs->size(); 
   int neqPerCPU = neq/this->numCPUs;
   if(this->subToSub->csize() == this->numCPUs) {
     int remainder = neq%this->numCPUs;
     firstAlpha = this->myCPU * neqPerCPU +
                  ((this->myCPU < remainder) ? this->myCPU : remainder);
     nRowAlpha  = neqPerCPU + ((this->myCPU < remainder) ? 1 : 0);
   }
#endif
   this->times.memoryGtGsky -= memoryUsed();
   // build coarse solver
   switch(this->fetiInfo->gtgSolver) {
     default:
     case FetiInfo::skyline: {
#ifdef DISTRIBUTED
       if(this->subToSub->csize() == this->numCPUs && 
          this->fetiInfo->type != FetiInfo::nonlinear &&
          !domain->solInfo().doEigSweep && !domain->solInfo().doFreqSweep) { // 1 subdomain per procesor
         sky = new GenDistSky<Scalar>(coarseConnectivity, c_cornerEqs, tolerance, 
                                      firstAlpha, nRowAlpha);     
         this->times.memoryGtGDelete = 8*sky->size();
       } else
#endif 
       sky = new GenSkyMatrix<Scalar>(coarseConnectivity, c_cornerEqs, tolerance, domain->solInfo().coarseScaled); 
       KccSparse = sky;
       KccSolver = sky;
     }
     break;
     case FetiInfo::pcg:
     case FetiInfo::blocksky: {
#ifdef DISTRIBUTED
       if(this->subToSub->csize() == this->numCPUs && 
          this->fetiInfo->type != FetiInfo::nonlinear && 
          !domain->solInfo().doEigSweep && !domain->solInfo().doFreqSweep) {
         blockSky = new GenDistBlockSky<Scalar>(coarseConnectivity, c_cornerEqs, tolerance,
                                                firstAlpha, nRowAlpha);   
         this->times.memoryGtGDelete = 8*blockSky->size();
       } else
#endif
       blockSky = new GenBlockSky<Scalar>(coarseConnectivity, c_cornerEqs, tolerance);  
       blockSky->zeroAll();
       KccSparse = blockSky;
       KccSolver = blockSky;
     }
     break;
     case FetiInfo::sparse: {
#ifdef DISTRIBUTED
       if(this->subToSub->csize() == this->numCPUs && 
          this->fetiInfo->type != FetiInfo::nonlinear && 
          !domain->solInfo().doEigSweep && !domain->solInfo().doFreqSweep) {
         BLKMatrix = new GenDistBLKSparse<Scalar>(coarseConnectivity, c_cornerEqs, 
                                                  tolerance, firstAlpha, nRowAlpha);   
         this->times.memoryGtGDelete = 8*BLKMatrix->size();
       } else
#endif
if(domain->solInfo().debug_icntl[2] == 1)
       BLKMatrix = new GenBLKSparseMatrix<Scalar>(coarseConnectivity, cornerEqs,
                                                  tolerance, domain->solInfo().sparse_renum, ngrbms);  // PJSA idea: GRBM for Kcc^* ??
else                       
       BLKMatrix = new GenBLKSparseMatrix<Scalar>(coarseConnectivity, c_cornerEqs, 
                                                  tolerance, domain->solInfo().sparse_renum); 
       BLKMatrix->zeroAll();
       KccSparse = BLKMatrix;
       KccSolver = BLKMatrix;
       KccSolver->setPrintNumGrbm(false);
       KccSolver->setPrintNumTrbm(false);
     }
     break;
#ifdef USE_SPOOLES
     case FetiInfo::spooles: {
       GenSpoolesSolver<Scalar> *spsolver = new GenSpoolesSolver<Scalar>(coarseConnectivity, c_cornerEqs);
       KccSparse = spsolver;
       KccSolver = spsolver;
     }
     break;
#endif
#ifdef USE_MUMPS
     case FetiInfo::mumps: {
       // Axel
       GenMumpsSolver<Scalar> *mumpsolver = new GenMumpsSolver<Scalar>(coarseConnectivity, c_cornerEqs, (int *)0, this->fetiCom);
       KccSparse = mumpsolver;
       KccSolver = mumpsolver;
     }
     break;
#endif
   }
   this->times.memoryGtGsky += memoryUsed();
   if(verboseFlag) filePrint(stderr, "."); 

   t5 -= getTime();
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::multKcc);
   t5 += getTime();
   if(verboseFlag) filePrint(stderr, ".");

   t0 -= getTime();
   // assemble the new coarse grid: Kcc -> Kcc - Krc^T Krr^-1 Krc
   for(iSub = 0; iSub < this->nsub; ++iSub) {
     GenAssembledFullM<Scalar> *kel = this->sd[iSub]->getKcc();
     int *dofs = this->sd[iSub]->getKccDofs(c_cornerEqs, glNumCorners, *this->subToEdge, mpcOffset); 
     if(KccSparse) KccSparse->add(*kel, dofs);
     //this->sd[iSub]->deleteKcc();  
   }

   t0 += getTime();
   //this->fetiCom->sync();
   //startTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
   if(verboseFlag) filePrint(stderr, ".\n");

   // Factor coarse solver
#ifdef DISTRIBUTED
   if(verboseFlag) filePrint(stderr, " ... Unify Kcc");
   KccSolver->unify(this->fetiCom);
   if(verboseFlag) filePrint(stderr, "                      ...\n");
#endif
   // HB: HARD CODED FOR SKYLINE
   bool printKcc = false;
   if(printKcc){
     filePrint(stderr, " ... write Kcc skyline in file Kcc.txt ...\n");
//     char MatlabKccFile[]="Kcc.txt";
     GenSkyMatrix<Scalar> *KccMat = dynamic_cast<GenSkyMatrix<Scalar>*>(KccSparse);
//     if(KccMat) KccMat->printMatlab(MatlabKccFile);
     if (KccMat) KccMat->print(stderr);
   }

   KccSolver->setPrintNumTrbm(false);
/*
   int n = KccSolver->neqs();
   long int nz = KccSolver->size();
   double dens = double(2*nz-n)/double(n*n);
   filePrint(stderr, " ... Kcc: n = %d, nz = %d, density = %f    ...\n", n, nz, dens);
   exit(-1); 
*/
   this->fetiCom->sync(); // without this the timings are difficult to interpret
   if(verboseFlag) filePrint(stderr, " ... Factorize Kcc solver");
   startTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
   KccSolver->parallelFactor();
   stopTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
   if(verboseFlag) filePrint(stderr, "           ...\n");

   //stopTimerMemory(this->times.pfactor, this->times.memoryGtGsky);

   if(KccSolver->numRBM() > 0)
      filePrint(stderr, " ... Kcc has %d singularities for tol %e ...\n", KccSolver->numRBM(), this->fetiInfo->crbm_tol);
   if(computeRbms && this->myCPU == 0) {
     if(!this->fetiInfo->constrain_kcc && KccSolver->numRBM() != ngrbms && domain->probType() != SolverInfo::Modal) {
       cerr << " *** WARNING: number of singularities in Kcc does not match the number of Geometric RBMs of the decomposed domain \n";
       cerr << " *** try adjusting global_cor_rbm_tol or set constrain_kcc on \n";
     }
   } 

 } else
   KccSolver = 0;

  if(coarseConnectivity) delete coarseConnectivity;
  stopTimerMemory(this->times.coarse1, this->times.memoryGtG);
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::KrrReSolve(int iSub, GenDistrVector<Scalar> &ur)
{
 this->sd[iSub]->Krr->reSolve(ur.subData(iSub));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::extractFr(int iSub, GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fr)
{
 this->sd[iSub]->getFr(f.subData(iSub), fr.subData(iSub));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::extractFc(int iSub, GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fc)
{
 this->sd[iSub]->getFc(f.subData(iSub), fc.subData(iSub));
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::extractFw(int iSub, GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fw)
{
 this->sd[iSub]->getFw(f.subData(iSub), fw.subData(iSub));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::getFc(GenDistrVector<Scalar> &f, GenVector<Scalar> &fc)
{
 GenDistrVector<Scalar> distFc(internalC);
 distFc.zero(); 

 execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::extractFc, f, distFc);

 // Assemble fc (force at the corners)
 fc.zero();
 for(int iSub=0; iSub<this->nsub; ++iSub) {
    int *subdCDofs = this->sd[iSub]->cornerEqNums;
    Scalar *subFc = distFc.subData(iSub);
    for(int iDof = 0; iDof < this->sd[iSub]->numCoarseDofs(); ++iDof)
      if(subdCDofs[iDof] > -1)
        fc[subdCDofs[iDof] ] += subFc[iDof];
 }
}

template<class Scalar>
bool
GenFetiDPSolver<Scalar>::updateActiveSet(GenDistrVector<Scalar> &v, int flag)
{
  // flag = 0 to add all constraints to active set if v_i > tol
  // flag = 1 to remove all constraints from active set if v_i > tol
  bool status_change = false;
  double tol = 0.0;
  //if(flag == 1) return updateActiveSet_one(v, tol, flag); // XXXX update one by one
  if(flag == 0 && ngrbms) paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::saveMpcStatus2);
  execParal4R(this->nsub, this, &GenFetiDPSolver<Scalar>::subUpdateActiveSet, v, tol, flag, status_change);
#ifdef DISTRIBUTED
  status_change = this->fetiCom->globalMax((int) status_change);
#endif
  if(status_change) {
    paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::sendMpcStatus, mpcPat, flag);
    mpcPat->exchange();
    paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::recvMpcStatus, mpcPat, flag);
    if(ngrbms && this->fetiInfo->project_g) rebuildGtGtilda();

    // enforce constraint qualification ... remove redundant inequality constraints from active set
    // notes: this is only necessary when there are equality constraints (ngrbms > 0) because the inequality constraints alone are always linearly independent
    //      : redundant inequality constraints still need to be enforced and this is not done unless equality constraints are projected (project_g == true)
    if(this->fetiInfo->cq_type == FetiInfo::crcq && flag == 0 && ngrbms && this->fetiInfo->project_g) {
      if(redundant()) { // eliminate redunant constraints by restoring active set and GtG to previous state then check/update one by one
        paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::restoreMpcStatus2);
        rebuildGtGtilda();
        status_change = false;
        while(updateActiveSet_one(v, tol, flag)) status_change = true; // update one by one
      }
    }
  }
  if(status_change) { 
    if(flag == 0) dualStatusChange = true; else primalStatusChange = true;
    if(this->fetiInfo->contactPrintFlag && this->myCPU == 0) cerr << " ";
  }
  return status_change;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subUpdateActiveSet(int iSub, GenDistrVector<Scalar> &lambda, double tol, int flag, bool &statusChange)
{
  this->sd[iSub]->updateActiveSet(lambda.subData(this->sd[iSub]->localSubNum()), tol, flag, statusChange);
}

template<class Scalar>
bool
GenFetiDPSolver<Scalar>::updateActiveSet_one(GenDistrVector<Scalar> &v, Scalar tol, int flag)
{
  // flag = 0: add one constraint to active set corresponding to max v[i] > tol
  // flag = +1: remove one constraint from active set corresponding to max v[i] > tol
  bool status_change = false;
  int mpcid = selectOne(v, tol, flag); 
  if(mpcid > -1) {
    paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::updateActiveSet_one, mpcid, flag); 
    if(ngrbms && this->fetiInfo->project_g) rebuildGtGtilda();
    status_change = true;
    if(this->fetiInfo->contactPrintFlag && this->myCPU == 0) cerr << " ";
  }
  return status_change;
}

template<class Scalar>
int
GenFetiDPSolver<Scalar>::selectOne(GenDistrVector<Scalar> &v, Scalar tol, int flag)
{
  // find maximum v[i] > tol and return corresponding mpcid to either add to active set (flag = 0) or remove from active set (flag = 1)

  // 1. find the maximum v[i] > tol for all the inactive inequality constraints
  //XXXX Scalar max = Max(v, tol, flag);
  Scalar max = Max(v, flag);
  if(!ScalarTypes::greaterThan(max, tol)) return -1;

  // 2. find the mpc with the largest global id which has this "largest positive value"
  int mpcid = Equal(v, max, flag);

  // 3. enforce linear independence constraint qualification
  if(this->fetiInfo->cq_type == FetiInfo::crcq && flag == 0 && ngrbms && this->fetiInfo->project_g) { 
    GenDistrVector<Scalar> n_u(this->interface);
    execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::getn_u, n_u, mpcid);
    project(n_u,n_u,3);
    if(n_u.sqNorm() <= this->fetiInfo->cq_tol) {
      paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::markRedundant, mpcid); // flag mpcid so it cannot be added to active set until after next primal status change
      return selectOne(v, tol, flag); // pick another mpc (recursive)
    }
  }
  return mpcid;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::getn_u(int iSub, GenDistrVector<Scalar> &n_u, int mpcid)
{
  this->sd[iSub]->getn_u(n_u.subData(this->sd[iSub]->localSubNum()), mpcid);
}

template<class Scalar>
Scalar
GenFetiDPSolver<Scalar>::Max(GenDistrVector<Scalar> &v, int flag)
{
  // find maximum v_i in subset defined by flag
  Scalar *partial = (Scalar *) alloca(sizeof(Scalar) * this->nsub);
  execParal3R(this->nsub, this, &GenFetiDPSolver<Scalar>::subMax, v, partial, flag);
  Scalar max = (this->nsub > 0) ? partial[0] : -numeric_limits<double>::max();
  for(int i = 1; i < this->nsub; i++)
    if(ScalarTypes::greaterThan(partial[i], max)) max = partial[i];
#ifdef DISTRIBUTED
  max = this->fetiCom->globalMax(max);
#endif
  return max; // note: max = tol if no v[i] > tol found in subset
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subMax(int iSub, GenDistrVector<Scalar> &v, Scalar *partial, int flag)
{
  this->sd[iSub]->Max(v.subData(this->sd[iSub]->localSubNum()), partial[iSub], flag);
}

template<class Scalar>
int
GenFetiDPSolver<Scalar>::Equal(GenDistrVector<Scalar> &v, Scalar val, int flag)
{
  // find unique mpcid with v_i = val in subset defined by flag
  int *ipartial = (int *) alloca(sizeof(int) * this->nsub);
  //execParal4R(this->nsub, this, &GenFetiDPSolver<Scalar>::subEqual, v, val, ipartial, flag);
  execParal4R<GenFetiDPSolver<Scalar>, GenDistrVector<Scalar>, Scalar, int, int>
          (this->nsub, this, &GenFetiDPSolver<Scalar>::subEqual, v, val, ipartial, flag);
  int mpcid = -1;
  for(int i = 0; i < this->nsub; i++)
    if(ipartial[i] > mpcid) mpcid = ipartial[i];
#ifdef DISTRIBUTED
  mpcid = this->fetiCom->globalMax(mpcid);
#endif
  return mpcid; // note return -1 if no v[i] = val found in subset
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subEqual(int iSub, GenDistrVector<Scalar> &v, Scalar val, int *ipartial, int flag)
{
  this->sd[iSub]->Equal(v.subData(this->sd[iSub]->localSubNum()), val, ipartial[iSub], flag);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::update(Scalar nu, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &p, 
                                GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fp, 
                                GenDistrVector<Scalar> &ur, GenDistrVector<Scalar> &dur, 
                                GenVector<Scalar> &uc, GenVector<Scalar> &duc, int l, double alphabar, double rho)
{
  // update dual (lambda)
  if(globalFlagCtc) expansionStep(lambda, nu, p, alphabar); // p = p_r, lambda += nu*p_r, Fp = F*p_r 
  else lambda.linAdd(nu,p);

  // Update residual (r)
/* XXXX compute r directly from lambda, may be more accurate due to roundoff 
  if(dualStatusChange) { 
    GenDistrVector<Scalar> &fr      = this->wksp->ret_fr(); 
    GenDistrVector<Scalar> &fw      = this->wksp->ret_fw(); 
    GenVector<Scalar> &fc           = this->wksp->ret_fc(); 
    localSolveAndJump(fr, lambda, ur, fc, uc, r, fw, rho, mu); 
    if(this->fetiInfo->linesearch) { // compute F*p_r
      GenDistrVector<Scalar> &r_k = this->wksp->ret_r_copy();
      Fp.linC(1.0/nu,r,-1.0/nu,r_k);
    }
  }
  else r.linAdd(nu,Fp);
*/
  if(dualStatusChange) { localSolveAndJump(p, dur, duc, Fp, rho);  nExtraFp++; }
  r.linAdd(nu,Fp); // note: r += nu*Fp is not thread-safe

  // check wolfe conditions, adjust step length and re-update if lagrangian is not sufficiently reduced
  if(this->fetiInfo->linesearch && linesearch(nu, p, l, alphabar)) 
    update(nu, lambda, p, r, Fp, ur, dur, uc, duc, l, alphabar, rho); // recursive 

  // Update primal (ur, uc)
  else { ur.linAdd(nu,dur); uc += nu*duc; }
}

template<class Scalar>
Scalar
GenFetiDPSolver<Scalar>::alpha_f()
{
  // returns the maximum feasible step length in direction p
  GenDistrVector<Scalar> &q = this->wksp->ret_q(), &lambda = this->wksp->ret_lambda(), &p = this->wksp->ret_p();
  if(newton_iter > 0) lambda += (*lambda_total);
  quotient(q, lambda, p);
  if(newton_iter > 0) lambda -= (*lambda_total);
  return Max(q, 0);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::expansionStep(GenDistrVector<Scalar> &lambda, Scalar nu, GenDistrVector<Scalar> &p, double alphabar)
{
  // step(nu,p_r): lambda^(k+1) = lambda^(k) + nu*p_r ... p_r is the reduced search direction which is also computed
  // notes: if step(nu,p) is feasible then p_r = p and dualStatusChange = false
  //        nu should always be negative since CG search direction p is a descent direction (hence r^T*p < 0) and F is positive (hence  p^T*Fp > 0)
  if(!this->fetiInfo->linesearch) this->wksp->save_lambda();
  Scalar nu1, nu2 = nu;

  // first 1/2 step
  bool expansion = false;
  //if(this->fetiInfo->expansion == 2 && ScalarTypes::greaterThan(nu1 = alpha_f()-this->fetiInfo->expansion_tol, nu)) {
  //if(this->fetiInfo->expansion == 2 && ScalarTypes::greaterThan(nu1 = alpha_f()*(1.0+this->fetiInfo->expansion_tol), nu)) {
  if(this->fetiInfo->expansion == 2) { 
    Scalar alphaf = alpha_f(); 
    nu1 = (ScalarTypes::Real(alphaf) == 0.0) ? -this->fetiInfo->expansion_tol : alpha_f()*(1.0+this->fetiInfo->expansion_tol);
    if(ScalarTypes::greaterThan(nu1, nu)) {
      //if(this->myCPU == 0) cerr << "expansion first half-step, nu = " << nu << ", alpha_f = " << nu1+this->fetiInfo->expansion_tol
      //                          << ", expansion_tol = " << this->fetiInfo->expansion_tol << ", nu1 = " << nu1 << endl;
      lambda.linAdd(nu1,p); 
      if(newton_iter > 0) lambda += (*lambda_total);
      if(updateActiveSet(lambda, 0)) { project(lambda, lambda, -1, true); }
      if(newton_iter > 0) lambda -= (*lambda_total);
      GenDistrVector<Scalar> &r = this->wksp->ret_r(), &Fp = this->wksp->ret_Fp();
      p.linC(1.0,r,nu1,Fp);
      project(p,p); // p^(k+1/2) = P^T*r^(k+1/2) = P^T*(r^(k)+alpha_f*Fp^(k))
      nu2 = -alphabar;
      expansion = true;
    }
  }

  // second 1/2 step
  //if(this->myCPU == 0 && expansion) cerr << "expansion second half-step, alphabar = " << alphabar << ", nu2 = " << nu2 
  //                                       << ", dualStatusChange after 1st half-step = " << dualStatusChange << endl;
  lambda.linAdd(nu2,p); 
  if(newton_iter > 0) lambda += (*lambda_total);
  if(updateActiveSet(lambda, 0)) { project(lambda, lambda, -1, true); } 
  if(newton_iter > 0) lambda -= (*lambda_total);
  if(expansion && !dualStatusChange) { dualStatusChange = true; if(this->myCPU == 0) cerr << "WARNING: no status change in expansion step! this shouldn't happen. Try increasing expansion_tol\n"; }

  // compute the reduced search direction: p_r = (lambda^(k+1) - lambda^(k))/nu and F*p_r
  if(dualStatusChange) {
    GenDistrVector<Scalar> &lambda_k = this->wksp->ret_lambda_copy();
    p.linC(1.0/nu,lambda,-1.0/nu,lambda_k);
  }
}

template<class Scalar>
bool
GenFetiDPSolver<Scalar>::linesearch(Scalar &nu, GenDistrVector<Scalar> &p, int &l, double &alphabar)
{
  // one linesearch iteration to find a better step length nu
  bool wolfe1, wolfe2;
  if(!checkWolfe(nu, p, this->fetiInfo->linesearch, wolfe1, wolfe2) && l++ < this->fetiInfo->linesearch_maxit) {
    restoreStep();
    Scalar nu_cg; if(this->fetiInfo->expansion == 2) { nu_cg = nu; nu = -alphabar; } // do linesearch on alphabar for 2-step expansion
    switch(this->fetiInfo->linesearch) {
      case 1 : { // backtracking linesearch iteration (armijo-goldstein)
        if(!wolfe1) nu *= this->fetiInfo->linesearch_tau; // default is 0.8
      } break;
      case 2 : { // bisection linesearch iteration based on weak wolfe conditions
        if(!wolfe1) { beta_l = nu; nu = 0.5*(alpha_l+beta_l); }
        else if(!wolfe2) { alpha_l = nu; nu = (beta_l == -1.0e256) ? 2.0*alpha_l : 0.5*(alpha_l+beta_l); }
      } break;
      case 3 : { // more-thuente linesearch based on strong wolfe conditions
        if(this->myCPU == 0) cerr << " *** ERROR: linesearch type 3 is not implemented\n";
        exit(-1);
      } break;
    }
    if(this->fetiInfo->contactPrintFlag >= 2 && this->myCPU == 0)
      cerr << " linesearch (type " << this->fetiInfo->linesearch << "): l = " << l << ", nu = " << nu << endl;
    nLinesearchIter++;
    stepLengthChange = true; // always need to reset orthoset if nu != nu_cg
    if(this->fetiInfo->expansion == 2) { alphabar = -ScalarTypes::Real(nu); nu = nu_cg; }
    return true;
  }
  else {
    if(l >= this->fetiInfo->linesearch_maxit && this->myCPU == 0) cerr << "linesearch did not converge in " << this->fetiInfo->linesearch_maxit << " iterations\n";
    return false;
  }
}

template<class Scalar>
bool
GenFetiDPSolver<Scalar>::checkWolfe(Scalar nu, GenDistrVector<Scalar> &p, int which, bool &wolfe1, bool &wolfe2)
{
  // check if wolfe conditions are satisfied for step(nu,p)
  // notes: 0 < c1 < c2 < 1  (defaults c1 = 1e-3, c2 = 0.9)
  //        which = 1: check first wolfe condition, which = 2: check first and second weak, which = 3: check first and second strong
  Scalar rp0, rp1, rp2, pFp;

  // check 1st wolfe condition
  if(which >= 1) {  
    GenDistrVector<Scalar> &p_copy = this->wksp->ret_p_copy(), &r_copy = this->wksp->ret_r_copy() /* r^(k) */, &Fp = this->wksp->ret_Fp();
    rp0 = r_copy*p_copy;
    rp1 = r_copy*p;
    pFp = p*Fp;
    wolfe1 = (ScalarTypes::Real(nu*nu/2.0*pFp + nu*rp1) <= this->fetiInfo->wolfe_c1*ScalarTypes::Real(nu*rp0));
    //if(!wolfe1 && this->myCPU == 0) cerr << "pFp = " << pFp << ", delta_lag = " << nu*nu/2.0*pFp + nu*rp1 << endl; // XXXX DEBUG
  }
  else wolfe1 = true;

  // check 2nd wolfe condition
  if(which >= 2) {
    GenDistrVector<Scalar> &r = this->wksp->ret_r(); // r^(k+1)
    rp2 = r*p;
    wolfe2 = (which == 2) ? (ScalarTypes::Real(-rp2) >= this->fetiInfo->wolfe_c2*ScalarTypes::Real(-rp0))  // weak
                          : (ScalarTypes::norm(-rp2) <= this->fetiInfo->wolfe_c2*ScalarTypes::norm(-rp0)); // strong
  }
  else wolfe2 = true;

  return (wolfe1 && wolfe2);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::saveStep()
{
  this->wksp->save();
  if(globalFlagCtc) paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::saveMpcStatus);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::restoreStep()
{
  if(dualStatusChange) { 
    this->wksp->restore(true);
    paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::restoreMpcStatus);
    if(ngrbms && this->fetiInfo->project_g) rebuildGtGtilda();
    dualStatusChange = false;
  }
  else this->wksp->restore(false);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::quotient(GenDistrVector<Scalar> &q, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &p)
{
  execParal3R(this->nsub, this, &GenFetiDPSolver<Scalar>::subQuotient, q, lambda, p);
}

template<class Scalar>
bool
GenFetiDPSolver<Scalar>::feasible(GenDistrVector<Scalar> &lambda, bool check_ieq, bool print_warning_eq)
{
  if(this->glNumMpc == 0) return true;
  // note: use fetiInfo->iequ_tol = 0.0 to check true feasibility of dual inequalities
  bool eq = true, ieq = true;
  double eq_error = 0.0, ieq_error = 0.0;
  if(newton_iter > 0) lambda += (*lambda_total);
  if(globalFlagCtc && check_ieq) {
    // check N^T*lambda <= 0
    GenDistrVector<Scalar> lambda_c(this->interface);
    execParal4R(this->nsub, this, &GenFetiDPSolver<Scalar>::chop, lambda, lambda_c, 0.0, 1);
    ieq = ((ieq_error = sqrt(lambda_c.sqNorm())) <= this->fetiInfo->iequ_tol);
  }
  if(ngrbms) {
    // check G^T*lambda+e = 0
    GenVector<Scalar> &e = this->wksp->ret_e();
    GenVector<Scalar> &gamma = this->wksp->ret_gamma();
    int flag = (this->fetiInfo->outerloop == FetiInfo::CGAL) ? -1 : 0;
    trMultG(lambda, gamma, 1.0, 0.0, flag); // gamma = G^T*lambda
    gamma += e;
    eq = ((eq_error = sqrt(gamma.sqNorm())) <= this->fetiInfo->equi_tol);
  }
  if(newton_iter > 0) lambda -= (*lambda_total);
  if(this->myCPU == 0) {
    if(!eq && print_warning_eq) 
      cerr << "warning: lambda is not feasible wrt equality constraints (error = " << eq_error << ", equi_tol = " << this->fetiInfo->equi_tol <<")\n";
    if(!ieq) 
      cerr << "warning: lambda is not feasible wrt inequality constraints, (error = " << ieq_error << ", iequ_tol = " << this->fetiInfo->iequ_tol << ")\n";
  }
  return (eq && ieq);
}

// declarations of non-templated functions, see FetiDPCore.C for implementation (double only)
template<>
void
GenFetiDPSolver<DComplex>::split(int, GenDistrVector<DComplex> &, GenDistrVector<DComplex> &, GenDistrVector<DComplex> &, GenDistrVector<DComplex> &);

template<>
void
GenFetiDPSolver<double>::split(int, GenDistrVector<double> &, GenDistrVector<double> &, GenDistrVector<double> &, GenDistrVector<double> &);

template<>
void
GenFetiDPSolver<DComplex>::chop(int, GenDistrVector<DComplex> &, GenDistrVector<DComplex> &, double, int);

template<>
void
GenFetiDPSolver<double>::chop(int, GenDistrVector<double> &, GenDistrVector<double> &, double, int);

template<>
void
GenFetiDPSolver<DComplex>::subQuotient(int, GenDistrVector<DComplex> &, GenDistrVector<DComplex> &, GenDistrVector<DComplex> &);

template<>
void
GenFetiDPSolver<double>::subQuotient(int, GenDistrVector<double> &, GenDistrVector<double> &, GenDistrVector<double> &);

template<> 
void
GenFetiDPSolver<double>::addRstar_gT(int, GenDistrVector<double> &, GenVector<double> &);

template<> 
void
GenFetiDPSolver<double>::subtractRstar_g(int, GenDistrVector<double> &, GenVector<double> &);

template<>
void
GenFetiDPSolver<DComplex>::computeProjectedDisplacement(GenDistrVector<DComplex> &);

template<>
void
GenFetiDPSolver<double>::computeProjectedDisplacement(GenDistrVector<double> &);

template<> 
void
GenFetiDPSolver<double>::assembleE(int, GenVector<double> &e, GenDistrVector<double> &);

template<>
bool
GenFetiDPSolver<DComplex>::inconsistent(GenVector<DComplex> &, bool);

template<>
bool
GenFetiDPSolver<double>::inconsistent(GenVector<double> &, bool);

template<>
void
GenFetiDPSolver<DComplex>::makeE(GenDistrVector<DComplex> &);

template<>
void
GenFetiDPSolver<double>::makeE(GenDistrVector<double> &);

template<>
void
GenFetiDPSolver<DComplex>::makeGtG();

template<> 
void
GenFetiDPSolver<double>::makeGtG();

template<>
void
GenFetiDPSolver<DComplex>::getRBMs(DComplex *globRBM);
                                                                                                                                             
template<>
void
GenFetiDPSolver<double>::getRBMs(double *globRBM);

template<>
void
GenFetiDPSolver<DComplex>::getRBMs(GenDistrVectorSet<DComplex> &globRBM);
                                                                                                                                             
template<>
void
GenFetiDPSolver<double>::getRBMs(GenDistrVectorSet<double> &globRBM);

template<>
void
GenFetiDPSolver<DComplex>::getGlobalRBM(int iSub, int &iRBM, GenDistrVector<DComplex> &R);

template<>
void
GenFetiDPSolver<double>::getGlobalRBM(int iSub, int &iRBM, GenDistrVector<double> &R);

template<>
void
GenFetiDPSolver<DComplex>::addRalpha(int, GenDistrVector<DComplex> &, GenVector<DComplex> &);

template<> 
void 
GenFetiDPSolver<double>::addRalpha(int, GenDistrVector<double> &, GenVector<double> &);

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::solve(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u)
{
 if (verboseFlag) filePrint(stderr,"\n ... Begin FETI-DP Solve            ...\n");

 init_solve();
 if(this->fetiInfo->type == FetiInfo::nonlinear) { // PJSA 9-18-2007
   if(newton_iter == 0 && this->glNumMpc > 0) {
     //cerr << "zeroing mpc forces \n";
     paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::zeroMpcForces);
   }
 }

 switch(this->fetiInfo->outerloop) {
   default:
   case 0:
       if(verboseFlag) {
         filePrint (stderr, " ... CG solver selected             ...\n\n");
       }
       if(this->myCPU == 0 && typeid(Scalar)==typeid(DComplex) && !this->fetiInfo->complex_hermitian)
         cerr << " *** WARNING: CG is not valid for all complex symmetric matrices even if positive definite.\n"
              << " *** If your matrix is complex hermitian include \"outerloop CG hermitian\" under FETI in your input file.\n"
              << " *** Otherwise, \"outerloop GMRES\" or \"outerloop GCR\" are safer choices.\n";
       if(verboseFlag) {
         filePrint (stderr, "----------------------------------------------\n");
         filePrint (stderr, "          Iterations loop monitoring          \n");
         filePrint (stderr, "----------------------------------------------\n");
       }
       solveCG(f, u);
       break;
   case 1:
       if(verboseFlag) {
         filePrint (stderr, " ... GMRES solver selected          ...\n");
         filePrint (stderr, "----------------------------------------------\n");
         filePrint (stderr, "          Iterations loop monitoring          \n");
         filePrint (stderr, "----------------------------------------------\n");
       }
       solveGMRES(f, u);
       break;
   case 2:
       if(verboseFlag) {
         filePrint (stderr, " ... GCR solver selected            ...\n\n");
         filePrint (stderr, "----------------------------------------------\n");
         filePrint (stderr, "          Iterations loop monitoring          \n");
         filePrint (stderr, "----------------------------------------------\n");
       }
       solveGCR(f, u);
       break;
   case 3:
       filePrint (stderr, " ... CG Augmented Lagrangian solver selected ...\n\n");
       filePrint (stderr, "----------------------------------------------\n");
       filePrint (stderr, "          Iterations loop monitoring          \n");
       filePrint (stderr, "----------------------------------------------\n");
       solveCG_augLag(f,u);
       break;
 }
 if(this->fetiInfo->type == FetiInfo::nonlinear) { // PJSA 9-18-2007
   if(newton_iter == 0 && this->glNumMpc > 0) { 
     if(lambda_total) delete lambda_total; 
     lambda_total = new GenDistrVector<Scalar>(this->interface);
     lambda_total->zero();
   }
   newton_iter++;
   if(this->glNumMpc > 0) {
     GenDistrVector<Scalar> &lambda = this->wksp->ret_lambda();
     *lambda_total += lambda;
   }
 }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::init_solve()
{
  // include here anything that needs to be initialized before each solve
  GenDistrVector<Scalar> &deltaU  = this->wksp->ret_deltaU(); 
  GenDistrVector<Scalar> &deltaF  = this->wksp->ret_deltaF();
  GenVector<Scalar> &alpha        = this->wksp->ret_alpha();
  deltaU.zero();
  deltaF.zero();
  alpha.zero();
  dualStatusChange = primalStatusChange = stepLengthChange = false;
  nExtraFp = nRebuildGtG = nRebuildCCt = nLinesearchIter = nSubIterDual = nSubIterPrimal = nStatChDual = nStatChPrimal = 0;
  if(globalFlagCtc && this->fetiInfo->linesearch == 0) { // activate linesearch if monotonicity not guaranteed
    if((this->fetiInfo->expansion == 1) || (this->fetiInfo->alphabar_cntl > 2.0) 
       || (this->fetiInfo->outerloop == FetiInfo::CGAL && this->fetiInfo->cgal_prec)  /* this shouldn't be necessary, just checking */
       || (this->fetiInfo->outerloop == FetiInfo::CG && (this->fetiInfo->precno > 0 || this->fetiInfo->mpc_precno > 0))) /* this shouldn't be necessary, just checking */ 
    {
      this->fetiInfo->linesearch = 1;
      if(this->myCPU == 0) cerr << " ... Activating linesearch          ... \n";
    }
  }
  if(globalFlagCtc && this->numSystems > 0) this->resetOrthoSet();
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::init_iter()
{
  // include here anything that needs to be initialized before each iteration
  dbg_alloca(0);
  if(dualStatusChange) nStatChDual++;
  if(primalStatusChange) { nStatChPrimal++; paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::unmarkRedundant); }
  if((dualStatusChange || primalStatusChange) && this->fetiInfo->contactPrintFlag && this->myCPU == 0) cerr << endl;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::init_linesearch()
{
  saveStep(); alpha_l = 0.0; beta_l = -1.0e256; stepLengthChange = false;
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::solveCG(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u)
{
 t7 = -getTime(); this->times.solve -= getTime();

 GenDistrVector<Scalar> &fr      = this->wksp->ret_fr();  
 GenDistrVector<Scalar> &ur      = this->wksp->ret_ur(); 
 GenDistrVector<Scalar> &dur     = this->wksp->ret_du(); 
 GenDistrVector<Scalar> &lambda  = this->wksp->ret_lambda();       // lagrange multipliers
 GenDistrVector<Scalar> &r       = this->wksp->ret_r();            // residual
 GenDistrVector<Scalar> &w       = this->wksp->ret_w();           // projected residual
 GenDistrVector<Scalar> &y       = this->wksp->ret_y();           // re-projected residual
 GenDistrVector<Scalar> &z       = this->wksp->ret_z();           // preconditioned residual
 GenDistrVector<Scalar> &p       = this->wksp->ret_p();           // search direction
 GenDistrVector<Scalar> &Fp      = this->wksp->ret_Fp(); 
 GenDistrVector<Scalar> &deltaU  = this->wksp->ret_deltaU(); 
 //GenDistrVector<Scalar> &deltaF  = this->wksp->ret_deltaF(); 
 GenDistrVector<Scalar> &fw      = this->wksp->ret_fw(); 
 GenVector<Scalar> &fc           = this->wksp->ret_fc(); 
 GenVector<Scalar> &uc           = this->wksp->ret_uc();  
 GenVector<Scalar> &duc          = this->wksp->ret_duc(); 

 int iter, l = 0;
 double ff, dd, relChError, lastError=0.0, error, dual_error, alphabar;    //CRW
 Scalar pFp, nu;
 if(this->glNumMpc == 0) this->fetiInfo->project_g = false;

 //cerr << "nonlinear = " << (this->fetiInfo->type == FetiInfo::nonlinear) << ", numSystems = " << this->numSystems
  //    << ", newton_iter = " << newton_iter << ", oSetCG->numDir() = " << this->oSetCG->numDir()
   //   << ", oSetCG->numOrthoSets() = " << this->oSetCG->numOrthoSets() << endl; // DEBUG NONLINEAR
#ifndef SALINAS
 if(this->fetiInfo->type == FetiInfo::nonlinear && /*this->numSystems > 0 &&*/ this->glNumMpc) { // DEBUG NONLINEAR
   paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::sendMpcRhs, mpcSPat);
   mpcSPat->exchange();
   paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::recvMpcRhs, mpcSPat);
 }
 if(domain->solInfo().debug_icntl[5] == 0) { // DEBUG NONLINEAR: rhs should be zero if first newton iter solved exactly
   if(this->fetiInfo->type == FetiInfo::nonlinear && this->numSystems > 0 && this->glNumMpc) paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::zeroFreeMpcRhs);
 }
#endif

 
 // extract fr, fc and fw from f
 ff = extractForceVectors(f, fr, fc, fw);

 if(globalFlagCtc && this->fetiInfo->expansion == 2) alphabar = this->fetiInfo->alphabar_cntl/computeFNorm();
 
 // Compute initial lagrange multipliers: lambda^0 = -G*(G^T*G)^-1 * e, where e = R^T*f ... also gamma = G^T*lambda+e
 // also for feti-dpc: expand active set and re-compute initial lagrange multipliers if lambda^0 is not feasible
 computeL0(lambda, f); 

 // Compute initial residual: r^0 = F*lambda^0 + d ... also uc = Kcc^-1(fc+Krr^-1*Krc^T*lambda), ur = Krr^-1*(fr-Br^T*lambda-Krc*uc)
 localSolveAndJump(fr, lambda, ur, fc, uc, r, fw); 

 // Compute initial projected residual and dual error: w^0 = P^T * r^0
 // also for feti-dpc: contract active set if w^0 is not proportional (primal planing)
 tProject(r, w, dual_error); 

 // Save and print initial dual error
 dd = dual_error;
 if(verboseFlag) filePrint(stderr," ... Initial residual norm %e\n", sqrt(dd));
 if(dd == 0.0) { mergeSolution(ur, uc, u, lambda); return; }

 // Multiple rhs prediction (note: not used for contact)
 if(predict(w, lambda)) {
   localSolveAndJump(fr, lambda, ur, fc, uc, r, fw);
   tProject(r, w, dual_error); 
   if(verboseFlag) filePrint(stderr," ... Initial residual norm after MRHS prediction %e\n", sqrt(dual_error));
 }

 if(verboseFlag) filePrint(stderr," Iteration  Relative Primal Error  Relative Dual Error\n");

 for(iter = 0; iter < this->maxiter; ++iter, ++iterTotal) {
   init_iter(); 

   // Precondition: z = M^-1 * w
   error = preCondition(w, z); 

   // Print errors
   if((this->fetiInfo->numPrint() > 0) && (iter % this->fetiInfo->numPrint() == 0) && verboseFlag)
     filePrint(stderr," %4d %23.6e %21.6e\n", iter, sqrt(error/ff), sqrt(dual_error/dd));

   // Test for convergence or maximum number of iterations
   if((fabs(error) < this->epsilon2 * ff) || (iter == this->maxiter-1)) {
     if(globalFlagCtc && domain->solInfo().debug_cntl[9] > this->fetiInfo->iequ_tol) domain->solInfo().debug_icntl[9] = 1; else { // XXXX make sure primal inequalities are feasible too!
       if(!((this->fetiInfo->numPrint() > 0) && (iter % this->fetiInfo->numPrint() == 0)) && verboseFlag)
         filePrint(stderr," %4d %23.6e %21.6e\n", iter, sqrt(error/ff), sqrt(dual_error/dd)); 
       this->times.iterations[this->numSystems].stagnated = 0;
     }
     break;
   }

   // Test for stagnation
   if(iter > 1 && (relChError = DABS((error-lastError)/lastError)) < this->fetiInfo->stagnation_tol) {
     filePrint(stderr, "STAGNATION: Relative change in primal error = %e\n", relChError);
     if(sqrt(error) < this->fetiInfo->absolute_tol) this->times.iterations[this->numSystems].stagnated = 0; else // XXXX continue if absolute error is sufficiently small 
     { this->times.setStagnate(this->numSystems); this->resetOrthoSet(); }
     break;
   } else lastError = error;

   // Krylov acceleration
   if(this->fetiInfo->nlPrecFlg) nlPreCondition(w, z);

   // Re-project: y = P * z
   project(z, y); 

   // Search direction
   orthogonalize(y, p);
   
   // Matrix vector product
   localSolveAndJump(p, dur, duc, Fp);

   // Step length
   if(!this->fetiInfo->complex_hermitian) // Note: CG by default uses the indefinite vector inner product i.e. x^y = x^T y
    { pFp = p^Fp; nu = -(w^p)/pFp; }      // CG in this form can be used for complex symmetric matrices but is NOT valid for
   else                                   // all such matrices, so we don't recommed using it. Try GCR or GMRES instead!!!
    { pFp = p*Fp; nu = -((w*p) / pFp); }  // CG is valid for positive definite Hermitian matrices but to set this up you need 
                                          // to use to the standard inner product x*y = y^H x for pFp, nu and beta (in Feti.d/CGOrthoset.C)
                                          // this option hasn't really been tested since our matrices are typically NOT Hermitian

   // Update solution: lambda += nu*p, r += nu*Fp, ur += nu*dur, uc += nu*duc 
   // optional linesearch to adjust step length nu if wolfe-armijo conditions are not satisfied
   // also for feti-dpc: reduce step and expand active set if lambda is not feasible (dual planing)
   if(globalFlagCtc) dualStatusChange = primalStatusChange = stepLengthChange = false;
   if(this->fetiInfo->linesearch) { l = 0; init_linesearch(); }
   update(nu, lambda, p, r, Fp, ur, dur, uc, duc, l, alphabar);

   // Project: w = P^T * r 
   // also for feti-dpc: contract active set if w is not proportional (primal planing)
   tProject(r, w, dual_error); 

   // add search direction to orthoset or reset if necessary
   orthoAdd(p, Fp, pFp);
 }

 // Assemble and store primal solution u
 ur += deltaU; // make solution compatible ur += deltaU
 mergeSolution(ur, uc, u, lambda);

 // Store number of iterations, errors, timings and memory used
 this->setAndStoreInfo(iter+1, (error/ff), sqrt(dual_error/dd));
 if(this->numSystems == 1) this->times.memoryFETI += memoryUsed();
 this->times.solve += getTime(); t7 += getTime();
 this->times.iterations[this->numSystems-1].cpuTime = this->times.solve;
 printSummary(iter);
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::solveCG_augLag(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u)
{
  t7 = -getTime(); this->times.solve -= getTime();
  if(this->oSetCG->numDir() > 0) this->oSetCG->reset(); // multiple rhs is not supported

  GenDistrVector<Scalar> &fr      = this->wksp->ret_fr();  
  GenDistrVector<Scalar> &ur      = this->wksp->ret_ur(); 
  GenDistrVector<Scalar> &dur     = this->wksp->ret_du(); 
  GenDistrVector<Scalar> &lambda  = this->wksp->ret_lambda();      // Lagrange multipliers
  GenDistrVector<Scalar> &r       = this->wksp->ret_r();           // residual
  GenDistrVector<Scalar> &w       = this->wksp->ret_w();           // projected residual
  GenDistrVector<Scalar> &y       = this->wksp->ret_y();           // re-projected residual
  GenDistrVector<Scalar> &z       = this->wksp->ret_z();           // preconditioned residual
  GenDistrVector<Scalar> &p       = this->wksp->ret_p();           // search direction
  GenDistrVector<Scalar> &Fp      = this->wksp->ret_Fp();
  GenDistrVector<Scalar> &deltaU  = this->wksp->ret_deltaU(); 
  //GenDistrVector<Scalar> &deltaF  = this->wksp->ret_deltaF();
  GenDistrVector<Scalar> &fw      = this->wksp->ret_fw();
  GenVector<Scalar> &fc           = this->wksp->ret_fc();
  GenVector<Scalar> &uc           = this->wksp->ret_uc(); 
  GenVector<Scalar> &duc          = this->wksp->ret_duc(); 
  GenVector<Scalar> &gamma        = this->wksp->ret_gamma(); 

  int iter = 0, l = 0;
  double dd, ff, gg, lastError = 0.0, error, dualLastError = 0.0, dual_error;
  Scalar pFp, lag=0, nu;    //CRW
  GenVector<Scalar> mu(ngrbms, 0.0);
  this->fetiInfo->project_g = false;

  // initialize rho, alphabar and M
  double Fnorm = computeFNorm();
  double rho = this->fetiInfo->rho_cntl*Fnorm; // rho = rho_cntl*||F|| 
  double alphabar = this->fetiInfo->alphabar_cntl/MAX(Fnorm,rho); // alphabar = alphabar_cntl/MAX(||F||,rho) 
                                                                  // if 0 <= alphabar_cntl <= 2 then innerloop convergence is semi-monotonic for 2-step expansion
  double M = this->fetiInfo->M;

  // extract fr, fc and fw from f
  ff = extractForceVectors(f, fr, fc, fw); if(ff == 0.0) { u.zero(); return; }

  // Compute initial lagrange multipliers: lambda^0 = -G*(G^T*G)^-1 * e, where e = R^T*f
  // also for feti-dpc: expand active set and re-compute initial lagrange multipliers if lambda^0 is not feasible wrt inequality constraints
  computeL0(lambda, f);
  this->fetiInfo->cq_type = FetiInfo::nocq;
  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::unmarkRedundant);

  // Compute initial residual r^0 =  F'*lambda^0 + d'
  localSolveAndJump(fr, lambda, ur, fc, uc, r, fw, rho, mu);

  // Compute initial projected residual w^0 = P^T * r^0
  tProject(r, w, dual_error);

  // Save and print initial errors
  dd = dual_error; gg = (ngrbms) ? gamma.sqNorm() : 0.0; double ww0 = w.sqNorm();
  if(verboseFlag) filePrint(stderr," ... Initial residual norm = %e\n", sqrt(dd));

  // Save initial state (used to compute lagrangian and update penalty parameter)
  if(ngrbms) this->wksp->save_initial();

  if(verboseFlag) filePrint(stderr," Iteration  Relative Primal Error  Relative Dual Error\n");

  for(int outerIter = 0; outerIter < this->fetiInfo->maxouterit; ++outerIter, ++iter, ++iterTotal) { // outer iteration on mu

    for(int innerIter = 0; innerIter < this->fetiInfo->maxinnerit; ++innerIter, ++iter, ++iterTotal) { // inner iteration on lambdabar

      init_iter();

      // Precondition and estimate primal error: z = M^-1 * w 
      bool error_flag = ((ngrbms == 0 && this->fetiInfo->cgal_prec) || 
                         ((sqrt(gg) <= this->fetiInfo->equi_tol) && (sqrt(dual_error) < this->fetiInfo->dual_tol*sqrt(dd))));
      error = preCondition(w, z, rho, error_flag); 

      // Print error/s and test for convergence and stagnation of primal error (iff outerloop and dual error converged)
      bool print_errors = ((this->fetiInfo->numPrint() > 0) && (innerIter % this->fetiInfo->numPrint() == 0) && verboseFlag);
      if(error_flag) { // outerloop converged and inner loop converged wrt dual error, print & check primal error
        if(print_errors) filePrint(stderr," %4d %23.6e %21.6e\n", iter, sqrt(error/ff), sqrt(dual_error/dd));
        if(sqrt(error) < this->fetiInfo->tol*sqrt(ff)) break;  // converged
        if((lastError != 0.0) && (DABS(error-lastError)/lastError < this->fetiInfo->stagnation_tol)) break; // stagnated (primal)
        lastError = error;
      }
      else {
        if(print_errors) filePrint(stderr," %4d %45.6e\n", iter, sqrt(dual_error/dd));
      }

      // Adaptive precision control (iff outerloop not converged)
      if(innerIter > 0 && ngrbms && (sqrt(gg) > this->fetiInfo->equi_tol) && (sqrt(w.sqNorm()) < MIN(M*sqrt(gg), this->fetiInfo->eta*sqrt(ww0)))) break;

      // Test for stagnation of inner loop using dual error
      if((dualLastError != 0.0) && (DABS(dual_error-dualLastError)/dualLastError < this->fetiInfo->dual_stagnation_tol)) break; 
      dualLastError = dual_error;

      // Test for maximum number of iterations
      if((innerIter == this->fetiInfo->maxinnerit) || (iter == this->maxiter)) break;

      // Krylov acceleration
      if(this->fetiInfo->nlPrecFlg) nlPreCondition(w, z);

      // Re-project: y = P * z
      project(z, y); 

      // Search direction
      orthogonalize(y, p);
   
      // Step length
      pFp = localSolveAndJump(p, dur, duc, Fp, rho);
      nu = -((w*p) / pFp);

      // Update solution: lambda += nu*p, r += nu*Fp, ur += nu*dur, uc += nu*duc
      // includes optional linesearch to adjust step length nu if wolfe-armijo conditions are not satisfied
      // also for feti-dpc: reduce step and expand active set if not feasible wrt inequality constraints (dual planing)
      if(globalFlagCtc) dualStatusChange = primalStatusChange = stepLengthChange = false;
      if(this->fetiInfo->linesearch) { l = 0; init_linesearch(); }
      update(nu, lambda, p, r, Fp, ur, dur, uc, duc, l, alphabar, rho);

      // Project: w = P^T * r
      // also for feti-dpc: contract active set if not proportional (primal planing)
      tProject(r, w, dual_error);

      // add search direction to orthoset or reset if necessary
      orthoAdd(p, Fp, pFp);

      // Compute feasibility error wrt equality constraints: gamma = G^T*lambda+e
      feasible(lambda, false, false); gg = (ngrbms) ? gamma.sqNorm() : 0.0;

    }

    // Check for maximum number of iterations
    if((outerIter == this->fetiInfo->maxouterit) || (iter == this->maxiter)) 
      { if(this->myCPU==0) cerr << " ... Failed to converge within maximum specified number of iterations \n"; break; }

    if(ngrbms) {
      // Check for convergence
      if(ngrbms && this->myCPU==0) cerr << " ... Outer iteration " << outerIter+1 << ": self-equil. error = " << sqrt(gg/ee) << endl;
      if((sqrt(gg) <= this->fetiInfo->equi_tol) && (sqrt(error) < this->fetiInfo->tol*sqrt(ff))) break; // converged

      // Compute lagrangian: L = 0.5*lambda^T*(r+d)+c 
      double lag_prev = ScalarTypes::Real(lag);
      lag = lagrangian(lambda, r, rho, mu); //cerr << "lag = " << setw(24) << setprecision(16) << lag << endl;

      // Update rho or M if necessary
      double rho_prev = rho;
      if((outerIter > 0) && (this->fetiInfo->beta != 1.0) && (ScalarTypes::Real(lag) < lag_prev+0.5*rho*gg)) {
        if(this->fetiInfo->cgal_adapt_icntl == 0) { // update penalty rho
          rho = MIN(this->fetiInfo->beta*rho, this->fetiInfo->maxrho);
          if(this->myCPU == 0 && rho != rho_prev) cerr << " ... Updating penalty parameter, rho = " << rho << endl;
          alphabar = this->fetiInfo->alphabar_cntl/MAX(Fnorm,rho);
        } else { // update precision M
          double M_prev = M;
          M = MAX(M/this->fetiInfo->beta, this->fetiInfo->minM);
          if(this->myCPU == 0 && M != M_prev) cerr << " ... Updating precision, M = " << M << endl;
        }
      }

      // Update mu: mu^(k+1) = mu^k + rho^k*(G^T*lambda+e)
      mu.linAdd(rho_prev,gamma);

      // Update r due to change in rho and/or mu: r += rho*G*(G^T*lambda+e)
      multG(gamma, r, rho, 1.0,-1); 

      // Update r0 due to change in mu and rho (used for lagrangian): r0 += G*(rho_prev*gamma + (rho-rho_prev)*gamma0)
      GenDistrVector<Scalar> &r0 = this->wksp->ret_r0(); GenVector<Scalar> &gamma0 = this->wksp->ret_gamma0();
      GenVector<Scalar> v(ngrbms); v.linC(rho_prev, gamma, rho-rho_prev, gamma0);
      multG(v, r0, 1.0, 1.0,-1);

      // Reset orthoset
      this->resetOrthoSet();

      // Update w: w = P^T*r
      if(globalFlagCtc) dualStatusChange = primalStatusChange = false;
      tProject(r, w, dual_error);

      lastError = dualLastError = 0.0;
    }
    else break;  // only need 1 outer iteration when ngrbms = 0
  } // end of outer iteration loop

  // Assemble and store primal solution u
  ur += deltaU; // make solution compatible ur += deltaU
  mergeSolution(ur, uc, u, lambda);

  // Store number of iterations, errors, timings, and memory used
  this->setAndStoreInfo(iter+1, (error/ff), sqrt(dual_error/dd));
  if(this->numSystems == 1) this->times.memoryFETI += memoryUsed();
  this->times.solve += getTime(); t7 += getTime();
  this->times.iterations[this->numSystems-1].cpuTime = this->times.solve;
  printSummary(iter);
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::solveGMRES(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u)
{
 this->times.solve -= getTime(); 
 t7 = -getTime();
 if(this->oSetGMRES->numDir() > 0) this->oSetGMRES->reInit(); // just in case someone tries to use MRHS with GMRES

 GenDistrVector<Scalar> &ur = this->wksp->ret_ur(); ur.zero();
 GenDistrVector<Scalar> &dur = this->wksp->ret_du(); dur.zero();
 GenDistrVector<Scalar> &lambda0 = this->wksp->ret_lambda();
 GenDistrVector<Scalar> &r = this->wksp->ret_r(); r.zero();
 GenDistrVector<Scalar> &z = this->wksp->ret_z(); z.zero();
 GenDistrVector<Scalar> &Fp = this->wksp->ret_Fp(); Fp.zero();
 GenVector<Scalar> &uc = this->wksp->ret_uc();
 GenVector<Scalar> &duc = this->wksp->ret_duc();
 GenDistrVector<Scalar> &deltaU  = this->wksp->ret_deltaU(); deltaU.zero();

 GenDistrVector<Scalar> &fr = this->wksp->ret_fr(); 
 GenVector<Scalar> &fc = this->wksp->ret_fc();
 GenDistrVector<Scalar> &fw = this->wksp->ret_fw();

 double ff = extractForceVectors(f, fr, fc, fw);

 computeL0(lambda0, f); // PJSA 1-23-08

 // Solve uc = Kcc^-1(  fc + Krr^-1 Krc^T lambda0) 
 //          = Kcc^-1(  fc )
 // Solve ur = Krr^-1 ( fr - Br^T lambda0 - Krc uc)
 //   and  r = Br ur
 localSolveAndJump(fr, lambda0, ur, fc, uc, r, fw);

 double r0Norm2 = r.sqNorm();  
 if(verboseFlag) filePrint(stderr," ... Initial residual norm %e\n", sqrt(r0Norm2) );
 if(r0Norm2 == 0.0) {  mergeSolution(ur, uc, u, lambda0); return; }

 // Precondition r, 
 double error = preCondition(r, z);  // z = F^-1 r

// if((this->fetiInfo->numPrint() > 0) && (this->fetiInfo->numPrint() < 10))
 if(verboseFlag)
   filePrint(stderr," Iteration  Relative Primal Error  Relative Preconditioned Dual Error\n");

 GenDistrVector<Scalar> *rzero = new GenDistrVector<Scalar>(this->interface);
 GenDistrVector<Scalar> *zzero = new GenDistrVector<Scalar>(this->interface);
 GenDistrVector<Scalar> *medvec = new GenDistrVector<Scalar>(this->interface);
 GenDistrVector<Scalar> *lambda = new GenDistrVector<Scalar>(this->interface);

 *rzero = r;
 *zzero = z; 

 initGMRES(z);

 bool primalresidual = this->fetiInfo->gmresResidual;
 int J = 0;
 const int ReStep = this->fetiInfo->maxortho; // PJSA 1-23-08 (for restarted GMRES)

 for(int iter = 0; true; ++iter) {
   // Arnoldi iteration (Algorithm see Saad SISC) 
   for (int j=0; j<ReStep; j++, J++) {

     localSolveAndJump(z, dur, duc, Fp); // Fp = F*z

     error = preCondition(Fp, *medvec);   // medvec = M^-1*Fp

     // Do Arnoldi step 
     double resGMRES = orthoAddGMRES(z, *medvec);   

     if((fabs(resGMRES)<=sqrt(this->epsilon2*ff)) || (J == this->maxiter-1) || primalresidual) {

       primalresidual = true; // Since now we compute the primal residual in each step
         
       GMRESSolution(*lambda);

       localSolveAndJump(*lambda, dur, duc, Fp);  // Fp = F*lambda

       r.linC(1.0,*rzero,1.0,Fp); // r = r0 + Fp

       error = preCondition(r, *medvec); // medvec = M^-1*r

       bool isConverged = true; // XXXX ((sqrt(error) < sqrt(this->epsilon2*ff)) || (J == this->maxiter-1));

       if(verboseFlag && (isConverged || ((this->fetiInfo->numPrint() > 0) && (J % this->fetiInfo->numPrint() == 0))))
         filePrint(stderr, "%4d %23.6e %23.6e\n", iter*ReStep+j,sqrt(error/ff),fabs(resGMRES)/sqrt(ff));

       // Determine convergence
       if (isConverged) {
         // For timing and records
         this->times.iterations[this->numSystems].stagnated = 0; 
         // Store number of iterations, primal error and dual error
         this->setAndStoreInfo(iter*ReStep+j, error/ff, resGMRES*resGMRES/ff);
         if(this->numSystems == 1) this->times.memoryFETI += memoryUsed();
         this->times.solve += getTime();
         this->times.iterations[this->numSystems-1].cpuTime = this->times.solve;

         t7 += getTime();
         printSummary(iter*ReStep+j+1); // PJSA

         *lambda += lambda0; // PJSA 1-23-08
         ur += dur;
         uc += duc;
         ur += deltaU; // make solution compatible u += deltaU
         mergeSolution(ur, uc, u, *lambda);

         if(rzero) delete rzero;
         if(zzero) delete zzero;
         if(medvec) delete medvec;
         if(lambda) delete lambda; 
         return;
       }
     } 
     else 
       if(verboseFlag && ((this->fetiInfo->numPrint() > 0) && (j % this->fetiInfo->numPrint() == 0)))
         filePrint(stderr, "%4d %47.6e\n", iter*ReStep+j, fabs(resGMRES)/sqrt(ff));
   }

   // PJSA 1-23-08 restart GMRES
   if(verboseFlag) filePrint(stderr, " *** Krylov Space Full - Restarting GMRES \n");
   if(!primalresidual) {
     GMRESSolution(*lambda);  // compute incremental solution lambda
     localSolveAndJump(*lambda, dur, duc, Fp); // Fp = F*lambda
     r.linC(1.0,*rzero,1.0,Fp); // r = rzero + Fp;
   }
   *rzero       = r;
   uc += duc;
   ur += dur;

   error = preCondition(Fp, *medvec);  // medvec = M^-1*Fp
   z.linC(1.0,*zzero,1.0,*medvec); // z = zzero + medvec;
   *zzero = z;

   lambda0 += (*lambda);
   primalresidual = this->fetiInfo->gmresResidual;  // primalresidual might not be reached after restart

   this->oSetGMRES->reInit(); // Reinitialize Krylov space and set z of last step as initial vector
   initGMRES(z);
 }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::solveGCR(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u)
{
 t7 = -getTime();
 this->times.solve -= getTime();

 GenDistrVector<Scalar> &fr      = this->wksp->ret_fr();
 GenDistrVector<Scalar> &ur      = this->wksp->ret_ur();
 GenDistrVector<Scalar> &dur     = this->wksp->ret_du();
 GenDistrVector<Scalar> &lambda0 = this->wksp->ret_lambda(); // Lagrange multipliers
 GenDistrVector<Scalar> &r       = this->wksp->ret_r();       // residual
 GenDistrVector<Scalar> &w       = this->wksp->ret_w();       // projected residual
 GenDistrVector<Scalar> &y       = this->wksp->ret_y();       // re-projected preconditioned residual
 GenDistrVector<Scalar> &z       = this->wksp->ret_z();       // preconditioned residual
 GenDistrVector<Scalar> &p       = this->wksp->ret_p();       // search direction
 GenDistrVector<Scalar> &Fp      = this->wksp->ret_Fp();
 GenDistrVector<Scalar> &Fr      = this->wksp->ret_Fr();
 GenDistrVector<Scalar> &deltaU  = this->wksp->ret_deltaU();
 GenDistrVector<Scalar> &deltaF  = this->wksp->ret_deltaF();
 GenDistrVector<Scalar> &fw      = this->wksp->ret_fw();
 GenVector<Scalar> &fc           = this->wksp->ret_fc();
 GenVector<Scalar> &uc           = this->wksp->ret_uc();
 GenVector<Scalar> &duc          = this->wksp->ret_duc();

 double ff, r0Norm2, ww, lastError, error = 0.0, dummy;
 int iter;

 deltaU.zero(); deltaF.zero();

 // extract fr, fc and fw from f
 ff = extractForceVectors(f, fr, fc, fw);

 // compute the initial guess for for the Lagrange Multipliers (zero if no projection)
 computeL0(lambda0, f);

 // compute the initial residual
 // Solve uc = Kcc^-1(  fc + Krr^-1 Krc^T lambda0)
 //       ur = Krr^-1 ( fr - Br^T lambda0 - Krc uc)
 //       r = Br ur
 localSolveAndJump(fr, lambda0, ur, fc, uc, r, fw);

 // project: w = Pstar^t r
 tProject(r, w, ww);

 // save and print initial Dual error
 r0Norm2 = ww;
 if(verboseFlag) filePrint(stderr," ... Initial residual norm %e\n", sqrt(r0Norm2));
 if(r0Norm2 == 0.0) { mergeSolution(ur, uc, u, lambda0); return; }

 // multiple rhs prediction
 if((r0Norm2 > 0.0) && predictGCR(w, lambda0)) {
   localSolveAndJump(fr, lambda0, ur, fc, uc, r, fw); 
   tProject(r, w, ww);  // project: w = Pstar^t r
   if(verboseFlag)
     filePrint(stderr," ... Initial residual norm after MRHS prediction %e\n", sqrt(ww));
 }

 if(verboseFlag)
   filePrint(stderr," Iteration  Relative Primal Error  Relative Dual Error\n");

 for(iter = 0; iter < this->fetiInfo->maxit; ++iter) {
   dbg_alloca(0);
   iterTotal++;
   lastError = error;

   // Precondition: z = F^-1 w
   error = preCondition(w, z);

   // Compute and print residuals
   if((this->fetiInfo->numPrint() > 0) && (iter % this->fetiInfo->numPrint() == 0) && verboseFlag) {
     filePrint(stderr," %4d %23.6e %21.6e\n", iter, sqrt(error/ff), sqrt(ww/r0Norm2));
   }

   // Test for convergence or maximum number of iterations
   if((error < this->epsilon2 * ff) || (iter == this->maxiter-1)) {
     this->times.iterations[this->numSystems].stagnated = 0;
     if(!((this->fetiInfo->numPrint() > 0) && (iter % this->fetiInfo->numPrint() == 0) && verboseFlag)){
       if(verboseFlag) filePrint(stderr," %4d %23.6e %21.6e\n", iter, sqrt(error/ff), sqrt(ww/r0Norm2)); 
     }
     break;
   }

   // Test for stagnation
   if(DABS(error-lastError) < (this->fetiInfo->stagnation_tol * lastError)) {
     this->times.setStagnate(this->numSystems);
     filePrint(stderr, "STAGNATION: Relative Primal Error Reached = "
               "%e %e %e %e\n", sqrt(error/ff), DABS(error-lastError), error, lastError);
     break;
   }

   // Re-project: y = P * z
   project(z, y);

   localSolveAndJump(y, dur, duc, Fr);
   tProject(Fr, Fr, dummy);

   orthogonalizeGCR(y, Fr, p, Fp);  // computes new p, Fp

   Scalar FpFp = Fp * Fp;
   orthoAddGCR(p, Fp, FpFp);

   Scalar rFp = w * Fp;
   Scalar alphaCR = -(rFp/FpFp);
   doubleUpdate(alphaCR, lambda0, p, r, Fp); // lambda += alphaCR*p; r += alphaCR*Fp

   // Project:  w = Pstar^t r
   tProject(r, w, ww);
 }

 // get final solution 
 ur.zero(); uc.zero(); r.zero();
 localSolveAndJump(fr, lambda0, ur, fc, uc, r, fw);
 tProject(r, w, ww); // ???? doesn't give same projected u as CG

 // merge ur and uc into u
 mergeSolution(ur, uc, u, lambda0);

 // Store number of iterations, errors, timings and memory used
 this->setAndStoreInfo(iter+1, error/ff, sqrt(ww/r0Norm2));
 if(this->numSystems == 1) this->times.memoryFETI += memoryUsed();
 this->times.solve += getTime(); t7 += getTime();
 this->times.iterations[this->numSystems-1].cpuTime = this->times.solve;
 printSummary(iter);
}

template<class Scalar>
double
GenFetiDPSolver<Scalar>::extractForceVectors(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fr,
                                             GenVector<Scalar> &fc, GenDistrVector<Scalar> &fw)
{
  // distribute force  for sfem inpc
  GenDistrVector<Scalar> *f_copy = 0;
  if(domain->solInfo().inpc) {
     f_copy = new GenDistrVector<Scalar>(f);
     execParal1R(this->nsub, this, &GenFetiDPSolver<Scalar>::fSplit, f);
  }

  // extract fr from f
  fr.zero();
  execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::extractFr, f, fr);

  // extract fw from f
  if(domain->solInfo().isCoupled) {
    fw.zero();
    execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::extractFw, f, fw);
  }

  // Assemble and split fr and fw on subdomain interface (note: f for first system is already split by topological scaling)
  if((this->numSystems == 0 && this->fetiInfo->scaling == FetiInfo::kscaling) || (this->numSystems > 0 && this->fetiInfo->rescalef)) {
    if(domain->solInfo().isCoupled) distributeForce(fr, fw);
    else distributeForce(fr);
  }
  double ffr = fr.sqNorm();
  double ffw = (domain->solInfo().isCoupled) ? fw.sqNorm() : 0.0;

  // extract fc from f
  getFc(f, fc);
#ifdef DISTRIBUTED
  GenVector<Scalar> fc_copy(fc);
  this->fetiCom->globalSum(fc_copy.size(), fc_copy.data());
  double ffc = fc_copy.sqNorm();
#else
  double ffc = fc.sqNorm();
#endif

  double ff = f.sqNorm();

  // print norms
  if((this->fetiInfo->numPrint() > 0) && (this->fetiInfo->numPrint() < 10) && verboseFlag) {
    filePrint(stderr, " ... f*f   %e             ...\n", ff);
    filePrint(stderr, " ... fr*fr %e             ...\n", ffr);
    filePrint(stderr, " ... fc*fc %e             ...\n", ffc);
    if(domain->solInfo().isCoupled)
      filePrint(stderr, " ... fw*fw %e             ...\n", ffw);
  }
  
  if(domain->solInfo().inpc) f = (*f_copy);

  if(ff == 0.0) filePrint(stderr, " *** WARNING: norm of rhs = 0 \n");
  return (ff == 0.0) ? 1.0 : ff;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::printSummary(int iter)
{
  if(verboseFlag) {
    if(globalFlagCtc) {
      //int tmp_count2 = 0;
      filePrint(stderr," ---------------------------------------------------------------\n");
      filePrint(stderr," number of main iter.                                = %5d ...\n", iter);
      filePrint(stderr," number of iter. with expansion (dual planing)       = %5d ...\n", nStatChDual); 
      if(this->fetiInfo->outerloop != FetiInfo::CGAL)
        filePrint(stderr," number of dual planing sub-iter.                    = %5d ...\n", nSubIterDual);
      filePrint(stderr," number of iter. with proportioning (primal planing) = %5d ...\n", nStatChPrimal);
      if(this->fetiInfo->outerloop != FetiInfo::CGAL)
        filePrint(stderr," number of primal planing sub-iter.                  = %5d ...\n", nSubIterPrimal-nStatChPrimal);
      if(this->fetiInfo->linesearch)
        filePrint(stderr," number of linesearch sub-iter.                      = %5d ...\n", nLinesearchIter);
      filePrint(stderr," total number of matrix-vector (Fp) operations       = %5d ...\n", iter+1+nExtraFp); 
      filePrint(stderr," number of GtG rebuilds                              = %5d ...\n", nRebuildGtG);
      filePrint(stderr," ---------------------------------------------------------------\n");
    }
    else {
      filePrint(stderr," --------------------------------------\n");
      filePrint(stderr," number of main iter.       = %5d ...\n", iter);
      filePrint(stderr," --------------------------------------\n");
    }
    filePrint(stderr," Assembly of Kcc*   %13.4f s ...\n",t0/1000.0);
    filePrint(stderr," Krc^T Krr^-1 f     %13.4f s ...\n",t1/1000.0);
    filePrint(stderr," Assembly of fc*    %13.4f s ...\n",t2/1000.0);
    filePrint(stderr," Kcc*^-1 fc*        %13.4f s ...\n",this->times.project/1000.0);
    filePrint(stderr," Local Solve        %13.4f s ...\n",t4/1000.0);
    filePrint(stderr," --------------------------------------\n");
    filePrint(stderr," Kcc - Krc^t Krr Krc %12.4f s ...\n",t5/1000.0);
    filePrint(stderr," Total Making Kcc    %12.4f s ...\n",this->times.coarse1/1000.0);
    filePrint(stderr," --------------------------------------\n");
    filePrint(stderr," Total Making FETI  %13.4f s ...\n",t6/1000.0);
    filePrint(stderr," Total Solve RHS#   %13.4f s ...\n",t7/1000.0);  // PJSA 2-10-05
  }
  //if(this->myCPU == 0) cerr << "iterTotal = " << iterTotal << endl;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::orthoAdd(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar pFp)
{
  if(dualStatusChange || primalStatusChange || stepLengthChange) this->resetOrthoSet(); 
  else GenFetiSolver<Scalar>::orthoAdd(p, Fp, pFp);
}

template<class Scalar>
double
GenFetiDPSolver<Scalar>::preCondition(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Mv, double rho, bool errorFlag)
{
  // top-level precondition function for CG with augmented lagrangian
  // compute Mv = (P^T*M^-1*P + 1/rho*G*G^T)*v
  if(this->fetiInfo->cgal_prec) {
    if(ngrbms == 0) return preCondition(v, Mv, errorFlag);
    else {
      project(v, Mv, 2);
      preCondition(Mv, Mv, false);
      project(Mv, Mv, 2);
      if(rho != 0.0) {
        GenVector<Scalar> beta(ngrbms);
        trMultG(v, beta, 1.0, 0.0,-1); // beta = G^T*v
        double r = 1.0/rho*this->fetiInfo->cgal_prec_cntl;
        multG(beta, Mv, r, 1.0,-1);  // Mv += 1/rho*G*beta
      }
    }
  }
  else Mv = v;
  return (errorFlag) ? primalError() : numeric_limits<double>::max();
}

template<class Scalar>
double
GenFetiDPSolver<Scalar>::preCondition(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Mv, bool errorFlag)
{
  // this function does dual mpc (CCt) preconditioning in addition to the usual feti preconditioning
  double error;
  if(mpcPrecon) {
    if((dualStatusChange || primalStatusChange) && this->fetiInfo->rebuildcct && CCtsolver) rebuildCCt();
    if(&Mv != &v) Mv = v; cctSolveMpc(Mv); // Mv = CCt^-1*v
    error = GenFetiSolver<Scalar>::preCondition(Mv, Mv, errorFlag); 
    cctSolveMpc(Mv); 
  }
  else error = GenFetiSolver<Scalar>::preCondition(v, Mv, errorFlag);
  return (errorFlag) ? error : numeric_limits<double>::max();
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::mergeSolution(GenDistrVector<Scalar> &ur, GenVector<Scalar> &uc,
                                       GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &lambda)
{
  execParal4R(this->nsub, this, &GenFetiDPSolver<Scalar>::mergeUr, ur, uc, u, lambda);
  if(ngrbms) { // add rigid body modes
    GenVector<Scalar> &alpha = this->wksp->ret_alpha(); // amplitudes of the rigid body modes
    execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::addRalpha, u, alpha);
    if(this->fetiInfo->uproj) computeProjectedDisplacement(u);
    else filePrint(stderr, " ... Do not project the displacement ...\n"); //HB
  }
  feasible(lambda); // check feasibility of all dual constraints
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::mergeUr(int iSub, GenDistrVector<Scalar> &ur, GenVector<Scalar> &uc, 
                                 GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &lambda) 
{
  this->sd[iSub]->mergeUr(ur.subData(iSub), uc.data(), u.subData(iSub), lambda.subData(iSub));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::localSolveAndJump(GenDistrVector<Scalar> &fr, GenDistrVector<Scalar> &lambda,
                                           GenDistrVector<Scalar> &ur, GenVector<Scalar> &fc,
                                           GenVector<Scalar> &uc, GenDistrVector<Scalar> &r,
                                           GenDistrVector<Scalar> &fw)
{
 startTimerMemory(this->times.sAndJ, this->times.memorySAndJ);
 r.zero();

 GenVector<Scalar> &FcStar(uc); 

 if(KccSolver) { 
   // Step 1: fc^*(s) = fc^(s) - (Krc^T Krr^-1)^(s) (fr^(s) - Br^(s)T lambda) - Bc^(s)tilde^T lambda
   t1 -= getTime();
   execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::makeFc, fr, lambda);
   FcStar = fc;
   t1 += getTime();

   // Step 2: Assemble local fc^* into global FcStar
   t2 -= getTime();
   assembleFcStar(FcStar);
   t2 += getTime();
#ifdef DISTRIBUTED
   this->fetiCom->globalSum(FcStar.size(), FcStar.data());
#endif

   // Step 3: Solve uc = Kcc^-1 FcStar
   this->times.project -= getTime();
   if(this->glNumMpc_primal > 0) execParal(this->mpcToSub_primal->csize(), this, &GenFetiDPSolver<Scalar>::addMpcRHS, FcStar.data()); 
   KccSolver->reSolve(FcStar);  // now Fcstar is uc;
   this->times.project += getTime();
 }
 else FcStar.zero();

 // Step 3.5: fr^(s) = fr^(s) - Krc^(s) uc^(s)
 t4 -= getTime();
 GenDistrVector<Scalar> &fr2 = this->wksp->ret_fr2(); fr2 = fr;
 execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::multKrc, fr2, FcStar);

 // Step 4:
 if(domain->solInfo().isCoupled && !this->fetiInfo->fsi_element) {
   timedParal5R(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolveCoupled1,
                ur, r, fr2, lambda, FcStar);
// JLchange 
//   this->wiPat->exchange();
   if (this->fetiInfo->fsi_corner == 0) this->wiPat->exchange();
   
   timedParal6R(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolveCoupled2,
                ur, r, fr2, lambda, FcStar, fw);
 }
 else 
   timedParal5R(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolve, 
                ur, r, fr2, lambda, FcStar); // FcStar is necessary for contact

 t4 += getTime();
 this->vPat->exchange();
 timedParal1R(this->times.solveAndJump, this->nsub, this, &GenFetiSolver<Scalar>::interfaceDiff, r);

 // Step 5:
 if(this->glNumMpc > 0) execParal1R(this->nsub, this, &GenFetiDPSolver<Scalar>::subtractMpcRhs, r);

 stopTimerMemory(this->times.sAndJ, this->times.memorySAndJ);
}

template<class Scalar> 
Scalar
GenFetiDPSolver<Scalar>::localSolveAndJump(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &dur, 
                                           GenVector<Scalar> &duc, GenDistrVector<Scalar> &Fp)
{
 startTimerMemory(this->times.sAndJ, this->times.memorySAndJ);

 GenVector<Scalar> &FcStar(duc);
 FcStar.zero();
 if(KccSolver) {
   // Step 1: fc^(s) = - (Krc^T Krr^-1)^(s) ( Br^(s)T * p) + Bc^(s)T * p
   t1 -= getTime();
   execParal1R(this->nsub, this, &GenFetiDPSolver<Scalar>::makeFcB, p);
   t1 += getTime();

   // Step 2: Assemble fc^(s) into FcStar
   t2 -= getTime();
   assembleFcStar(FcStar);
   t2 += getTime();
#ifdef DISTRIBUTED
   this->fetiCom->globalSum(FcStar.size(), FcStar.data());
#endif

   // Step 3: Solve uc = Kcc^-1 FcStar
   this->times.project -= getTime();
   KccSolver->reSolve(FcStar);
   this->times.project += getTime();
 }

 // Step 3.5: fr2^(s) = - Krc^(s) uc^(s)
 t4 -= getTime();
 GenDistrVector<Scalar> &fr2 = this->wksp->ret_fr2(); fr2.zero();
 execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::multKrc, fr2, FcStar);

 // Step 4:
 // CKT : fr=Krc (-uc) 
 // CKT : Fp = B Krr^-1 [B^T p + Krc (-uc)] - Bc (-uc)
 // B contains both usual Br u = +- u and ctc restriction B u = +- u.n
 // Bc uc < > 0 if some corners nodes are in ctc. Then Bc uc = +- uc.n
 if(domain->solInfo().isCoupled && !this->fetiInfo->fsi_element) {
   timedParal5R(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolveCoupled1,
                dur, Fp, fr2, p, FcStar);
// JLchange 
//   this->wiPat->exchange();
   if (this->fetiInfo->fsi_corner == 0) this->wiPat->exchange();
   timedParal5R(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolveCoupled2,
                dur, Fp, fr2, p, FcStar);
 }
 else
   timedParal5R(this->times.solveAndJump, this->nsub, this, &GenFetiDPSolver<Scalar>::subdomainSolve,
                dur, Fp, fr2, p, FcStar);
 t4 += getTime();
 this->vPat->exchange();
 timedParal1R(this->times.solveAndJump, this->nsub, this, &GenFetiSolver<Scalar>::interfaceDiff, Fp);
 Scalar ret = p*Fp;
#ifdef DEBUG_FETI
 Scalar pHFp = ScalarTypes::conj(ret); // note: p*Fp = (Fp)^H p therefore p^H Fp = conj(p*Fp)
 if(this->myCPU == 0 && (this->fetiInfo->outerloop == FetiInfo::CG || this->fetiInfo->outerloop == FetiInfo::CGAL)
    && (ScalarTypes::Real(pHFp) < 0.0 || fabs(ScalarTypes::Imag(pHFp)) > 1.0e-10))
   cerr << " *** WARNING: x^H F x = " << pHFp << ", must be positive and real for any x when F is Hermitian and positive definite. CG may not work \n"; // XXXX
#endif
 stopTimerMemory(this->times.sAndJ, this->times.memorySAndJ);
 return ret;
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::assembleFcStar(GenVector<Scalar> &FcStar)
{
 int i, iSub;
 for(iSub = 0; iSub < this->nsub; ++iSub) {
   int numCornerDofs = this->sd[iSub]->numCoarseDofs();
   int *dofs = this->sd[iSub]->cornerEqNums;
   Scalar *fc = this->sd[iSub]->getfc();  // returns sub fcstar, not condensed
   for(i = 0; i < numCornerDofs; ++i) {  // assemble global condensed fcstar
     if(dofs[i] != -1) {
       FcStar[dofs[i]] += fc[i];
     }
   }
 }
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::multKrc(int iSub, GenDistrVector<Scalar> &fr, GenVector<Scalar> &uc)
{
  this->sd[iSub]->multKrc(fr.subData(iSub), uc.data());
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::makeFc(int iSub, GenDistrVector<Scalar> &fr, /*GenVector<Scalar> &fc,*/ 
                                GenDistrVector<Scalar> &lambda)
{
  this->sd[iSub]->multfc(fr.subData(this->sd[iSub]->localSubNum()), /*fc.data(),*/
                         lambda.subData(this->sd[iSub]->localSubNum()));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::makeFcB(int iSub, GenDistrVector<Scalar> &p)
{
  this->sd[iSub]->multFcB(p.subData(this->sd[iSub]->localSubNum()));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::makeEdgeConnectivity()
{ 
  int iSub;
  int glNumSub = this->subToSub->csize();

  // First count number of edges per subdomain
  int *cx = new int[glNumSub+1];
  int *cxx = new int[glNumSub+1];
  for(iSub=0; iSub<glNumSub+1; ++iSub) { cx[iSub] = 0; cxx[iSub] = 0; }

  // count edges in parallel
  execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::countEdges, cx);

#ifdef DISTRIBUTED
  this->fetiCom->globalSum(glNumSub+1, cx);
#endif

  // Sum each subdomain's edge count to compute total number of edges
  // We modify 'edges' at the same time so that it has the first edge number
  // for each subdomain
  int numEdges      = 0;
  int totalNumEdges = 0;
#ifdef DISTRIBUTED
  int i;
  int *numEdgesPerSub = new int[glNumSub];
  for(iSub=0; iSub<glNumSub; ++iSub)
    numEdgesPerSub[iSub] = 0;

  for(iSub=0; iSub<this->nsub; ++iSub)  
    numEdgesPerSub[this->sd[iSub]->subNum()] = this->sd[iSub]->numEdgeNeighbors();  // PJSA: don't include virtual neighbors

  this->fetiCom->globalSum(glNumSub, numEdgesPerSub);

  for(iSub=0; iSub<glNumSub; ++iSub) {
    int tmp = numEdges;
    numEdges += cx[iSub];
    cx[iSub] = tmp;

    cxx[iSub] = totalNumEdges;
    totalNumEdges += numEdgesPerSub[iSub];
  }
  // make sure I can delete this
  delete [] numEdgesPerSub;
#else
  for(iSub=0; iSub<glNumSub; ++iSub) {
    int tmp = numEdges;
    numEdges += cx[iSub];
    cx[iSub] = tmp;

    cxx[iSub] = totalNumEdges;
    totalNumEdges += this->sd[iSub]->numEdgeNeighbors(); // PJSA: don't include virtual neighbors
  }
#endif
  // Find total number of edges including duplicated ones.
  // totalNumEdges = 2 * numEdges (i.e. totalNumEdges counts the duplicate
  // edges since 2 subdomains see the same edge)
  // filePrint(stderr,"nE %d tnE %d\n",numEdges, totalNumEdges);

  int *connect = new int[totalNumEdges];
  cxx[glNumSub] = totalNumEdges;

#ifdef DISTRIBUTED
  for(i=0; i<totalNumEdges; ++i) connect[i] = -1;
#endif

  execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::numberEdges, cx, cxx, connect, this->sPat);
  delete [] cx;

 // Due to the -1, how could I do a globalSum?
#ifdef DISTRIBUTED
  for(i=0; i<totalNumEdges; ++i) connect[i] += 1;
  this->fetiCom->globalSum(totalNumEdges, connect);
  for(i=0; i<totalNumEdges; ++i) connect[i] -= 1;
#endif

  this->sPat->exchange();
  execParal(this->nsub, this, &GenFetiDPSolver<Scalar>::receiveNeighbEdgeNums, cxx, connect, this->sPat);

#ifdef DISTRIBUTED
  int *connect2 = new int[totalNumEdges];
  for(i=0; i<totalNumEdges; ++i) {
    if(connect[i] == -1) 
      connect2[i] = 0;
    else
      connect2[i] = connect[i]; 
  }
  this->fetiCom->globalSum(totalNumEdges, connect2);
  for(i=0; i<totalNumEdges; ++i) {
    if(connect[i] == -1) 
     connect[i] = connect2[i];
  }
  delete [] connect2;
#endif
  this->subToEdge = new Connectivity(glNumSub, cxx, connect);
  // create the edge to subdomain connectivity
  this->edgeToSub = this->subToEdge->reverse();
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::countEdges(int iSub, int *edges)
{
  // get my global subdomain number
  int myNum = this->sd[iSub]->subNum();

  // count a subdomain's edges
  edges[myNum] = 0;
  int j;
  int numNeighbor = this->sd[iSub]->numNeighbors();
  for(j=0; j<numNeighbor; ++j) {
    int subI = this->sd[iSub]->getSComm()->subNums[j];
    if(this->sd[iSub]->isEdgeNeighbor(j) && (myNum < subI))  // PJSA
      edges[myNum] += 1;
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::numberEdges(int iSub, int *eP, int *ep2, int* edges, FSCommPattern<int> *sPat)
{
  // first compute each subdomain's starting edge number
  int myNum = this->sd[iSub]->subNum();
  int startI = eP[myNum];
  int fP = ep2[myNum];

  // number all subdomain's edges
  int jSub;
  int numNeighbor = this->sd[iSub]->numNeighbors();
  for(jSub=0; jSub<numNeighbor; ++jSub) {
    int subJ = this->sd[iSub]->getSComm()->subNums[jSub];
    FSSubRecInfo<int> sInfo = this->sPat->getSendBuffer(myNum, subJ);
    if(this->sd[iSub]->isEdgeNeighbor(jSub)) { // PJSA
      if(myNum < subJ) {
        edges[fP] = startI;
        startI++;
      }
      else edges[fP] = -1;
      // Send the numbered edges to all neighbors so that they can complete
      // the edge vector
      sInfo.data[0] = edges[fP];
      fP++;
    }
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::receiveNeighbEdgeNums(int iSub, int *eP, int* edges, FSCommPattern<int> *sPat)
{
  int myNum = this->sd[iSub]->subNum();
  // number all subdomain's edges
  int jSub;
  int numNeighbor = this->sd[iSub]->numNeighbors();
  int jEdgeN = 0;
  for(jSub=0; jSub<numNeighbor; ++jSub) {
    if(this->sd[iSub]->isEdgeNeighbor(jSub)) { // PJSA
      FSSubRecInfo<int> rInfo = this->sPat->recData(this->sd[iSub]->getSComm()->subNums[jSub], myNum);
      int en =  rInfo.data[0];
      if(en >= 0) edges[eP[myNum]+jEdgeN] = en;
      jEdgeN++;
    }
  }
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::factorLocalMatrices(int iSub)
{
  GenSolver<Scalar> *K = this->sd[iSub]->Krr;
  //int numCornerDofs = this->sd[iSub]->numCornerDofs();
  if(K) {
    K->setPrintNumGrbm(false);
    K->factor(); 
    if(K->numRBM() != 0) {
      filePrint(stderr," ... Subdomain %3d found %3d ZEMs   ...\n",
                this->sd[iSub]->localSubNum()+1,K->numRBM());
    }
  }
  this->sd[iSub]->factorKii();
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::subdomainSolve(int iSub, GenDistrVector<Scalar> &v1, 
                                        GenDistrVector<Scalar> &v2, GenDistrVector<Scalar> &v3,
                                        GenDistrVector<Scalar> &v4, GenVector<Scalar> &v5)
{
  int sn = this->sd[iSub]->localSubNum();

  Scalar *localvec  = v1.subData(sn);
  Scalar *interfvec = v2.subData(sn);
  Scalar *localsrc  = v3.subData(sn);
  Scalar *interfsrc = v4.subData(sn);
  Scalar *uc        = v5.data(); // Necessary for contact

  int localLen = this->sd[iSub]->localRLen();
  int i;

  for(i = 0; i < localLen; ++i) 
     localvec[i] = localsrc[i];

  int interfaceLen = this->sd[iSub]->interfLen();
  for(i = 0; i < interfaceLen; ++i) 
     interfvec[i] = interfsrc[i];

  this->sd[iSub]->fetiBaseOp(uc, this->sd[iSub]->Krr, localvec, interfvec);
  
  this->sd[iSub]->sendInterf(interfvec, this->vPat);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subdomainSolveCoupled1(int iSub, GenDistrVector<Scalar> &v1,
                                                GenDistrVector<Scalar> &v2, GenDistrVector<Scalar> &v3,
                                                GenDistrVector<Scalar> &v4, GenVector<Scalar> &v5)
{
  int sn = this->sd[iSub]->localSubNum();

  Scalar *localvec  = v1.subData(sn);
  Scalar *interfvec = v2.subData(sn);
  Scalar *localsrc  = v3.subData(sn);
  Scalar *interfsrc = v4.subData(sn);
  //Scalar *uc        = v5.data();

  int localLen = this->sd[iSub]->localRLen();
  int i;

  for(i = 0; i < localLen; ++i)
     localvec[i] = localsrc[i];

  int interfaceLen = this->sd[iSub]->interfLen();
  for(i = 0; i < interfaceLen; ++i)
     interfvec[i] = interfsrc[i];

  this->sd[iSub]->fetiBaseOpCoupled1(this->sd[iSub]->Krr, localvec, interfvec, this->wiPat);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subdomainSolveCoupled2(int iSub, GenDistrVector<Scalar> &v1,
                                                GenDistrVector<Scalar> &v2, GenDistrVector<Scalar> &v3,
                                                GenDistrVector<Scalar> &v4, GenVector<Scalar> &v5,
                                                GenDistrVector<Scalar> &fw)
{
  int sn = this->sd[iSub]->localSubNum();

  Scalar *localvec  = v1.subData(sn);
  Scalar *interfvec = v2.subData(sn);
  //Scalar *localsrc  = v3.subData(sn);
  //Scalar *interfsrc = v4.subData(sn);
  Scalar *uc        = v5.data();
  Scalar *localfw   = fw.subData(sn);

  this->sd[iSub]->fetiBaseOpCoupled2(uc, localvec, interfvec, this->wiPat, localfw);
  this->sd[iSub]->sendInterf(interfvec, this->vPat);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subdomainSolveCoupled2(int iSub, GenDistrVector<Scalar> &v1,
                                                GenDistrVector<Scalar> &v2, GenDistrVector<Scalar> &v3,
                                                GenDistrVector<Scalar> &v4, GenVector<Scalar> &v5)
{
  int sn = this->sd[iSub]->localSubNum();

  Scalar *localvec  = v1.subData(sn);
  Scalar *interfvec = v2.subData(sn);
  //Scalar *localsrc  = v3.subData(sn);
  //Scalar *interfsrc = v4.subData(sn);
  Scalar *uc        = v5.data();
  Scalar *localfw   = 0;

  this->sd[iSub]->fetiBaseOpCoupled2(uc, localvec, interfvec, this->wiPat, localfw);
  this->sd[iSub]->sendInterf(interfvec, this->vPat);
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::clean_up()
{
}

template<class Scalar> 
double
GenFetiDPSolver<Scalar>::getFNormSq(GenDistrVector<Scalar> &f)
{
  // this is used by nonlinear analysis, need to assemble force residual on subdomain interface for error norm
  GenDistrVector<Scalar> &fr = this->wksp->ret_fr();
  fr.zero();
  execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::extractFr, f, fr);
  distributeForce(fr);
  GenVector<Scalar> &fc  = this->wksp->ret_fc();
  getFc(f, fc);
#ifdef DISTRIBUTED
  this->fetiCom->globalSum(fc.size(), fc.data());
#endif
  return (fr.sqNorm() + fc.sqNorm());
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::subtractMpcRhs(int iSub, GenDistrVector<Scalar> &dv1) 
{
  this->sd[iSub]->subtractMpcRhs(dv1.subData(this->sd[iSub]->localSubNum()));
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::subtractKGap(int iSub, GenDistrVector<Scalar> &referenceRHS)
{
  // GR: compute f - K * contactGap
  Scalar *refRHS = referenceRHS.subData(this->sd[iSub]->localSubNum());
  this->sd[iSub]->subtractKGap(refRHS);
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::singularValueDecomposition(FullM &A, FullM &U, int ncol, int nrow, int &rank, double tol, FullM *V)
{
  int info = 0;
  int mindim = myMin(nrow,ncol);
  int maxdim = myMax(nrow,ncol);
  double max_value = A.maxAbs();
#ifdef FILERING
  for(int i=0; i<A.numCol()*A.numRow(); i++) //HB
    if(fabs((A.data())[i])<1.E-16*max_values+1.E-20) (A.data())[i] = 0.0;
#endif
  double *w    = (double*) dbg_alloca(sizeof(double)* maxdim);
  double *e    = (double*) dbg_alloca(sizeof(double)* maxdim);
  double *work = (double*) dbg_alloca(sizeof(double)* maxdim);

  for(int i=0; i<maxdim; i++) { w[i] = 0.0; e[i] = 0.0; work[i] = 0.0; }
                                                                                                  
  if(V == 0)
    Tsvdc(A.data(), ncol, ncol, nrow, w, e, U.data(), ncol, U.data(), ncol, work, 10, info);
  else 
    Tsvdc(A.data(), ncol, ncol, nrow, w, e, U.data(), ncol, V->data(), nrow, work, 11, info);

  //  ... CHECK RETURN STATUS OF DSVDC:
  if(info <  0)
    filePrint(stderr," *** ERROR: Illegal value in argument #%d in FetiDPSolver::singularValueDecomposition()\n",info);
  if(info >  0)
    filePrint(stderr," *** WARNING: %d diagonals did not converge in FetiDPSolver::singularValueDecomposition(). Check your result.\n",info);

  //  ... DETERMINE RANK
  rank = 0;
  double tolerance = max_value*tol;
  for(int i=0; i<mindim; ++i) {
    if(fabs(w[i]) > tolerance) rank += 1;
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::rebuildGtGtilda()
{
  nRebuildGtG++; nStatChDual_gtg = nStatChDual; nStatChPrimal_gtg = nStatChPrimal;
  startTimerMemory(this->times.coarse1, this->times.memoryGtG);

  if(GtGtilda == GtG) {
    GtGtilda = newSolver(this->fetiInfo->auxCoarseSolver, coarseConnectGtG, eqNumsGtG, this->fetiInfo->grbm_tol, GtGsparse);
    GtGtilda->setPrintNumTrbm(false);
  } else
  GtGtilda->zeroAll();
  execParal(nGroups1, this, &GenFetiDPSolver<Scalar>::assembleGtG, 1);
#ifdef DISTRIBUTED
  GtGtilda->unify(this->fetiCom);
#endif
  startTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
  GtGtilda->parallelFactor();
  stopTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
  if(GtGtilda->numRBM() > GtG->numRBM()) 
    filePrint(stderr, " ... GtGtilda has %d singularities for tol %e ...\n", GtGtilda->numRBM(), this->fetiInfo->grbm_tol);

  stopTimerMemory(this->times.coarse1, this->times.memoryGtG);
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::assembleGtG(int iGroup, int flag)
{
  // flag = 0: assemble G^T*G, 1: Gtilda^T*Gtilda
  // assembles groups in parallel, subdomains with same group sequentially 
  // threadsafe implementation - avoids simultaneous writing to same memory
  // note: the distributed version will work for shared memory too, but the 
  // alternative code is a bit more efficient
  int i;
#ifdef DISTRIBUTED
  for(i = 0; i < this->nsub; ++i) {
    if(this->sd[i]->group == groups[iGroup]) 
      this->sd[i]->assembleGtGsolver(GtGsparse, flag);
  }
#else
  int *grsubs = (*groupToSub)[iGroup];
  for(i = 0; i < groupToSub->num(iGroup); ++i) {
    int iSub = grsubs[i];
    this->sd[iSub]->assembleGtGsolver(GtGsparse, flag);
  }
#endif
}

template<class Scalar>
void 
GenFetiDPSolver<Scalar>::buildCCt()
{
  // PJSA: build CC^t for preconditioning mpc residual
  startTimerMemory(this->times.buildCCt, this->times.memoryBuildCCt);
  // find subs with mpcs
  mpcSubMap = new int[this->nsub];
  subsWithMpcs = new ResizeArray<GenSubDomain<Scalar> *>(0);
  numSubsWithMpcs = 0;
  for(int i=0; i<this->nsub; ++i) {
    if(this->sd[i]->numMPC > 0) {
      mpcSubMap[i] = numSubsWithMpcs;
      (*subsWithMpcs)[numSubsWithMpcs++] = this->sd[i];
    }
    else mpcSubMap[i] = -1;
  }

  switch(this->fetiInfo->mpc_precno) {
    case (FetiInfo::globalCCt) :
      CCtsolver = new GlobalCCtSolver<Scalar>(mpcToMpc, mpcToCpu, numSubsWithMpcs, subsWithMpcs->data(), 
                                              this->fetiInfo, this->fetiCom);
      break;
    case (FetiInfo::blockDiagCCt) : {
      Connectivity *blockToMpc = getBlockToMpc();
      CCtsolver = new BlockCCtSolver<Scalar>(blockToMpc, mpcToMpc, this->mpcToSub, mpcToCpu, this->numSubsWithMpcs, subsWithMpcs->data(),
                                             mpcSubMap, this->fetiInfo, this->fetiCom);
      } break;
    case (FetiInfo::subBlockDiagCCt) :
      CCtsolver = new SubBlockCCtSolver<Scalar>(mpcToMpc, this->mpcToSub, numSubsWithMpcs, subsWithMpcs->data(), 
                                                this->fetiCom, this->cpuToSub);
      break;
    case (FetiInfo::superBlockDiagCCt) : {
      Connectivity *blockToMpc = getBlockToMpc();
      bool super_flag = (this->fetiInfo->mpc_block == FetiInfo::subBlock) ? false : true;
      bool sub_flag = (this->fetiInfo->mpc_block == FetiInfo::mortarBlock) ? true : false;
      CCtsolver = new SuperBlockCCtSolver<Scalar>(blockToMpc, mpcToMpc, this->mpcToSub, mpcToCpu, numSubsWithMpcs, subsWithMpcs->data(),
                                                  this->fetiInfo, this->fetiCom, super_flag, sub_flag);
    } break;
    default :
      cerr << " *** ERROR: don't know mpc_precno = " << this->fetiInfo->mpc_precno << endl;
      break;
  }
  CCtsolver->assemble();
  CCtsolver->factor();
  stopTimerMemory(this->times.buildCCt, this->times.memoryBuildCCt);
}

template<class Scalar>
Connectivity *
GenFetiDPSolver<Scalar>::getBlockToMpc()
{
 // return a decomposition of global mpcs into blocks for block solver
 Connectivity *blockToMpc;
 switch(this->fetiInfo->mpc_block) {
   case (FetiInfo::topoBlock) :
     {
       compStruct renumber = mpcToMpc->renumByComponent(-1);  // -1 = no renumbering (optimized implementation)
       int nMpcBlocks = renumber.numComp;
       renumber.order = new int[this->glNumMpc];
       for(int i=0; i<this->glNumMpc; ++i) renumber.order[renumber.renum[i]] = i;
       int *target = new int[this->glNumMpc];
       for(int i=0; i<this->glNumMpc; ++i) target[renumber.renum[i]] = i;
       int *pointer = new int[nMpcBlocks+1];
       for(int i=0; i<=nMpcBlocks; ++i) pointer[i] = renumber.xcomp[i];
       renumber.clearMemory();
       blockToMpc = new Connectivity(nMpcBlocks, pointer, target);
     } 
     break;
   case (FetiInfo::subBlock) :
     { 
       //HB: WARNING: HARD CODED FOR SHARED MEMORY !!!
       Connectivity *subToMpc = this->mpcToSub->reverse();
       int nMpcBlocks = 0; 
       int ntarget = 0;
       for(int s=0;s<subToMpc->csize();s++){
         if(subToMpc->num(s)) { 
           nMpcBlocks++; 
           ntarget += subToMpc->num(s);
         }
       }
       int *ptr = new int[nMpcBlocks+1]; ptr[0] = 0;
       int *target = new int[ntarget];
       int iblk = 0;
       for(int s=0;s<subToMpc->csize();s++){
         if(subToMpc->num(s)) { 
           ptr[iblk+1] = ptr[iblk];
           for(int i=0;i<subToMpc->num(s);i++){
             target[ptr[iblk]+i] = (*subToMpc)[s][i];
             ptr[iblk+1]++;
           } 
           iblk++;
         }
       }
       blockToMpc = new Connectivity(nMpcBlocks,ptr,target);
       delete subToMpc;
     } 
     break;
   case (FetiInfo::mortarBlock) :
     {
       // make initial blockToMpc connectivity: make one block for each set of MPCs representing 
       // a mortar face, and put all remaining (non-mortar) MPCs together in one block 
       int numMortarMpcs = domain->GetnMortarLMPCs();
       int numOtherMpcs = this->glNumMpc - numMortarMpcs;
       if(numOtherMpcs > 0) {
         int *pointer = new int[2]; pointer[0] = 0; pointer[1] = numOtherMpcs;
         int *target = new int[numOtherMpcs];
         for(int i=0; i<numOtherMpcs; ++i) target[i] = i;
         Connectivity *otherToMpc = new Connectivity(1, pointer, target);
         if(numMortarMpcs == 0) blockToMpc = otherToMpc;
         else blockToMpc = otherToMpc->merge(domain->GetMortarToMPC());
       }
       else blockToMpc = domain->GetMortarToMPC()->copy();
     } 
     break;
 }

 // block overlapping if requested
 if(blockToMpc->csize() > 1) {
   int blockOverlapLevel = this->fetiInfo->mpcBlkOverlap;
   filePrint(stderr, " ... Block Overlap Level %d ... \n", blockOverlapLevel);
   for(int i = 0; i < blockOverlapLevel; ++i) {
     Connectivity *tmpBlockToMpc = blockToMpc->transcon(mpcToMpc);
     delete blockToMpc; blockToMpc = tmpBlockToMpc;
   }
 }

 return blockToMpc;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::rebuildCCt()
{
  CCtsolver->zeroAll();
  CCtsolver->assemble();
  CCtsolver->factor();
  nRebuildCCt++;
}

template<class Scalar>
void 
GenFetiDPSolver<Scalar>::cctSolveMpc(GenDistrVector<Scalar> &v)
{
  // PJSA: this function is used when applying the generalized preconditioner
  // to the mpc part of the residual
  // computes:  v = (CC^t)^-1 v
  startTimerMemory(this->times.precond, this->times.memoryPrecond);
  startTimerMemory(this->times.solveCCt, this->times.memorySolveCCt);
  
  CCtsolver->reSolve(v); 

  stopTimerMemory(this->times.solveCCt, this->times.memorySolveCCt);
  if(this->times.memorySolveCCt < 0) this->times.memorySolveCCt = 0;
  stopTimerMemory(this->times.precond, this->times.memoryPrecond); 
}

template<class Scalar>
GenFetiDPSolver<Scalar>::~GenFetiDPSolver() 
{
  //if(this->sPat) { delete this->sPat; this->sPat = 0; }  
  //if(this->wiPat) { delete this->wiPat; this->wiPat = 0; }

  if(ngrbmGr) { delete [] ngrbmGr; ngrbmGr = 0; }
  if(groups) { delete [] groups; groups = 0; }
  if(groupToSub != bodyToSub) { delete groupToSub; groupToSub = 0; }
  if(subToGroup && (subToGroup != subToBody)) { delete subToGroup; subToGroup = 0; }
  if(subToBody) { delete subToBody; subToBody = 0; }
  if(GtGtilda && GtGtilda != GtG) { delete GtGtilda; GtGtilda = 0; }
  if(GtG) { delete GtG; GtG = 0; }
  if(glCrnGroup) { delete [] glCrnGroup; glCrnGroup = 0; }
  
  if(KccSparse) { delete KccSparse; KccSolver = 0; KccSparse = 0; }
  if(cornerToSub) { delete cornerToSub; cornerToSub = 0; }
  if(cornerEqs) { delete cornerEqs; cornerEqs = 0; }
  if(c_cornerEqs) { delete c_cornerEqs; c_cornerEqs = 0; }
  if(coarseConnectGtG) { delete coarseConnectGtG; coarseConnectGtG = 0; }
  if(eqNumsGtG) { delete eqNumsGtG; eqNumsGtG = 0; }
  
  if(CCtsolver) { delete CCtsolver; CCtsolver = 0; }
  if(subsWithMpcs)  { delete subsWithMpcs ; subsWithMpcs = 0; }
  if(mpcSubMap) { delete [] mpcSubMap; mpcSubMap = 0; }
  // don't delete: this->sd, this->subToSub, this->mpcToSub, mpcToMpc, this->cpuToSub, this->fetiInfo, bodyToSub 
  // interfPattern, this->glSubToLoc, mpcToCpu 
  if(lambda_total) { delete lambda_total; lambda_total = 0; }
  if(mpcPat) { delete mpcPat; mpcPat = 0; }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::getLocalMpcForces(int iSub, double *mpcLambda)
{
  // mpcLambda is the local MPC forces for subdomain iSub (required for salinas)
  GenVector<Scalar> &uc = this->wksp->ret_uc();
  this->sd[iSub]->getLocalMpcForces(mpcLambda, c_cornerEqs, mpcOffset, uc);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::addMpcRHS(int iMPC, Scalar *fcstar)
{
 int dof = c_cornerEqs->firstdof(mpcOffset+iMPC);
 int myNum = this->glSubToLoc[(*this->mpcToSub_primal)[iMPC][0]];
 if(myNum >= 0)
   fcstar[dof] += this->sd[myNum]->getMpcRhs_primal(iMPC);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::setBodyRBMoffset(int iSub, int *zColOffset)
{
  int subBody = (*subToBody)[this->sd[iSub]->subNum()][0];
  this->sd[iSub]->setBodyRBMoffset(zColOffset[subBody]);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::initialize()
{
  KccSolver = 0; KccSparse = 0; glNumCorners = 0;
  cornerToSub = 0; cornerEqs = 0; c_cornerEqs = 0;
  mpcOffset = 0; ngrbm = 0;
  ngrbmGr = 0; nGroups = 0; nGroups1 = 0; groups = 0;
  groupToSub = 0; bodyToSub = 0; subToBody = 0; subToGroup = 0;
  GtG = 0; GtGtilda = 0; coarseConnectGtG = 0; eqNumsGtG = 0;
  glCrnGroup = 0; 
  ngrbms = 0; firstAlpha = 0; nRowAlpha = 0;
  mpcToMpc = 0; CCtsolver = 0; 
  subsWithMpcs = 0; numSubsWithMpcs = 0; mpcSubMap = 0;
  nExtraFp = nRebuildGtG = nRebuildCCt = nLinesearchIter = nSubIterDual = nSubIterPrimal = nStatChDual = nStatChPrimal = 0;
  newton_iter = 0; lambda_total = 0;
  mpcPat = 0;
}

template<class Scalar>
int
GenFetiDPSolver<Scalar>::numRBM()
{
 if(GtGtilda) return GtGtilda->numRBM();
 else return (KccSolver) ? KccSolver->numRBM() : 0;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::wetInterfaceComms()
{
  // send and receive numWIdof
  FSCommPattern<int> *wiOnePat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
  for(int iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setWIoneCommSize(wiOnePat);
  wiOnePat->finalize();
  paralApply(this->nsub, this->sd, &BaseSub::sendNumWIdof, wiOnePat);
  wiOnePat->exchange();
  paralApply(this->nsub, this->sd, &BaseSub::recvNumWIdof, wiOnePat);
  delete wiOnePat;

  // send and receive glToLocalWImaps
  FSCommPattern<int> *wiMapPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend,
                                                        FSCommPattern<int>::NonSym);
  for(int iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setWImapCommSize(wiMapPat);
  wiMapPat->finalize();
  paralApply(this->nsub, this->sd, &BaseSub::sendWImap, wiMapPat);
  wiMapPat->exchange();
  paralApply(this->nsub, this->sd, &BaseSub::recvWImap, wiMapPat);
  delete wiMapPat;

  // create this->wiPat FSCommPattern object, used to send/receive wet this->interface interaction vectors
  this->wiPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<Scalar>::CopyOnSend,
                                    FSCommPattern<Scalar>::NonSym);
  for(int iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setWICommSize(this->wiPat);
  this->wiPat->finalize();
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::rebuildLHSfreq()
{
  // PJSA 12-15-04: this function is used to rebuild the LHS in the case of a change in the frequncy or waveNumber
  // assuming nothing else is changed (ie no new mpcs or contact for example)

  reconstruct();

  startTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);
  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::zeroLocalMatrices);
  timedParal(this->times.consMatrix, this->nsub, this, &GenFetiSolver<Scalar>::assembleMatrices);
  stopTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);

  refactor();
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::reconstruct()
{
  // 1. reset the orthoset
  this->resetOrthoSet();
 
  if(domain->solInfo().isCoupled) domain->computeCoupledScaleFactors();
  if(domain->probType() == SolverInfo::Modal) { computeRbms = false; this->fetiInfo->project_g = false; }

  if(this->fetiInfo->isEdgeAugmentationOn() && this->fetiInfo->numdir > 0) { // augmentation depends on freq/k
    // 2. rebuild augmentation Q
    filePrint(stderr," ... Rebuild Edge Augmentation (Q)  ... \n");
    if(this->fetiInfo->waveMethod != FetiInfo::averageMat) computeLocalWaveNumbers();
    paralApplyToAll(this->nsub, this->sd, &BaseSub::zeroEdgeDofSize);
    paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::makeQ);  // rebuild augmentation matrix

    // 3. reconstruct local Kcc etc. since size of Kcc may have changed
    startTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);
    paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::constructKcc);
    if(domain->solInfo().isCoupled)  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::constructKcw);
    stopTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::refactor()
{
  // 4.  stiffness scaling && edge weighting -> check !!
  if(this->fetiInfo->scaling == FetiInfo::kscaling) {
    execParal(this->nsub, this, &GenFetiSolver<Scalar>::sendScale);
    this->vPat->exchange();
    paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::collectScaling, this->vPat);
  }

  if(domain->solInfo().isCoupled) {
    paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::reScaleAndReSplitKww); // HB
  }

  if(this->fetiInfo->augment == FetiInfo::WeightedEdges)
    paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::weightEdgeGs); // PJSA: W*Q
 
  // 5. factor local matrices
  filePrint(stderr," ... Factor Subdomain Matrices      ... \n");
  startTimerMemory(this->times.factor, this->times.memoryFactor);
  if(this->fetiInfo->solvertype != FetiInfo::spooles)
    timedParal(this->times.factorMat, this->nsub, this, &GenFetiDPSolver<Scalar>::factorLocalMatrices);
  else // PJSA: unresolved problem when factoring local spooles in parallel
    for(int iSub=0; iSub<this->nsub; ++iSub) factorLocalMatrices(iSub);
  stopTimerMemory(this->times.factor, this->times.memoryFactor);

  // 6. zero and rebuild and refactor Kcc 
  if(this->fetiInfo->isEdgeAugmentationOn() && this->fetiInfo->numdir > 0) {
    filePrint(stderr, " ... Reconstruct Kcc solver         ... \n");
    delete KccSparse; KccSparse = 0; KccSolver = 0;
    makeKcc();

    delete this->wksp;
    this->times.memoryDV -= memoryUsed();
    int numC = (KccSolver) ? KccSolver->neqs() : 0;
    this->wksp = new GenFetiWorkSpace<Scalar>(this->interface, internalR, internalWI, ngrbms, numC, globalFlagCtc);
    this->times.memoryDV += memoryUsed();

    int tLocalCLen = 0;
    for(int iSub = 0; iSub < this->nsub; ++iSub) {
      internalC.domLen[iSub] = this->sd[iSub]->numCoarseDofs();
      tLocalCLen += internalC.domLen[iSub];
    }
    internalC.len = tLocalCLen;
  }
  else {
    filePrint(stderr, " ... Re-assemble Kcc solver         ... \n");
    startTimerMemory(this->times.coarse1, this->times.memoryGtG);
    if(KccSolver) KccSolver->zeroAll();
    paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::multKcc);
    for(int iSub = 0; iSub < this->nsub; ++iSub) {
      GenAssembledFullM<Scalar> *kel = this->sd[iSub]->getKcc();
      int *dofs = this->sd[iSub]->getKccDofs();
      if(KccSparse) KccSparse->add(*kel, dofs);
    }
    stopTimerMemory(this->times.coarse1, this->times.memoryGtG);
    filePrint(stderr, " ... Factorize Kcc solver           ... \n");
    startTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
#ifdef DISTRIBUTED
    if(KccSolver) KccSolver->unify(this->fetiCom);
#endif
    if(KccSolver) KccSolver->parallelFactor();
    stopTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
  }

}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::localSolveAndJump(GenDistrVector<Scalar> &fr, GenDistrVector<Scalar> &lambda,
                                           GenDistrVector<Scalar> &ur, GenVector<Scalar> &fc,
                                           GenVector<Scalar> &uc, GenDistrVector<Scalar> &r,
                                           GenDistrVector<Scalar> &fw, double rho, GenVector<Scalar> &mu)
{
  // augmented lagrangian matrix vector multiply plus rhs
  // F'*lambda + d' = (P^T*F*P + rho*G*G^T)*(lambda-lambda^e) + P^T*(d+F*lambda^e) + G*mu = P^T*(F*(P*lambda-G*(G^T*G)^-1*e)+d) + G*(mu+rho*gamma)
  if(this->fetiInfo->outerloop == FetiInfo::CGAL && ngrbms) {
    GenDistrVector<Scalar> Plambda(this->interface);
    project(lambda, Plambda, 2, true); // Plambda = P*lambda - G(G^T*G)^-1 * e
    localSolveAndJump(fr, Plambda, ur, fc, uc, r, fw); // r = F*Plambda+d
    project(r, r, 2); // r = P^T*r (note: P^T = P)
    GenVector<Scalar> beta(mu);
    GtG->forward(beta);
    if(rho != 0.0) {
      GenVector<Scalar> &gamma = this->wksp->ret_gamma(); // gamma = G^T*lambda+e (already computed)
      beta.linAdd(rho,gamma); // beta = mu + rho*gamma
    }
    multG(beta, r, 1.0, 1.0,-1); // r += G*beta
  }
  else {
    localSolveAndJump(fr, lambda, ur, fc, uc, r, fw); // r = F*lambda+d
  }
}

template<class Scalar>
Scalar
GenFetiDPSolver<Scalar>::localSolveAndJump(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &dur, GenVector<Scalar> &duc,
                                           GenDistrVector<Scalar> &Fp, double rho)
{
  // augmented lagrangian matrix vector multiply: 
  // F'p = (P^T*F*P + rho*G*G^T)*p
  Scalar pFp;
  if(this->fetiInfo->outerloop == FetiInfo::CGAL && ngrbms) {
    GenDistrVector<Scalar> Pp(this->interface);
    project(p, Pp, 2); // Pp = P*p
    localSolveAndJump(Pp, dur, duc, Fp); // Fp = F*P*p
    project(Fp, Fp, 2); // Fp = P^T*F*P*p (note: P^T = P)
    if(rho != 0.0) {
      GenVector<Scalar> beta(ngrbms);
      trMultG(p, beta, rho, 0.0,-1); // beta = rho*G^T*p
      multG(beta, Fp, 1.0, 1.0,-1);  // Fp += G*beta
    }
    pFp = p*Fp;
  }
  else {
    pFp = localSolveAndJump(p, dur, duc, Fp); // p = F*p
  }
  return pFp;
}

template<class Scalar>
Scalar
GenFetiDPSolver<Scalar>::lagrangian(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &d)
{
  // L(lambda) = 0.5*lambda^T*(r+d)
  GenDistrVector<Scalar> x(r); x += d;
  return 0.5*(lambda*x);
}

template<class Scalar>
Scalar
GenFetiDPSolver<Scalar>::lagrangian(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &r, double rho, GenVector<Scalar> &mu)
{
  GenDistrVector<Scalar> &lambda0  = this->wksp->ret_lambda0(); 
  GenDistrVector<Scalar> &r0 = this->wksp->ret_r0();
  GenVector<Scalar> &gamma0 = this->wksp->ret_gamma0();

  // L(lambda,mu,rho) = 0.5*(lambda-lambda0)^T*(r+r0) + gamma0^T(mu+0.5*rho*gamma0);
  GenDistrVector<Scalar> l(this->interface);  l.linC(1.0, lambda, -1.0, lambda0);
  GenDistrVector<Scalar> x(this->interface);  x.linC(1.0, r, 1.0, r0);
  GenVector<Scalar> b(ngrbms);                b.linC(1.0, mu, 0.5*rho, gamma0);
  return 0.5*(l*x)+gamma0*b;
}

template<class Scalar>
bool
GenFetiDPSolver<Scalar>::redundant(bool print_flag) 
{
  // check for redundancies in the active dual inequality and equality constraints
  int nr = GtGtilda->numRBM() - GtG->numRBM();
  if(nr > 0 && print_flag && this->myCPU == 0) 
    cerr << "warning: redundant constraints in active set (number = " << nr << ", grbm_tol = " << this->fetiInfo->grbm_tol << ")\n";
  return (nr > 0);
}

template<class Scalar>
double
GenFetiDPSolver<Scalar>::primalError()
{
  GenDistrVector<Scalar> &fr      = this->wksp->ret_fr();
  GenDistrVector<Scalar> &ur      = this->wksp->ret_ur();
  GenDistrVector<Scalar> &lambda  = this->wksp->ret_lambda(); 
  GenDistrVector<Scalar> &fw      = this->wksp->ret_fw();
  GenVector<Scalar> &fc           = this->wksp->ret_fc();
  GenVector<Scalar> &uc           = this->wksp->ret_uc();
  GenVector<Scalar> &alpha        = this->wksp->ret_alpha(); alpha.zero();

  GenDistrVector<Scalar> v(this->interface); double vv;
  localSolveAndJump(fr, lambda, ur, fc, uc, v, fw); // v = d + F*lambda
  nExtraFp++;
  if(ngrbms > 0) { 
    if((nStatChDual > nStatChDual_gtg) || (nStatChPrimal > nStatChPrimal_gtg)) rebuildGtGtilda(); // status change has occurred since last GtG rebuild
    tProject(v, v, vv, 3, false); // now v is the jump
  }
  double error = preCondition(v,v); // v = M^-1 * v and error is usual estimate of Ku-f
  return error; 
}

template<class Scalar>
double
GenFetiDPSolver<Scalar>::dualError(GenDistrVector<Scalar> &w, bool &proportional)
{
  if(globalFlagCtc) {
    GenDistrVector<Scalar> w_f(this->interface);
    GenDistrVector<Scalar> w_c(this->interface);
    GenDistrVector<Scalar> w_p(this->interface);
    execParal4R(this->nsub, this, &GenFetiDPSolver<Scalar>::split, w, w_f, w_c, w_p); // w_f = free gradient, w_c = chopped gradient, w_p = w_f+w_c
    //Scalar n1 = w*w, n2 = w_f*w_f, n3 = w_c*w_c, n4 = w_p*w_p;
    //cerr << "w*w = " << n1 << ", w_f*w_f = " << n2 << ", w_c*w_c = " << n3 << ", w_p*w_p = " << n4 << endl;
    proportional = (w_c.sqNorm() <= this->fetiInfo->gamma*this->fetiInfo->gamma*w_f.sqNorm());
    if(domain->solInfo().debug_icntl[3] == 1 && dualStatusChange) proportional = true; // XXXX don't do proportioning after expansion step
    if(domain->solInfo().debug_icntl[9] == 1) proportional = false; // XXXX
    domain->solInfo().debug_cntl[9] = sqrt(w_c.sqNorm()); domain->solInfo().debug_icntl[9] = 0; // XXXX
    return w_p.sqNorm();
  }
  else {
    proportional = true;
    return w.sqNorm();
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::projectActiveIneq(GenDistrVector<Scalar> &v)
{
  execParal1R(this->nsub, this, &GenFetiDPSolver<Scalar>::subProjectActiveIneq, v);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subProjectActiveIneq(int iSub, GenDistrVector<Scalar> &v)
{
  this->sd[iSub]->projectActiveIneq(v.subData(this->sd[iSub]->localSubNum()));
}

template<class Scalar>
bool
GenFetiDPSolver<Scalar>::checkInequalities(GenDistrVector<Scalar> &v, int flag)
{
  // flag = 0 : return true if N^T*lambda = 0 (ie active set constraints are satisfied) 
  // flag = 1 : return true if lambda_i <= 0 for all inequality constraints. default
  bool ret = true;
  execParal3R(this->nsub, this, &GenFetiDPSolver<Scalar>::subCheckInequalities, v, ret, flag);
#ifdef DISTRIBUTED
  ret = this->fetiCom->globalMin(int(ret));
#endif
  return ret;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subCheckInequalities(int iSub, GenDistrVector<Scalar> &v, bool &ret, int flag)
{
  bool print_debug = true;
  this->sd[iSub]->checkInequalities(v.subData(this->sd[iSub]->localSubNum()), ret, flag, print_debug);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::multG(GenVector<Scalar> &x, GenDistrVector<Scalar> &y, double alpha, double beta, int flag)
{
  // flag = 0:  y = alpha*G*x + beta*y
  // flag = 1:  y = alpha*Gtilda*x + beta*y (default)
  // flag = -1: y = alpha*G*L*x + beta*y
  if(beta == 0) y.zero(); else y *= beta;
  GenVector<Scalar> v(x);
  if(flag == -1) GtG->backward(v);
  execParal4R(this->nsub, this, &GenFetiDPSolver<Scalar>::subMultG, v, y, alpha, int(flag==1)); // y += alpha*G*v
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subMultG(int iSub, GenVector<Scalar> &x, GenDistrVector<Scalar> &y, double alpha, int flag)
{
  this->sd[iSub]->multG(x, y.subData(this->sd[iSub]->localSubNum()), alpha, flag); 
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::trMultG(GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha, double beta, int flag)
{
  // flag = 0:  y = alpha*G^T*x + beta*y
  // flag = 1:  y = alpha*Gtilda^T*x + beta*y (default)
  // flag = -1: y = alpha*L^T*G^T*x + beta*y
  if(beta == 0) y.zero(); else y *= beta;
  GenVector<Scalar> v(ngrbms, 0.0);
  execParal4R(nGroups1, this, &GenFetiDPSolver<Scalar>::subTrMultG, x, v, alpha, int(flag==1)); // v += alpha*G^T*x
#ifdef DISTRIBUTED
  this->fetiCom->globalSum(ngrbms, v.data());
#endif
  if(flag == -1) GtG->forward(v);
  y += v;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subTrMultG(int iGroup, GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha, int flag)
{
#ifdef DISTRIBUTED
  for(int i=0; i<this->nsub; ++i) {
    if(this->sd[i]->group == groups[iGroup])
      this->sd[i]->trMultG(x.subData(this->sd[i]->localSubNum()), y, alpha, flag); 
  }
#else
  int *grsubs = (*groupToSub)[iGroup];
  for(int i = 0; i < groupToSub->num(iGroup); ++i) {
    int iSub = grsubs[i];
    this->sd[iSub]->trMultG(x.subData(this->sd[iSub]->localSubNum()), y, alpha, flag); 
  }
#endif
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::reSolveGtG(GenVector<Scalar> &x, int flag)
{
  // flag = 0:  x = (G^T*G)^-1 * x
  // flag = 1:  x = (Gtilda^T*Gtilda)^-1 * x
  if(flag == 0) GtG->reSolve(x);
  else if(flag == 1) GtGtilda->reSolve(x);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::project(GenDistrVector<Scalar> &z, GenDistrVector<Scalar> &y, int project_level, bool eflag)
{
  // project_level = 1 (i), 2 (g), 3 (i+g) ... i is the inequality constraints from contact, g is the equality constraints from grbms
  startTimerMemory(this->times.project, this->times.memoryProject1);
  bool gproj = (ngrbms) && (this->fetiInfo->project_g || project_level > 1);
  bool iproj = (globalFlagCtc && (project_level != 2));
  // y = P*z + y0
  // gproj && iproj :  P = P_i*(I - G*(G^T*P_i*G)^(-1)*G^T)*P_i = P_i*(I - G*(Gtilda^T*Gtilda)^(-1)*Gtilda^T), y0 = (eflag) ? Gtilda*(Gtilda^T*Gtilda)^(-1)*e : 0 
  // gproj && !iproj:  P = I - G*(G^T*G)^(-1)*G^T, y0 = (eflag) ? G*(G^T*G)^(-1)*e : 0
  // iproj && !gproj:  P = P_i, y0 = 0
  // note: P_i depends on active set which is a function of y(ie lambda) hence dual planing subiterations required for feasible projection

  // 0. y = z
  if(&y != &z) y = z;

  if(gproj) {
    // 1. alpha = -Gtilda^T*z
    GenVector<Scalar> alpha(ngrbms);
    trMultG(z, alpha, -1.0, 0.0, int(iproj));
    if(eflag) alpha -= e_copy;

    // 2. alpha = (Gtilda^T*Gtilda)^(-1)*alpha
    reSolveGtG(alpha, int(iproj));

    // 3. y = z + G*alpha
    multG(alpha, y, 1.0, 1.0);
  }

  // 4. y = P_i*y
  if(iproj) projectActiveIneq(y);

  // 5. dual planing: if y (lambda) is not feasible then expand active set and re-project (recursive)
  if(iproj && eflag && updateActiveSet(y, 0))
    { stopTimerMemory(this->times.project, this->times.memoryProject1); project(y, y, project_level, true); nSubIterDual++; return; }

  stopTimerMemory(this->times.project, this->times.memoryProject1);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::tProject(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &w, double &dual_error, int project_level, bool pflag)
{
  // project_level = 1 (i), 2 (g), 3 (i+g) ... i is the inequality constraints from contact, g is the equality constraints from grbms
  startTimerMemory(this->times.project, this->times.memoryProject1);
  bool gproj = (ngrbms) && (this->fetiInfo->project_g || project_level > 1);
  bool iproj = (globalFlagCtc && (project_level != 2));
  bool proportional;
  // w = P^T*r
  // gproj && iproj :  P^T = P = P_i*(I - G*(G^T*P_i*G)^(-1)*G^T)*P_i = P_i*(I - G*(Gtilda^T*Gtilda)^(-1)*Gtilda^T)
  // gproj && !iproj:  P^T = P = I - G*(G^T*G)^(-1)*G^T 
  // iproj && !gproj :  P^T = P = P_i 
  // note: P_i depends on active set which is a function of r+G*alpha, hence primal planing subiterations required if not proportional
  // this function also saves alpha in the feti workspace

  // 0. w = r
  if(&w != &r) w = r;

  if(gproj) {
    // 1. alpha = -Gtilda^t*r
    GenVector<Scalar> &alpha = this->wksp->ret_alpha();
    trMultG(r, alpha, -1.0, 0.0, int(iproj));

    // 2. alpha = (Gtilda^T*Gtilda)^(-1)*alpha
    reSolveGtG(alpha, int(iproj));

    // 3. w = r + G*alpha
    multG(alpha, w, 1.0, 1.0);
  }

  // 4. primal planing: if not r+G*alpha is not proportional then contract active set and re-project (recursive)
  if(pflag) {
    dual_error = dualError(w, proportional);
    if(iproj && !proportional && updateActiveSet(w,1))
      { stopTimerMemory(this->times.project, this->times.memoryProject1); tProject(w, w, dual_error, project_level, pflag); nSubIterPrimal++; return; }
  }

  // 5. w = P_i*w
  if(iproj) projectActiveIneq(w);

  stopTimerMemory(this->times.project, this->times.memoryProject1);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::computeL0(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f)
{
  lambda.zero();
  if(newton_iter > 0) return; // nonlinear
  if(this->fetiInfo->initial_lambda == 1) { // PJSA 8-23-2007 new initialization for lambda from Gosselet, Rey and Rixen (CMAME 2003)
    GenDistrVector<Scalar> &fr = this->wksp->ret_fr();
    execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::computeL00, lambda, fr);
    this->vPat->exchange();
    execParal1R(this->nsub, this, &GenFetiSolver<Scalar>::interfaceDiff, lambda);
    if(mpcPrecon) cctSolveMpc(lambda);
  }
  if(this->glNumMpc) {
    if(ngrbms > 0 /*&& newton_iter == 0*/) makeE(f); // compute e = R^T*f (first time only for nonlinear)
    int project_level = (domain->solInfo().debug_icntl[1] == 1 && ngrbms) ? 3 : -1;  
    //if(newton_iter > 0) lambda += (*lambda_total); // nonlinear
    project(lambda, lambda, project_level, true); // lambda = P*lambda + G*(G^T*G)^(-1) *(G^T*lambda+e)  (recursive)
    //if(newton_iter > 0) lambda -= (*lambda_total); // nonlinear
    if(globalFlagCtc) feasible(lambda, true, (this->fetiInfo->outerloop != 3)); // check if lambda is feasible
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::computeL00(int iSub, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &fr)
{ 
  this->sd[iSub]->computeL00(lambda.subData(this->sd[iSub]->localSubNum()), fr.subData(this->sd[iSub]->localSubNum()));
  this->sd[iSub]->sendInterf(lambda.subData(this->sd[iSub]->localSubNum()), this->vPat);
} 


template<class Scalar>
double
GenFetiDPSolver<Scalar>::computeFNorm()
{
  int i;
  double Fnorm = 1.0, Fnorm_prev=1.0;    //CRW
  GenDistrVector<Scalar> &p = this->wksp->ret_p(), &Fp = this->wksp->ret_Fp(), &dur = this->wksp->ret_du();
  GenVector<Scalar> &duc = this->wksp->ret_duc();
  for(i=0; i<this->fetiInfo->max_power_iter; ++i) {
    if(i==0) initRand(p); else p = Fp;
    Scalar pp = p*p;
    Scalar pFp = localSolveAndJump(p, dur, duc, Fp, 0.0); nExtraFp++;
    Fnorm = ScalarTypes::Real(pFp/pp);
    if(i>0 && (Fnorm-Fnorm_prev)/Fnorm < this->fetiInfo->power_iter_tol) break; else Fnorm_prev = Fnorm;
  }
  if(verboseFlag && this->myCPU == 0) cerr << "||F|| estimate = " << Fnorm << " after power iteration = " << i+1 
                                           << " of " << this->fetiInfo->max_power_iter << endl;
  return Fnorm;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::initRand(GenDistrVector<Scalar> &p)
{
  p.initRand(); 
  execParal1R(this->nsub, this, &GenFetiSolver<Scalar>::interfSend, p); 
  this->vPat->exchange(); 
  execParal1R(this->nsub, this, &GenFetiSolver<Scalar>::interfaceDiff, p);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::normalizeC()
{
  GenVector<Scalar> cnorm(this->glNumMpc, 0.0);
  for(int i=0; i<this->nsub; ++i) this->sd[i]->normalizeCstep1(cnorm.data());
#ifdef DISTRIBUTED
  this->fetiCom->globalSum(this->glNumMpc, cnorm.data());
#endif
  for(int i=0; i<this->glNumMpc; ++i) cnorm[i] = ScalarTypes::sqrt(cnorm[i]);
  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::normalizeCstep2, cnorm.data());
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::trMultC(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f)
{
  // compute f = C^T*lambda
  f.zero();
  execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::subTrMultC, lambda, f);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subTrMultC(int iSub, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f)
{ 
  this->sd[iSub]->multAddCT(lambda.subData(this->sd[iSub]->localSubNum()), f.subData(this->sd[iSub]->localSubNum())); 
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::multC(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &cu)
{
  // compute cu = C*u
  cu.zero();
  execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::subMultC, u, cu);

  execParal1R(this->nsub, this, &GenFetiSolver<Scalar>::interfSend, cu);
  this->vPat->exchange();
  execParal1R(this->nsub, this, &GenFetiSolver<Scalar>::interfaceDiff, cu);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subMultC(int iSub, GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &cu)
{
  this->sd[iSub]->multC(u.subData(this->sd[iSub]->localSubNum()), cu.subData(this->sd[iSub]->localSubNum()));
}

