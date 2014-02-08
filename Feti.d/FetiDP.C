#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <stdexcept>

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
                  Rbm **, bool _rbmFlag, bool _geometricRbms)
 : GenFetiSolver<Scalar>(_nsub, threadManager->numThr()), internalR(_nsub), internalC(_nsub), internalWI(_nsub) 
{
/*
 if(geoSource->isShifted())
   filePrint(stderr, " ... FETI-DPH Solver is Selected    ...\n");
 else
   filePrint(stderr, " ... FETI-DP Solver is Selected     ...\n");
*/
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

 rbmFlag = _rbmFlag; // if this is true then check for singularities in Kcc^* and each Krr^{(s)}
                     // and deal with the rigid body modes using projection method
 geometricRbms = _geometricRbms; // if this is true then use geometric method
                                 // to compute the rigid body modes of Kcc^* and Krr^{(s)}
                                 // when rbmFlag is true, otherwise use algebraic null space   

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

 if(sysMatrices != 0) {
   for(iSub = 0; iSub < this->nsub; ++iSub) {
     this->sd[iSub]->Krr = sysMatrices[iSub];
     this->sd[iSub]->KrrSparse = sysSparse[iSub];
     if(this->fetiInfo->printMatLab) {
       std::stringstream filename;
       filename << "localmat" << this->sd[iSub]->subNum();
       this->sd[iSub]->KrrSparse->printSparse(filename.str());
     }
     this->fetiOps[iSub]->setSysMatrix(this->sd[iSub]->Krr, sysSparse[iSub]);
   }
 }

 if(verboseFlag) filePrint(stderr," ... Build Edge Augmentation (Q)    ... \n");
 computeLocalWaveNumbers();
 paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::makeQ);  // build augmentation matrix
 if(this->fetiInfo->augment == FetiInfo::Gs) {
   // exchange number of each neighbors rbms
   paralApply(this->nsub, this->sd, &BaseSub::sendNumNeighbGrbm, this->sPat);
   this->sPat->exchange();
   paralApply(this->nsub, this->sd, &BaseSub::recvNumNeighbGrbm, this->sPat);
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
 if(verboseFlag) filePrint(stderr," ... Factor Subdomain Matrices      ... \n");
 startTimerMemory(this->times.factor, this->times.memoryFactor);
 if(this->fetiInfo->solvertype != FetiInfo::spooles && this->fetiInfo->solvertype != FetiInfo::mumps) // spooles/mumps factor is apparently not thread-safe
   timedParal(this->times.factorMat, this->nsub, this, &GenFetiDPSolver<Scalar>::factorLocalMatrices);
 else 
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
 if(ngrbms) makeGtG();  // currently G = C^T*R (ie restriction of R to mpc interface)

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

 paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::initMpcStatus);

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
 startTimerMemory(this->times.coarse1, this->times.memoryGtG);

 int i, iSub;
 int glNumSub = this->subToSub->csize();

 // STEP 1. count number of corner nodes and make subToCorner connectivity
 Connectivity *subToCorner = 0;
 if(!cornerToSub) {
   int *pointer = new int[glNumSub+1];
   for(i=0; i<glNumSub+1; ++i) pointer[i] = 0;
   for(iSub=0; iSub<this->nsub; ++iSub)
     pointer[this->sd[iSub]->subNum()] = this->sd[iSub]->numCorners();
#ifdef DISTRIBUTED
   this->fetiCom->globalSum(glNumSub, pointer);
#endif
   int total = 0;
   for(iSub=0; iSub < glNumSub; ++iSub) {
     int tmp = pointer[iSub];
     pointer[iSub] = total;
     total += tmp;
   }
   pointer[glNumSub] = total;

   int *glCornerNodes = new int[total];
   for(i=0; i<total; ++i) glCornerNodes[i] = 0;
   for(iSub=0; iSub<this->nsub; ++iSub) {
     int numCorner = this->sd[iSub]->numCorners();
     int *localCornerNodes = this->sd[iSub]->getLocalCornerNodes();
     int *glN = this->sd[iSub]->getGlNodes();
     for(int iCorner=0; iCorner<numCorner; ++iCorner)
       glCornerNodes[pointer[this->sd[iSub]->subNum()]+iCorner] = glN[localCornerNodes[iCorner]];
   }
#ifdef DISTRIBUTED
   this->fetiCom->globalSum(total, glCornerNodes);
#endif
   glNumCorners = 0;
#ifdef TFLOP
  map<int, int, less<int> > glCornerMap;
#else
   map<int, int> glCornerMap;
#endif
   for(int iCorner=0; iCorner<total; ++iCorner)
     if(glCornerMap.find(glCornerNodes[iCorner]) == glCornerMap.end() )
        glCornerMap[ glCornerNodes[iCorner] ] = glNumCorners++;
   delete [] glCornerNodes;
   if(verboseFlag) filePrint(stderr," ... Total Number of Corners %5d  ...\n", glNumCorners);

   int *target = new int[total];
   for(i=0; i<total; ++i) target[i] = 0;
   for(iSub=0; iSub<this->nsub; ++iSub) {
     int numCorner    = this->sd[iSub]->numCorners();
     int *cornerNodes = this->sd[iSub]->getCornerNodes();
     int *localCornerNodes = this->sd[iSub]->getLocalCornerNodes();
     int *glN = this->sd[iSub]->getGlNodes();
     for(int iCorner=0; iCorner<numCorner; ++iCorner) {
        cornerNodes[iCorner] = glCornerMap[glN[localCornerNodes[iCorner]]];
        target[iCorner+pointer[this->sd[iSub]->subNum()]] = cornerNodes[iCorner];
     }
   }
#ifdef DISTRIBUTED
   this->fetiCom->globalSum(total, target);
#endif

   subToCorner = new Connectivity(glNumSub, pointer, target);
   cornerToSub = subToCorner->reverse();
 }
 else {
   subToCorner = cornerToSub->reverse();
 }

 // STEP 2. make subToCoarse connectivity (includes primal mpcs and augmentation) and coarseConnectivity
 Connectivity *coarseToSub = cornerToSub;
 mpcOffset = coarseToSub->csize();
 if(this->glNumMpc_primal > 0) {
   coarseToSub = coarseToSub->merge(this->mpcToSub_primal);
 }
 augOffset = coarseToSub->csize();
 switch(this->fetiInfo->augment) {
   case FetiInfo::Gs: {
     Connectivity *augcoarseToSub = coarseToSub->merge(this->subToSub);
     if(coarseToSub != cornerToSub) delete coarseToSub;
     coarseToSub = augcoarseToSub;
   } break;
   case FetiInfo::WeightedEdges:
   case FetiInfo::Edges: {
     if(!this->edgeToSub) makeEdgeConnectivity();
     Connectivity *augcoarseToSub = coarseToSub->merge(this->edgeToSub);
     if(coarseToSub != cornerToSub) delete coarseToSub;
     coarseToSub = augcoarseToSub;
   } break;
 }
 Connectivity *subToCoarse = (coarseToSub != cornerToSub) ? coarseToSub->reverse() : subToCorner;
 Connectivity *coarseConnectivity = coarseToSub->transcon(subToCoarse);

 // STEP 3. make the coarse problem equation numberer
 int renumFlag = (this->fetiInfo->gtgSolver == FetiInfo::skyline || this->fetiInfo->gtgSolver == FetiInfo::blocksky) ? domain->solInfo().renum : 0;
 compStruct renumber = coarseConnectivity->renumByComponent(renumFlag);
 if(cornerEqs) delete cornerEqs; // PJSA
 cornerEqs = new DofSetArray(coarseConnectivity->csize(), renumber.renum, 1);
 delete [] renumber.xcomp;

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
 delete subToCorner;

 int *glCornerDofs = new int[glNumCorners];
 for(i = 0; i < glNumCorners; ++i) glCornerDofs[i] = 0;
 for(iSub = 0; iSub < this->nsub; ++iSub) this->sd[iSub]->markCornerDofs(glCornerDofs);
#ifdef DISTRIBUTED
 this->fetiCom->globalMpiOp(glNumCorners, glCornerDofs, MPI_BOR); // MPI_BOR is an mpi bitwise or
#endif
 for(i = 0; i < glNumCorners; ++i) cornerEqs->mark(i, glCornerDofs[i]);
 delete [] glCornerDofs;

 if(this->glNumMpc_primal > 0) {
   for(int iMPC=0; iMPC<this->glNumMpc_primal; ++iMPC)
     cornerEqs->setWeight(mpcOffset+iMPC, 1);
   this->times.numMPCs = this->glNumMpc_primal;
 }

 if(this->fetiInfo->augment == FetiInfo::Gs) {
   this->times.numRBMs = 0;
   int *numRBMPerSub = new int[glNumSub];
   for(iSub=0; iSub<glNumSub; ++iSub) numRBMPerSub[iSub] = 0;
   for(iSub=0; iSub<this->nsub; ++iSub)
     numRBMPerSub[this->sd[iSub]->subNum()] = this->sd[iSub]->numRBM();
#ifdef DISTRIBUTED
   this->fetiCom->globalSum(glNumSub,numRBMPerSub);
#endif
   for(iSub=0; iSub<glNumSub; ++iSub) {
     cornerEqs->setWeight(augOffset+iSub, numRBMPerSub[iSub]);
     this->times.numRBMs += numRBMPerSub[iSub];
   }
   delete [] numRBMPerSub;
 }

 if(this->fetiInfo->isEdgeAugmentationOn()) {
   this->times.numEdges = 0;
   int *edgeWeights = new int[this->edgeToSub->csize()];
   for(i=0; i<this->edgeToSub->csize(); ++i) edgeWeights[i] = 0;
   for(iSub=0; iSub<this->nsub; ++iSub) {
     int numNeighbor = this->sd[iSub]->numNeighbors();
     int myNum = this->sd[iSub]->subNum();
     int jEdgeN = 0;
     for(int jSub=0; jSub<numNeighbor; ++jSub) {
       int subJ = this->sd[iSub]->getSComm()->subNums[jSub];
       if(this->sd[iSub]->isEdgeNeighbor(jSub)) {
         if(myNum < subJ) { // PJSA
           int nEdge = this->sd[iSub]->numEdgeDofs(jSub);
           edgeWeights[(*this->subToEdge)[this->sd[iSub]->subNum()][jEdgeN]] = nEdge;
         }
         jEdgeN++;
       }
     }
   }
#ifdef DISTRIBUTED
   this->fetiCom->globalSum(this->edgeToSub->csize(), edgeWeights);
#endif
   for(int iEdge=0; iEdge<this->edgeToSub->csize(); ++iEdge) {
     cornerEqs->setWeight(augOffset+iEdge, edgeWeights[iEdge]);
     this->times.numEdges += edgeWeights[iEdge];
   }
   delete [] edgeWeights;
 }

 cornerEqs->finish();
 this->times.numCRNs = cornerEqs->size();

 if(verboseFlag) {
   filePrint(stderr, " ... Size of Interface  %7d     ...\n", this->interface.len);
   filePrint(stderr, " ... Size of Global Kcc %7d     ...\n", cornerEqs->size());
   filePrint(stderr, " ... global_rbm_tol     = %9.3e ...\n", this->fetiInfo->grbm_tol);
   filePrint(stderr, " ... global_cor_rbm_tol = %9.3e ...\n", this->fetiInfo->crbm_tol);
 }

 // STEP 4. make the rigid body modes and condensed equation numberer
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
 if(this->sd && this->nsub > 0) {
   groups[0] = (*subToGroup)[this->sd[0]->subNum()][0];
   int n = 1;
   for(int i = 1; i < this->nsub; ++i) {
     int group = (*subToGroup)[this->sd[i]->subNum()][0];
     int j;
     for(j = 0; j < n; ++j) if(group == groups[j]) j = n+1;
     if(j == n) groups[n++] = group;
   }
   nGroups1 = n;  // number of groups represented on this processor
 }
 else nGroups1 = 0;
#else
 for(int i = 0; i < nGroups; ++i) groups[i] = i; 
 nGroups1 = nGroups;  
#endif
 if(ngrbmGr) delete [] ngrbmGr;
 ngrbmGr = new int[nGroups];
 for(int i = 0; i < nGroups; ++i) ngrbmGr[i] = 0;
 if(verboseFlag) filePrint(stderr, " ... Number of bodies = %3d         ...\n", nGroups);
 
 int ngrbm = 0;
 ngrbms = 0;
 if(rbmFlag && geometricRbms) {
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
   for(int i = 0; i < nGroups; ++i) {
     nNodes[i] = 0.0;
     for(int j = 0; j < 3; ++j) centroid[3*i+j] = 0.0;
   }
   for(int i = 0; i < this->nsub; ++i) this->sd[i]->addNodeXYZ(centroid, nNodes);  // groups could be done in parallel
#ifdef DISTRIBUTED
   this->fetiCom->globalSum(nGroups, nNodes);
   this->fetiCom->globalSum(nGroups*3, centroid);
#endif
   for(int i = 0; i < nGroups; ++i) {
     if(nNodes[i] > 0) 
       for(int j = 0; j < 3; ++j) 
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
   int zColDim1 = (this->sd && this->nsub > 0) ? this->sd[0]->zColDim() : 0;  // (6 for 3D, 3 for 2D)
#ifdef DISTRIBUTED
   zColDim1 = this->fetiCom->globalMax(zColDim1);  // enforce it to be the same
#endif
   for(int i = 0; i < nGroups; ++i) {
     zRowDim[i] = 0; 
     if(mbgflag) {
       int cOff = 0;
       for(int j = 0; j < groupToBody->num(i); ++j) {
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
   for(int iSub = 0; iSub < this->nsub; ++iSub) {
     int subGroup = (*subToGroup)[this->sd[iSub]->subNum()][0];
     zRowDim[subGroup] += this->sd[iSub]->zRowDim();
   }
   if(this->glNumMpc_primal > 0) {
     for(int i = 0; i < nGroups; ++i) zRowDim[i] += groupToMpc->num(i);
   }
   if(groupToBody) delete groupToBody;

#ifdef DISTRIBUTED
   int *zRowOffset = new int[this->numCPUs*nGroups];
   for(int i = 0; i < this->numCPUs*nGroups; ++i) zRowOffset[i] = 0;
   for(int i = 0; i < nGroups1; ++i) {
     int iGroup = groups[i];
     for(int j = this->myCPU+1; j < this->numCPUs; ++j) zRowOffset[iGroup*this->numCPUs +j] = zRowDim[iGroup];
   }
   this->fetiCom->globalSum(nGroups, zRowDim);
   this->fetiCom->globalSum(this->numCPUs*nGroups, zRowOffset);
   for(int i = 0; i < nGroups; ++i) zRow[i] = zRowOffset[i*this->numCPUs + this->myCPU];
   delete [] zRowOffset;
#else
   for(int i = 0; i < nGroups; ++i) zRow[i] = 0;
#endif
   for(int i = 0; i < nGroups; ++i) {
     globalZstar[i] = new FullM(zRowDim[i], zColDim[i]);
     globalZstar[i]->zero();
   }
   // could do this in parallel (by groups)
   for(int iSub = 0; iSub < this->nsub; ++iSub) {
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
 
   int *groupProc = new int[nGroups];
#ifdef DISTRIBUTED
   for(int i = 0; i < nGroups; ++i) {
     this->fetiCom->globalSum(zRowDim[i]*zColDim[i], globalZstar[i]->data());
     groupProc[i] = -1;
   }
   for(int i = 0; i < nGroups1; ++i) groupProc[groups[i]] = this->myCPU;
   for(int i = 0; i < nGroups; ++i) groupProc[i] = this->fetiCom->globalMax(groupProc[i]);
#else
   for(int i = 0; i < nGroups; ++i) groupProc[i] = this->myCPU;
#endif
 
   // now do svd on globalZstar for each group to get globalQ for each group
   ngrbm = 0;  // total of all groups
   FullM  **Qtranspose;
   Qtranspose = new FullM * [nGroups];
   for(int i = 0; i < nGroups1; ++i) {
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
 
   delete [] groupProc;
   delete [] zRow;
   delete [] zRowDim;
   delete [] zColDim;
   for(int i = 0; i < nGroups; ++i) delete globalZstar[i];
   delete [] globalZstar;
 
   // make local Rstar
   paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::makeLocalRstar, Qtranspose);
   for(int i = 0; i < nGroups1; ++i) delete Qtranspose[groups[i]];
   delete [] Qtranspose;
 }
 if(rbmFlag && !geometricRbms && this->fetiInfo->corners == FetiInfo::noCorners) {
   paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::useKrrNullspace);
   ngrbmGr = new int[nGroups];
   ngrbm = 0;
   for(int i = 0; i < nGroups; ++i) ngrbmGr[i] = 0;
   for(int iSub = 0; iSub < this->nsub; ++iSub) {
     int subGroup = (*subToGroup)[this->sd[iSub]->subNum()][0];
     ngrbmGr[subGroup] = this->sd[iSub]->Krr->numRBM();
     ngrbm += ngrbmGr[subGroup];
   }
#ifdef DISTRIBUTED
   ngrbms = this->fetiCom->globalSum(ngrbm);  // total number of rigid body modes for all processes
#else
   ngrbms = ngrbm;
#endif
 }
 
 if(ngrbms) {
#ifdef DISTRIBUTED
   this->fetiCom->globalSum(nGroups, ngrbmGr);
#endif
   paralApply(this->nsub, this->sd, &BaseSub::setNumGroupRBM, ngrbmGr);
   paralApply(this->nsub, this->sd, &BaseSub::deleteLocalRBMs);
 }

 // end new code ********************************************************

 KccSparse = 0;

 if(cornerEqs->size() > 0) {

   double tolerance = this->fetiInfo->crbm_tol;  // this is set in input file using global_cor_rbm_tol
   GenBLKSparseMatrix<Scalar> *BLKMatrix = 0;
   GenBlockSky<Scalar> *blockSky = 0;
   GenSkyMatrix<Scalar> *sky = 0;
#ifdef DISTRIBUTED
   int firstAlpha,nRowAlpha;  
   int neq       = cornerEqs->size(); 
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
         sky = new GenDistSky<Scalar>(coarseConnectivity, cornerEqs, tolerance, 
                                      firstAlpha, nRowAlpha);     
         this->times.memoryGtGDelete = 8*sky->size();
       } else
#endif 
       // TODO: pass grbms...
       sky = new GenSkyMatrix<Scalar>(coarseConnectivity, cornerEqs, tolerance, domain->solInfo().coarseScaled); 
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
         blockSky = new GenDistBlockSky<Scalar>(coarseConnectivity, cornerEqs, tolerance,
                                                firstAlpha, nRowAlpha);   
         this->times.memoryGtGDelete = 8*blockSky->size();
       } else
#endif
       blockSky = new GenBlockSky<Scalar>(coarseConnectivity, cornerEqs, tolerance);  
       blockSky->zeroAll();
       KccSparse = blockSky;
       KccSolver = blockSky;
     }
     break;
     case FetiInfo::sparse: {
       int sparse_ngrbms = (geometricRbms) ? ngrbms : 0; // TODO pass Rbm object, not just ngrbms
#ifdef DISTRIBUTED
       if(this->subToSub->csize() == this->numCPUs && 
          this->fetiInfo->type != FetiInfo::nonlinear && 
          !domain->solInfo().doEigSweep && !domain->solInfo().doFreqSweep) {
         BLKMatrix = new GenDistBLKSparse<Scalar>(coarseConnectivity, cornerEqs, 
                                                  tolerance, firstAlpha, nRowAlpha);   
         this->times.memoryGtGDelete = 8*BLKMatrix->size();
       } else
#endif
       BLKMatrix = new GenBLKSparseMatrix<Scalar>(coarseConnectivity, cornerEqs,
                                                  tolerance, domain->solInfo().sparse_renum, sparse_ngrbms);
       BLKMatrix->zeroAll();
       KccSparse = BLKMatrix;
       KccSolver = BLKMatrix;
     }
     break;
#ifdef USE_EIGEN3
     case FetiInfo::ldlt: {
       GenEiSparseMatrix<Scalar, Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> > *cholsolver = new
         GenEiSparseMatrix<Scalar,Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >(coarseConnectivity, cornerEqs);
       KccSparse = cholsolver;
       KccSolver = cholsolver;
     }
     break;
#endif
     case FetiInfo::spooles: {
       GenSpoolesSolver<Scalar> *spsolver = new GenSpoolesSolver<Scalar>(coarseConnectivity, cornerEqs);
       KccSparse = spsolver;
       KccSolver = spsolver;
     }
     break;
     case FetiInfo::mumps: {
       GenMumpsSolver<Scalar> *mumpsolver = new GenMumpsSolver<Scalar>(coarseConnectivity, cornerEqs, (int *)0, this->fetiCom);
       KccSparse = mumpsolver;
       KccSolver = mumpsolver;
     }
     break;
   }
   this->times.memoryGtGsky += memoryUsed();

   // assemble the coarse problem: Kcc^* -> Kcc - Krc^T Krr^-1 Krc
   if(verboseFlag) filePrint(stderr, " ... Assemble Kcc solver            ...\n");
   t5 -= getTime();
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::multKcc); // create the local Kcc^*
   t5 += getTime();

   t0 -= getTime();
   paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::makeKccDofs, cornerEqs, augOffset, this->subToEdge, mpcOffset);
   if(KccSparse) for(iSub = 0; iSub < this->nsub; ++iSub) this->sd[iSub]->assembleKccStar(KccSparse); // assemble local Kcc^* into global Kcc^*
   t0 += getTime();

   // Factor coarse solver
   startTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
#ifdef DISTRIBUTED
   if(verboseFlag) filePrint(stderr, " ... Unify Kcc                      ...\n");
   KccSolver->unify(this->fetiCom);
#endif
   if(this->fetiInfo->printMatLab) {
     KccSparse->printSparse("coarsemat");
   }

   if(verboseFlag) filePrint(stderr, " ... Factor Kcc solver              ...\n");
   KccSolver->setPrintNullity(this->fetiInfo->contactPrintFlag && this->myCPU == 0);
   KccSolver->parallelFactor();
   stopTimerMemory(this->times.pfactor, this->times.memoryGtGsky);

   if(rbmFlag && geometricRbms && this->myCPU == 0) {
     if(KccSolver->neqs() > 0 && KccSolver->numRBM() != ngrbms) {
       cerr << " *** WARNING: number of singularities in Kcc (" << KccSolver->numRBM() << ")" << endl
            << "     does not match the number of Geometric RBMs (" << ngrbms << ")" << endl
            << " *** try adjusting global_cor_rbm_tol or use TRBM method" << endl;
     }
   } 
 
   if(rbmFlag && !geometricRbms && (ngrbms = KccSolver->numRBM()) > 0) {
     kccrbms = new Scalar[KccSolver->neqs()*KccSolver->numRBM()];
     KccSolver->getNullSpace(kccrbms);
     if(this->fetiInfo->nullSpaceFilterTol > 0.0) {
       for(int i= 0; i<KccSolver->numRBM(); ++i)
         for(int j=0; j<KccSolver->neqs(); ++j)
           if(ScalarTypes::norm(kccrbms[i*KccSolver->neqs()+j]) < this->fetiInfo->nullSpaceFilterTol) kccrbms[i*KccSolver->neqs()+j] = 0.0; // FILTER
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
GenFetiDPSolver<Scalar>::updateActiveSet(GenDistrVector<Scalar> &v, int flag, double tol)
{
  // flag = 0 dual planing
  // flag = 1 primal planing
  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::saveMpcStatus2);

  bool *local_status_change = new bool[this->nsub];
  execParal4R(this->nsub, this, &GenFetiDPSolver<Scalar>::subUpdateActiveSet, v, tol, flag, local_status_change);
  bool status_change1 = false;
  for(int i=0; i<this->nsub; ++i) if(local_status_change[i]) { status_change1 = true; break; }
#ifdef DISTRIBUTED
  status_change1 = this->fetiCom->globalMax((int) status_change1);
#endif

  if(status_change1) {
    paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::sendMpcStatus, mpcPat, flag);
    mpcPat->exchange();
    //paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::recvMpcStatus, mpcPat, flag);
    execParal3R(this->nsub, this, &GenFetiDPSolver<Scalar>::subRecvMpcStatus, mpcPat, flag, local_status_change);
    bool status_change2 = false;
    for(int i=0; i<this->nsub; ++i) if(local_status_change[i]) { status_change2 = true; break; }
#ifdef DISTRIBUTED
    status_change2 = this->fetiCom->globalMax((int) status_change2);
#endif
    if(status_change2 && ngrbms) rebuildGtGtilda();
    if(this->fetiInfo->contactPrintFlag && this->myCPU == 0) cerr << " ";
    delete [] local_status_change;
    return status_change2;
  }
  else {
    delete [] local_status_change;
    return false;
  }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subUpdateActiveSet(int iSub, GenDistrVector<Scalar> &lambda, double tol, int flag, bool *statusChange)
{
  this->sd[iSub]->updateActiveSet(lambda.subData(this->sd[iSub]->localSubNum()), tol, flag, statusChange[iSub]);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subRecvMpcStatus(int iSub, FSCommPattern<int> *mpcPat, int flag, bool *statusChange)
{
  this->sd[iSub]->recvMpcStatus(mpcPat, flag, statusChange[iSub]);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::update(Scalar nu, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &p, 
                                GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fp, 
                                GenDistrVector<Scalar> &ur, GenDistrVector<Scalar> &dur, 
                                GenVector<Scalar> &uc, GenVector<Scalar> &duc)
{
  if(globalFlagCtc) saveStep();
  GenDistrVector<Scalar> &lambda_k = this->wksp->ret_lambda_copy();
  GenDistrVector<Scalar> &r_k = this->wksp->ret_r_copy();

  int i;
  for(i = 0; i < this->fetiInfo->linesearch_maxit+1; ++i, ++nLinesearchIter) {

    // update lambda
    lambda.linAdd(nu,p);
    if(globalFlagCtc) project(lambda, lambda, true);

    // Update residual (r)
    if(dualStatusChange) {
      p.linC(1.0/nu, lambda, -1.0/nu, lambda_k); // reduced search direction p = (lambda-lambda_copy)/nu
      localSolveAndJump(p, dur, duc, Fp); // recompute Fp using reduced search direction
    }
    r.linAdd(nu, Fp); 

    if(i == 0 && !dualStatusChange) break; // CG step
    else { // gradient projection step
      Scalar rp = r_k*p;
      Scalar pFp = p*Fp; 
      //if(ScalarTypes::lessThan(pFp, 0)) throw std::runtime_error("FETI operator is not positive semidefinite");
      Scalar delta_f = nu*nu/2.0*pFp + nu*rp;
      if(this->fetiInfo->contactPrintFlag >= 2 && this->myCPU == 0)
        cerr << " linesearch: iteration = " << i << ", delta_f = " << delta_f << ", pFp = " << pFp << ", nu = " << nu << ", rp = " << rp << endl;
      if(ScalarTypes::lessThanEq(delta_f, 0)) break; // sequence is monotonic (note: check for gcr and gmres)
      else {
        if(i < this->fetiInfo->linesearch_maxit) { 
          restoreStep();
          nu *= this->fetiInfo->linesearch_tau;
        }
        else {
          throw std::runtime_error("linesearch did not converge");
        }
      }
    }
  }

  stepLengthChange = (i > 0);
  if(stepLengthChange) nLinesearch++;
  
  // Update primal (ur, uc)
  ur.linAdd(nu,dur); 
  uc += nu*duc;
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
  this->wksp->restore();
  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::restoreMpcStatus);
  if(ngrbms) rebuildGtGtilda();
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::solve(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u)
{
 if (verboseFlag) filePrint(stderr," ... Begin FETI-DP Solve            ...\n");

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
 }
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::solveCG(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u)
{
 t7 = -getTime(); this->times.solve -= getTime();

 // vectors
 GenDistrVector<Scalar> &fr      = this->wksp->ret_fr();  
 GenDistrVector<Scalar> &ur      = this->wksp->ret_ur(); 
 GenDistrVector<Scalar> &dur     = this->wksp->ret_du(); 
 GenDistrVector<Scalar> &lambda  = this->wksp->ret_lambda();      // lagrange multipliers
 GenDistrVector<Scalar> &r       = this->wksp->ret_r();           // residual
 GenDistrVector<Scalar> &w       = this->wksp->ret_w();           // projected residual
 GenDistrVector<Scalar> &y       = this->wksp->ret_y();           // re-projected residual
 GenDistrVector<Scalar> &z       = this->wksp->ret_z();           // preconditioned residual
 GenDistrVector<Scalar> &p       = this->wksp->ret_p();           // search direction
 GenDistrVector<Scalar> &Fp      = this->wksp->ret_Fp(); 
 GenDistrVector<Scalar> &deltaU  = this->wksp->ret_deltaU(); 
 GenDistrVector<Scalar> &fw      = this->wksp->ret_fw(); 
 GenVector<Scalar> &fc           = this->wksp->ret_fc(); 
 GenVector<Scalar> &uc           = this->wksp->ret_uc();  
 GenVector<Scalar> &duc          = this->wksp->ret_duc(); 

 double ww, ww0, ff, error;
 int iter = 0;
 Scalar pFp, nu;

 nMatVecProd = nRebuildGtG = nRebuildCCt = nLinesearch = nLinesearchIter = nSubIterDual = nSubIterPrimal = nStatChDual = nStatChPrimal = 0;
 if(globalFlagCtc && this->numSystems > 0) this->resetOrthoSet();

 // extract fr, fc and fw from f
 ff = extractForceVectors(f, fr, fc, fw);
 
 // Compute initial lagrange multipliers: lambda^0 = P * lambda^00 - G*(G^T*G)^-1 * e, where e = R^T*f ... also gamma = G^T*lambda+e
 // also for feti-dpc: expand active set and re-compute initial lagrange multipliers if lambda^0 is not feasible
 computeL0(lambda, f); 

 // Compute initial residual: r^0 = F*lambda^0 + d ... also uc = Kcc^-1(fc+Krr^-1*Krc^T*lambda), ur = Krr^-1*(fr-Br^T*lambda-Krc*uc)
 localSolveAndJump(fr, lambda, ur, fc, uc, r, fw); 

 // Compute initial projected residual: w^0 = P^T * r^0 
 // also for feti-dpc: contract active set if w^0 is not proportional (primal planing)
 ww0 = ww = tProject(r, w); 
 if(verboseFlag) filePrint(stderr," ... Initial residual norm %e\n", sqrt(ww0));

 //if(ww == 0.0) { 
 //  error = 0.0;
 //  this->times.iterations[this->numSystems].stagnated = 0;
 //} 
 //else { 
   // Multiple rhs prediction (note: not used for contact)
   if(this->predict(w, lambda)) {
     localSolveAndJump(fr, lambda, ur, fc, uc, r, fw);
     ww = tProject(r, w); 
     if(verboseFlag) filePrint(stderr," ... Initial residual norm after MRHS prediction %e\n", sqrt(ww));
   }

   if(verboseFlag) filePrint(stderr," Iteration  Relative Primal Error  Relative Dual Error\n");

   for(iter = 0; iter < this->maxiter; ++iter, ++iterTotal) {
     dbg_alloca(0);

     // Precondition: z = M^-1 * w
     error = preCondition(w, z);

     // Check stopping criteria
     bool stop = checkStoppingCriteria(iter, error, ff);

     // Print errors
     if(verboseFlag && (stop || ((this->fetiInfo->numPrint() > 0) && (iter % this->fetiInfo->numPrint() == 0))))
       filePrint(stderr," %4d %23.6e %21.6e\n", iter, (ff == 0.0) ? sqrt(error) : sqrt(error/ff), sqrt(ww/ww0));

     // Exit CG iterations loop if one of the stopping criteria has been satisfied
     if(stop) break;

     // Krylov acceleration
     if(this->fetiInfo->nlPrecFlg) this->nlPreCondition(w, z);

     // Re-project: y = P * z
     project(z, y); 

     // Search direction
     this->orthogonalize(y, p);
   
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
     // optional linesearch to adjust step length nu if sequence is not monotonic
     // also for feti-dpc: reduce step and expand active set if lambda is not feasible (dual planing)
     update(nu, lambda, p, r, Fp, ur, dur, uc, duc);

     // Project: w = P^T * r 
     // also for feti-dpc: contract active set if w is not proportional (primal planing)
     ww = tProject(r, w); 

     // add search direction to orthoset or reset if necessary
     if(globalFlagCtc && (dualStatusChange || primalStatusChange || stepLengthChange)) this->resetOrthoSet();
     else this->orthoAddCG(p, Fp, pFp);
   }

   ur += deltaU; // make solution compatible ur += deltaU
 //}

 // Assemble and store primal solution u
 mergeSolution(ur, uc, u, lambda);

 // Store number of iterations, errors, timings and memory used
 this->setAndStoreInfo(iter+1, (ff == 0.0) ? error : (error/ff), (ww0 == 0.0) ? ww : (ww/ww0));
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

 this->initGMRES(z);

 bool primalresidual = this->fetiInfo->gmresResidual;
 int J = 0;
 const int ReStep = this->fetiInfo->maxortho; // PJSA 1-23-08 (for restarted GMRES)

 for(int iter = 0; true; ++iter) {
   // Arnoldi iteration (Algorithm see Saad SISC) 
   for (int j=0; j<ReStep; j++, J++, ++iterTotal) {

     localSolveAndJump(z, dur, duc, Fp); // Fp = F*z

     error = preCondition(Fp, *medvec);   // medvec = M^-1*Fp

     // Do Arnoldi step 
     double resGMRES = this->orthoAddGMRES(z, *medvec);   

     if((fabs(resGMRES)<=sqrt(this->epsilon2*ff)) || (J == this->maxiter-1) || primalresidual) {

       primalresidual = true; // Since now we compute the primal residual in each step
         
       this->GMRESSolution(*lambda);

       localSolveAndJump(*lambda, dur, duc, Fp);  // Fp = F*lambda

       r.linC(1.0,*rzero,1.0,Fp); // r = r0 + Fp

       error = preCondition(r, *medvec); // medvec = M^-1*r

       bool isConverged = ((sqrt(error) < sqrt(this->epsilon2*ff)) || (J == this->maxiter-1));

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
     this->GMRESSolution(*lambda);  // compute incremental solution lambda
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
   this->initGMRES(z);
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
 GenDistrVector<Scalar> &work1   = this->wksp->ret_du();
 GenDistrVector<Scalar> &lambda  = this->wksp->ret_lambda();  // Lagrange multipliers
 GenDistrVector<Scalar> &r       = this->wksp->ret_r();       // residual
 GenDistrVector<Scalar> &z       = this->wksp->ret_z();       // preconditioned residual
 GenDistrVector<Scalar> &p       = this->wksp->ret_p();       // search direction
 GenDistrVector<Scalar> &Fp      = this->wksp->ret_Fp();
 GenDistrVector<Scalar> &Fz      = this->wksp->ret_y();
 GenDistrVector<Scalar> &fw      = this->wksp->ret_fw();
 GenDistrVector<Scalar> &deltaU  = this->wksp->ret_deltaU();
 GenVector<Scalar> &fc           = this->wksp->ret_fc();
 GenVector<Scalar> &uc           = this->wksp->ret_uc();
 GenVector<Scalar> &work2        = this->wksp->ret_duc();

 double rr, rr0, ff, error;
 int iter = 0;

 // extract fr, fc and fw from f
 ff = extractForceVectors(f, fr, fc, fw);

 // compute the initial guess for for the Lagrange Multipliers (zero if no projection)
 computeL0(lambda, f);

 // compute the initial residual
 // Solve uc = Kcc^-1(  fc + Krr^-1 Krc^T lambda)
 //       ur = Krr^-1 ( fr - Br^T lambda - Krc uc)
 //       r = Br ur
 localSolveAndJump(fr, lambda, ur, fc, uc, r, fw);
 rr0 = rr = r.sqNorm();
 if(verboseFlag) filePrint(stderr," ... Initial residual norm %e\n", sqrt(rr0));

 // multiple rhs prediction
 if(this->predictGCR(r, lambda)) {
   localSolveAndJump(fr, lambda, ur, fc, uc, r, fw); 
   rr = r.sqNorm();
   if(verboseFlag) filePrint(stderr," ... Initial residual norm after MRHS prediction %e\n", sqrt(rr));
 }

 if(verboseFlag) filePrint(stderr," Iteration  Relative Primal Error  Relative Dual Error\n");

 for(iter = 0; iter < this->maxiter; ++iter, ++iterTotal) {
   dbg_alloca(0);

   // Precondition: z = F^-1 r
   error = preCondition(r, z);

   // Check stopping criteria
   bool stop = checkStoppingCriteria(iter, error, ff);

   // Print errors
   if(verboseFlag && (stop || ((this->fetiInfo->numPrint() > 0) && (iter % this->fetiInfo->numPrint() == 0))))
     filePrint(stderr," %4d %23.6e %21.6e\n", iter, (ff == 0.0) ? sqrt(error) : sqrt(error/ff), sqrt(rr/rr0));

   // Exit CG iterations loop if one of the stopping criteria has been satisfied
   if(stop) break;

   localSolveAndJump(z, work1, work2, Fz);

   this->orthogonalizeGCR(z, Fz, p, Fp);  // computes new p, Fp

   Scalar FpFp = Fp * Fp;

   Scalar rFp = r * Fp;
   Scalar nu = -(rFp/FpFp);

   // Update solution: lambda += nu*p, r += nu*Fp
   lambda.linAdd(nu, p);
   r.linAdd(nu, Fp);
   rr = r.sqNorm();

   this->orthoAddGCR(p, Fp, FpFp);
 }

 // get primal solution 
 localSolveAndJump(fr, lambda, ur, fc, uc, r, fw);

 ur += deltaU;

 // merge ur and uc into u
 mergeSolution(ur, uc, u, lambda);

 // Store number of iterations, errors, timings and memory used
 this->setAndStoreInfo(iter+1, error/ff, sqrt(rr/rr0));
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
    if(domain->solInfo().isCoupled) this->distributeForce(fr, fw);
    else this->distributeForce(fr);
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

// RT 05/08/2010: bug in the g++ compiler
  if(ff == 0.0) {
     //filePrint(stderr, " *** WARNING: norm of rhs = 0 \n");
     return 1.0;
  }
  else return ff;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::printSummary(int iter)
{
  if(verboseFlag) {
    if(globalFlagCtc) {
      filePrint(stderr," -----------------------------------------------\n");
      filePrint(stderr," number of main iter.                = %5d ...\n", iter);
      filePrint(stderr," number of line searches             = %5d ...\n", nLinesearch);
      filePrint(stderr," number of line search sub-iter.     = %5d ...\n", nLinesearchIter);
      filePrint(stderr," number of dual planings             = %5d ...\n", nStatChDual); 
      filePrint(stderr," number of dual planing sub-iter.    = %5d ...\n", nSubIterDual);
      filePrint(stderr," number of primal planings           = %5d ...\n", nStatChPrimal);
      filePrint(stderr," number of primal planing sub-iter.  = %5d ...\n", nSubIterPrimal);
      filePrint(stderr," total number of matrix-vector op.s  = %5d ...\n", nMatVecProd); 
      filePrint(stderr," number of GtG rebuilds              = %5d ...\n", nRebuildGtG);
      filePrint(stderr," -----------------------------------------------\n");
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
}

template<class Scalar>
double
GenFetiDPSolver<Scalar>::preCondition(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Mv, bool errorFlag)
{
  // this function does dual mpc (CCt) preconditioning in addition to the usual feti preconditioning
  double error;
  if(mpcPrecon) {
    if(globalFlagCtc && (dualStatusChange || primalStatusChange) && this->fetiInfo->rebuildcct && CCtsolver) rebuildCCt();
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
  if(ngrbms) {
    if(!geometricRbms && this->fetiInfo->corners != FetiInfo::noCorners) {
      GenVector<Scalar> &alpha = this->wksp->ret_alpha();
      for(int i= 0; i<KccSolver->numRBM(); ++i) 
        for(int j=0; j<KccSolver->neqs(); ++j) 
          uc[j] += alpha[i]*kccrbms[i*KccSolver->neqs()+j];
      for(int i=0; i<this->nsub;++i) this->sd[i]->addTrbmRalpha(kccrbms, KccSolver->numRBM(), KccSolver->neqs(), alpha.data(), ur.subData(i));
    }
  }

  execParal4R(this->nsub, this, &GenFetiDPSolver<Scalar>::mergeUr, ur, uc, u, lambda);
  if(ngrbms) { // add rigid body modes
    if(geometricRbms || this->fetiInfo->corners == FetiInfo::noCorners) {
      GenVector<Scalar> &alpha = this->wksp->ret_alpha(); // amplitudes of the rigid body modes
      execParal2R(this->nsub, this, &GenFetiDPSolver<Scalar>::addRalpha, u, alpha);
    }
    if(this->fetiInfo->uproj) computeProjectedDisplacement(u);
    else filePrint(stderr, " ... Do not project the displacement ...\n"); //HB
  }
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

 nMatVecProd++;
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
 if(this->myCPU == 0 && (this->fetiInfo->outerloop == FetiInfo::CG)
    && (ScalarTypes::Real(pHFp) < 0.0 || fabs(ScalarTypes::Imag(pHFp)) > 1.0e-10))
   cerr << " *** WARNING: x^H F x = " << pHFp << ", must be positive and real for any x when F is Hermitian and positive definite. CG may not work \n";
#endif
 nMatVecProd++;
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

  paralApply(this->nsub, this->sd, &BaseSub::findEdgeNeighbors);  // PJSA

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
  if(K) {
    K->setPrintNullity(false);
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
  this->distributeForce(fr);
  GenVector<Scalar> &fc  = this->wksp->ret_fc();
  getFc(f, fc);
#ifdef DISTRIBUTED
  this->fetiCom->globalSum(fc.size(), fc.data());
#endif
  double mpcerr = 0.0;

  for(int i=0; i<this->nsub; ++i) mpcerr += this->sd[i]->getMpcError();
#ifdef DISTRIBUTED
  mpcerr = this->fetiCom->globalSum(mpcerr);
#endif
  //cerr << "mpc error = " << mpcerr << endl;

  return (fr.sqNorm() + fc.sqNorm() + mpcerr);
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::subtractMpcRhs(int iSub, GenDistrVector<Scalar> &dv1) 
{
  this->sd[iSub]->subtractMpcRhs(dv1.subData(this->sd[iSub]->localSubNum()));
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
  nRebuildGtG++;
  startTimerMemory(this->times.coarse1, this->times.memoryGtG);

  if(GtGtilda == NULL) {
    GtGtilda = this->newSolver(this->fetiInfo->auxCoarseSolver, coarseConnectGtG, eqNumsGtG, this->fetiInfo->grbm_tol, GtGsparse);
    GtGtilda->setPrintNullity(this->fetiInfo->contactPrintFlag && this->myCPU == 0);
  } else
  GtGtilda->zeroAll();
  execParal(nGroups1, this, &GenFetiDPSolver<Scalar>::assembleGtG);
#ifdef DISTRIBUTED
  GtGtilda->unify(this->fetiCom);
#endif
  startTimerMemory(this->times.pfactor, this->times.memoryGtGsky);
  GtGtilda->parallelFactor();
  //int* pivnull_list = ((GenMumpsSolver<Scalar>*) GtGtilda)->getPivnull_list();
  stopTimerMemory(this->times.pfactor, this->times.memoryGtGsky);

  stopTimerMemory(this->times.coarse1, this->times.memoryGtG);
}

template<class Scalar> 
void
GenFetiDPSolver<Scalar>::assembleGtG(int iGroup)
{
  // assembles groups in parallel, subdomains with same group sequentially 
  // threadsafe implementation - avoids simultaneous writing to same memory
  // note: the distributed version will work for shared memory too, but the 
  // alternative code is a bit more efficient
  int i;
#ifdef DISTRIBUTED
  for(i = 0; i < this->nsub; ++i) {
    if(this->sd[i]->group == groups[iGroup]) 
      this->sd[i]->assembleGtGsolver(GtGsparse);
  }
#else
  int *grsubs = (*groupToSub)[iGroup];
  for(i = 0; i < groupToSub->num(iGroup); ++i) {
    int iSub = grsubs[i];
    this->sd[iSub]->assembleGtGsolver(GtGsparse);
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
  if(ngrbmGr) { delete [] ngrbmGr; ngrbmGr = 0; }
  if(groups) { delete [] groups; groups = 0; }
  if(groupToSub != bodyToSub) { delete groupToSub; groupToSub = 0; }
  if(subToGroup && (subToGroup != subToBody)) { delete subToGroup; subToGroup = 0; }
  if(subToBody) { delete subToBody; subToBody = 0; }
  if(GtGtilda) { delete GtGtilda; GtGtilda = 0; }
  //if(glCrnGroup) { delete [] glCrnGroup; glCrnGroup = 0; }
  
  if(KccSparse) { delete KccSparse; KccSolver = 0; KccSparse = 0; }
  if(cornerToSub) { delete cornerToSub; cornerToSub = 0; }
  if(cornerEqs) { delete cornerEqs; cornerEqs = 0; }
  if(coarseConnectGtG) { delete coarseConnectGtG; coarseConnectGtG = 0; }
  if(eqNumsGtG) { delete eqNumsGtG; eqNumsGtG = 0; }
  
  if(CCtsolver) { delete CCtsolver; CCtsolver = 0; }
  if(subsWithMpcs)  { delete subsWithMpcs ; subsWithMpcs = 0; }
  if(mpcSubMap) { delete [] mpcSubMap; mpcSubMap = 0; }
  // don't delete: this->sd, this->subToSub, this->mpcToSub, mpcToMpc, this->cpuToSub, this->fetiInfo, bodyToSub 
  // interfPattern, this->glSubToLoc, mpcToCpu 
  if(mpcPat) { delete mpcPat; mpcPat = 0; }
  if(kccrbms) { delete [] kccrbms; kccrbms = 0; }
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::getLocalMpcForces(int iSub, double *mpcLambda)
{
  // mpcLambda is the local MPC forces for subdomain iSub (required for salinas)
  GenVector<Scalar> &uc = this->wksp->ret_uc();
  this->sd[iSub]->getLocalMpcForces(mpcLambda, cornerEqs, mpcOffset, uc);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::addMpcRHS(int iMPC, Scalar *fcstar)
{
 int dof = cornerEqs->firstdof(mpcOffset+iMPC);
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
  cornerToSub = 0; cornerEqs = 0;
  mpcOffset = 0;
  ngrbmGr = 0; nGroups = 0; nGroups1 = 0; groups = 0;
  groupToSub = 0; bodyToSub = 0; subToBody = 0; subToGroup = 0;
  GtGtilda = 0; coarseConnectGtG = 0; eqNumsGtG = 0;
  ngrbms = 0;
  mpcToMpc = 0; CCtsolver = 0; 
  subsWithMpcs = 0; numSubsWithMpcs = 0; mpcSubMap = 0;
  dualStatusChange = primalStatusChange = stepLengthChange = false;
  nMatVecProd = nRebuildGtG = nRebuildCCt = nLinesearch = nLinesearchIter = nSubIterDual = nSubIterPrimal = nStatChDual = nStatChPrimal = 0;
  mpcPat = 0;
  kccrbms = 0;
}

template<class Scalar>
int
GenFetiDPSolver<Scalar>::numRBM()
{
  bool useKccSolver = (this->glNumMpc == 0 && !geometricRbms);
  if(GtGtilda && !useKccSolver) {
    return GtGtilda->numRBM();
  }
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
GenFetiDPSolver<Scalar>::reconstruct()
{
  // 1. reset the orthoset
  this->resetOrthoSet();
 
  if(domain->solInfo().isCoupled) domain->computeCoupledScaleFactors();

  if(this->fetiInfo->isEdgeAugmentationOn() && this->fetiInfo->numdir > 0) { // augmentation depends on freq/k
    // 2. reconstruct local Kcc etc. since size of Kcc may have changed
    startTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);
    paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::constructKcc);
    if(domain->solInfo().isCoupled)  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::constructKcw);
    stopTimerMemory(this->times.constructMatrices, this->times.memorySubMatrices);

    // 3. rebuild augmentation Q
    if(verboseFlag) filePrint(stderr," ... Rebuild Edge Augmentation (Q)  ... \n");
    if(this->fetiInfo->waveMethod != FetiInfo::averageMat) computeLocalWaveNumbers();
    paralApplyToAll(this->nsub, this->sd, &BaseSub::zeroEdgeDofSize);
    paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::makeQ);  // rebuild augmentation matrix
  }

  geometricRbms = 0;
  ngrbms = 0;
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

 // MPCs 
 mpcPrecon = false;
 if(this->glNumMpc > 0) {
   if(this->fetiInfo->mpc_scaling == FetiInfo::kscaling) { // MPC stiffness scaling
     FSCommPattern<Scalar> *mpcDiagPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU,
                                                                   FSCommPattern<Scalar>::CopyOnSend);
     for(int iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setMpcDiagCommSize(mpcDiagPat);
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
     deleteCCt();
     buildCCt();
   }
 }

  // 5. factor local matrices
  if(verboseFlag) filePrint(stderr," ... Factor Subdomain Matrices      ... \n");
  startTimerMemory(this->times.factor, this->times.memoryFactor);
  if(this->fetiInfo->solvertype != FetiInfo::spooles && this->fetiInfo->solvertype != FetiInfo::mumps)
    timedParal(this->times.factorMat, this->nsub, this, &GenFetiDPSolver<Scalar>::factorLocalMatrices);
  else // PJSA: spooles and mumps are not thread safe
    for(int iSub=0; iSub<this->nsub; ++iSub) factorLocalMatrices(iSub);
  stopTimerMemory(this->times.factor, this->times.memoryFactor);

  if(verboseFlag) filePrint(stderr, " ... Reconstruct Kcc solver         ... \n");
  delete KccSparse; KccSparse = 0; KccSolver = 0;
  makeKcc();

  if(GtGtilda) { delete GtGtilda; GtGtilda = 0; } 
  if(ngrbms) makeGtG();

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

  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::initMpcStatus);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::projectActiveIneq(GenDistrVector<Scalar> &x, GenDistrVector<Scalar> &y)
{
  if(&x != &y) y = x;
  execParal1R(this->nsub, this, &GenFetiDPSolver<Scalar>::subProjectActiveIneq, y);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subProjectActiveIneq(int iSub, GenDistrVector<Scalar> &y)
{
  this->sd[iSub]->projectActiveIneq(y.subData(this->sd[iSub]->localSubNum()));
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::multG(GenVector<Scalar> &x, GenDistrVector<Scalar> &y, double alpha, double beta)
{
  // y = alpha*G*x + beta*y
  if(beta == 0) y.zero(); else y *= beta;
  execParal3R(this->nsub, this, &GenFetiDPSolver<Scalar>::subMultG, x, y, alpha); // y += alpha*G*x
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subMultG(int iSub, GenVector<Scalar> &x, GenDistrVector<Scalar> &y, double alpha)
{
  this->sd[iSub]->multG(x, y.subData(this->sd[iSub]->localSubNum()), alpha); 
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::trMultG(GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha, double beta)
{
  // y = alpha*G^T*x + beta*y
  if(beta == 0) y.zero(); else y *= beta;
  GenVector<Scalar> v(ngrbms, 0.0);
  execParal3R(nGroups1, this, &GenFetiDPSolver<Scalar>::subTrMultG, x, v, alpha); // v += alpha*G^T*x
#ifdef DISTRIBUTED
  this->fetiCom->globalSum(ngrbms, v.data());
#endif
  y += v;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::subTrMultG(int iGroup, GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha)
{
#ifdef DISTRIBUTED
  for(int i=0; i<this->nsub; ++i) {
    if(this->sd[i]->group == groups[iGroup])
      this->sd[i]->trMultG(x.subData(this->sd[i]->localSubNum()), y, alpha); 
  }
#else
  int *grsubs = (*groupToSub)[iGroup];
  for(int i = 0; i < groupToSub->num(iGroup); ++i) {
    int iSub = grsubs[i];
    this->sd[iSub]->trMultG(x.subData(this->sd[iSub]->localSubNum()), y, alpha); 
  }
#endif
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::project(GenDistrVector<Scalar> &z, GenDistrVector<Scalar> &y, int eflag)
{
  startTimerMemory(this->times.project, this->times.memoryProject1);
  // if eflag is true this function computes y as the projection of z on to the feasible subset in which all active constraints remain active
  // if eflag is false this function computes y as the projection of z on to the tangent subspace

  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::saveMpcStatus1);

  GenDistrVector<Scalar> x(z);
  bool status_change = false;

  int i;
  for(i = 0; i < this->fetiInfo->dual_plan_maxit+1; ++i) {
    // y = P_i*x
    projectActiveIneq(x, y);

    // res = -e-G^T*x
    GenVector<Scalar> res(ngrbms);
    if(ngrbms) { 
      trMultG(y, res, -1.0, 0.0);
      if(eflag) { 
        GenVector<Scalar> &e = this->wksp->ret_e();
        res -= e; 
      }
    }

    // check stopping criteria
    if(i > 0) {
      double resnorm = (eflag && ngrbms) ? res.norm() : 0;
      if(this->fetiInfo->contactPrintFlag && this->myCPU == 0) cerr << "dual planing: iteration = " << i << ", residual = " << resnorm << endl;
      if(/*resnorm <= this->fetiInfo->dual_proj_tol ||*/ !status_change) break;
      else if(i == MAX(1,this->fetiInfo->dual_plan_maxit)) {
        if(this->myCPU == 0) cerr << "warning: dual planing did not converge after " << i << " iterations. Error = " << resnorm << endl;
        // note: if we break the loop here then y will not be feasible wrt the equality constraints (i.e. G^T*y != e)
        // if we don't break here then y will not be feasible wrt the inequality constraints
        break;
      }
    }
  
    // res = (G^T*P_i*G)^(-1)*res
    if(GtGtilda) GtGtilda->reSolve(res);

    // x = x + G*res
    if(ngrbms) multG(res, x, 1.0, 1.0);

    // update the active set
    if(eflag && globalFlagCtc)
      status_change = updateActiveSet(x, 0, -this->fetiInfo->dual_plan_tol);
  }

  if(dualStatusChange = (i > 1)) {
    nSubIterDual += (i-1);
    nStatChDual++;
    if(this->fetiInfo->contactPrintFlag && this->myCPU == 0) cerr << endl;
  }

  stopTimerMemory(this->times.project, this->times.memoryProject1);
}

template<class Scalar>
double
GenFetiDPSolver<Scalar>::tProject(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &w)
{
  startTimerMemory(this->times.project, this->times.memoryProject1);
  // this function computes w as the projection of r on to the tangent cone at a feasible point
  // unless the "chopped gradient" error term is proportional, in which case this function
  // computes w as the projection of r on to the tangent subspace
  // note: this function also computes alpha which is stored in the feti workspace

  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::saveMpcStatus1);

  GenDistrVector<Scalar> x(r);
  GenVector<Scalar> &alpha = this->wksp->ret_alpha();
  GenDistrVector<Scalar> &gc = this->wksp->ret_gc(); // chopped gradient
  GenDistrVector<Scalar> &gf = this->wksp->ret_gf(); // free gradient
  alpha.zero();
  bool status_change = false;
  bool proportional = false;

  int i;
  for(i = 0; i < this->fetiInfo->primal_plan_maxit+1; ++i) {
    // w = P_i*x
    projectActiveIneq(x, w);

    // res = -G^T*x
    GenVector<Scalar> res(ngrbms);
    if(ngrbms) trMultG(w, res, -1.0, 0.0);

    // check stopping criteria
    if(i > 0) {
      double resnorm = (ngrbms) ? res.norm() : 0;
      if(this->fetiInfo->contactPrintFlag && this->myCPU == 0) cerr << "primal planing: iteration " << i << ", residual = " << resnorm << endl;
      if(/*resnorm <= this->fetiInfo->primal_proj_tol ||*/ !status_change) break;
      else if(i == MAX(1,this->fetiInfo->primal_plan_maxit)) {
        if(this->myCPU == 0) cerr << "warning: primal planing did not converge after " << i << " iterations. " << endl;
        break;
      }
    }
  
    // res = (G^T*P_i*G)^(-1)*res
    if(GtGtilda) GtGtilda->reSolve(res);

    // x += G*res, alpha += res
    if(ngrbms) multG(res, x, 1.0, 1.0);
    alpha += res;

    // update active set unless error is proportional
    if(globalFlagCtc) {
      execParal3R(this->nsub, this, &GenFetiDPSolver<Scalar>::split, x, gf, gc);
      proportional = (i == 0 && (gc.norm() <= this->fetiInfo->gamma*gf.norm()));
      if(!proportional)
        status_change = updateActiveSet(x, 1, -this->fetiInfo->primal_plan_tol);
    }
  }

  if(primalStatusChange = (i > 1)) {
    nSubIterPrimal += (i-1);
    nStatChPrimal++;
    if(this->fetiInfo->contactPrintFlag && this->myCPU == 0) cerr << endl;
  }

  stopTimerMemory(this->times.project, this->times.memoryProject1);

  return w.sqNorm() + (proportional ? gc.sqNorm() : 0);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::computeL0(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f)
{
  lambda.zero(); // XXXX consider starting from previous iteration or time (load) step for nonlinear dynamics (statics)
                 // currently lambda is initialized to zero in workspace constructor for feti-dp and a 
                 // new workspace is constructed every time the solver is refactored
  if(this->glNumMpc) {
    if(ngrbms > 0) makeE(f); // compute e = R^T*f
    project(lambda, lambda, true); 
  }
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

template<class Scalar>
void
GenFetiDPSolver<Scalar>::reconstructMPCs(Connectivity *_mpcToSub, Connectivity *_mpcToMpc, Connectivity *_mpcToCpu)
{
 this->mpcToSub   = _mpcToSub;    // MPC to subdomain connectivity
 this->glNumMpc = (this->mpcToSub) ? this->mpcToSub->csize() : 0;
 mpcToMpc   = _mpcToMpc;    // MPC to MPC connectivity (PJSA: used for CC^t preconditioner
 mpcToCpu   = _mpcToCpu;
 globalFlagCtc = domain->getNumCTC();
#ifdef DISTRIBUTED
 globalFlagCtc = this->fetiCom->globalMax((int) globalFlagCtc);
#endif
 int iSub;

 // create this->vPat FSCommPattern object, used to send/receive a scalar vector (interfaceDOFs)
 delete this->vPat;
 this->vPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<Scalar>::CopyOnSend);
 for(iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setDofCommSize(this->vPat);
 this->vPat->finalize();
 
 // create this->sPat FSCommPattern objects, used to send/receive a single integer
 delete this->sPat;
 this->sPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
 for(iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setCommSize(this->sPat, 1);
 this->sPat->finalize();

 if(mpcPat) { delete mpcPat; mpcPat = 0; }
 if(globalFlagCtc) {
   mpcPat = new FSCommPattern<int>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<int>::CopyOnSend);
   for(iSub=0; iSub<this->nsub; ++iSub) this->sd[iSub]->setMpcCommSize(mpcPat);
   mpcPat->finalize();
 }
 
 int tInterfLen    = 0;
 this->halfSize = 0;
 for(iSub = 0; iSub < this->nsub; ++iSub) {
   this->interface.domLen[iSub] = this->sd[iSub]->interfLen();

   this->sd[iSub]->computeMasterFlag(this->mpcToSub);
   this->fetiOps[iSub]->setHalfOffset(this->halfSize);
   this->halfSize      += this->sd[iSub]->halfInterfLen();
   tInterfLen    += this->interface.domLen[iSub];
 }
 this->interface.len = tInterfLen;

 // PJSA: compute the masterFlags
 bool *interfaceMasterFlag = new bool[tInterfLen];
 this->interface.recomputeOffsets();
 for(iSub = 0; iSub < this->nsub; ++iSub) {
   bool *subMasterFlag = this->sd[iSub]->getMasterFlag();
   int subOffset = this->interface.subOffset[iSub];
   int j;
   for(j=0; j<this->interface.domLen[iSub]; ++j)
     interfaceMasterFlag[subOffset+j] = subMasterFlag[j];
 }
 this->interface.setMasterFlag(interfaceMasterFlag);
 // don't delete interfaceMasterFlag

 // Allocate space for reorthogonalization set
 this->times.memoryOSet -= memoryUsed();
 if(this->fetiInfo->outerloop == 0) {
   delete this->oSetCG;
   this->oSetCG = (this->fetiInfo->maxortho > 0) ? new GenCGOrthoSet<Scalar>(this->halfSize, this->fetiInfo->maxortho, this->fetiCom) : 0;
 }
 else if(this->fetiInfo->outerloop == 1) {
   delete this->oSetGMRES;
   this->oSetGMRES = new GenGMRESOrthoSet<Scalar>(this->halfSize, this->fetiInfo->maxortho, this->fetiCom);
 }
 else {
   delete this->oSetGCR;
   this->oSetGCR = new GenGCROrthoSet<Scalar>(this->halfSize, this->fetiInfo->maxortho, this->fetiCom);
 }
 this->times.memoryOSet += memoryUsed();

 paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::rebuildKbb);

 paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::initScaling);
 if((this->fetiInfo->scaling == FetiInfo::kscaling) || ((this->fetiInfo->mpc_scaling == FetiInfo::kscaling) && (this->glNumMpc_primal > 0))
    || (this->fetiInfo->augment == FetiInfo::WeightedEdges)) {
   execParal(this->nsub, this, &GenFetiSolver<Scalar>::sendScale);
   this->vPat->exchange();
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::collectScaling, this->vPat);
 }

 deleteCCt(); mpcPrecon = false;
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
/* this is called in refactor
     buildCCt();
*/
   }
   if(subsWithMpcs) { delete subsWithMpcs; subsWithMpcs = 0; }
   if(mpcSubMap) { delete [] mpcSubMap; mpcSubMap = 0; }
 }

 paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::cleanMpcData);

/*
 if(cornerEqs->size() > 0) 
   paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::reMultKcc);

 //paralApplyToAll(this->nsub, this->sd, &GenSubDomain<Scalar>::cleanMpcData);

 if(ngrbms) makeGtG();  // currently G = C^T*R (ie restriction of R to mpc interface)

 delete this->wksp;
 this->times.memoryDV -= memoryUsed();
 int numC = (KccSolver) ? KccSolver->neqs() : 0;
 this->wksp = new GenFetiWorkSpace<Scalar>(this->interface, internalR, internalWI, ngrbms, numC, globalFlagCtc);
 this->times.memoryDV += memoryUsed();
*/
}

template<class Scalar>
bool
GenFetiDPSolver<Scalar>::checkStoppingCriteria(int iter, double error, double ff)
{
  // 1. check if maximum number of iterations reached
  if(iter == this->maxiter) {
    this->times.iterations[this->numSystems].stagnated = 0;
    return true;
  }
 
  // 2. check for convergence
  if(sqrt(error) <= std::max(this->fetiInfo->tol*sqrt(ff), this->fetiInfo->absolute_tol)) {
    this->times.iterations[this->numSystems].stagnated = 0;
    return true;
  }

  // 3. check for stagnation
  if(iter > 0 && (std::fabs(sqrt(error)-sqrt(lastError)) < std::max(this->fetiInfo->stagnation_tol*sqrt(lastError), this->fetiInfo->absolute_stagnation_tol))) {
     this->times.iterations[this->numSystems].stagnated = 1;
     return true;
  }

  lastError = error;
  return false;
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::zeroG()
{
  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::zeroG);
}

template<class Scalar>
void
GenFetiDPSolver<Scalar>::deleteG()
{
  paralApply(this->nsub, this->sd, &GenSubDomain<Scalar>::deleteG);
}

