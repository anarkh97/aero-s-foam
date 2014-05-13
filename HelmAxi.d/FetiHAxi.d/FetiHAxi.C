#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/SubDomain.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Math.d/BLKSparseMatrix.h>
#include <HelmAxi.d/AxiHElem.h>
#include <HelmAxi.d/FourierHelmBCs.h>
#include <Feti.d/FetiInfo.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <HelmAxi.d/MDAxiData.h>
#include <Threads.d/PHelper.h>
#include <Timers.d/Timing.h>
#include <Timers.d/GetTime.h>
#include <Feti.d/DistrVector.h>
#include <HelmAxi.d/FetiHAxi.d/DistrComplexVector.h>
#include <HelmAxi.d/FetiHAxi.d/GCRC.h>
#include <HelmAxi.d/FetiHAxi.d/FetiHAxi.h>
#include <Utils.d/Memory.h>


FetiHAxiSolver::FetiHAxiSolver(int _nsub, MDAxiData **_md, 
                Connectivity *_subToSub, int _MPCsize, FetiInfo *_fetiInfo) :
                internal(_nsub), interface(_nsub), 
                times(threadManager->numThr(), _nsub)
{

 // Compute memory used by FETI Solver
 times.memoryFETI -= threadManager->memoryUsed();

 nsub    = _nsub;       // Number of subdomains
 mdAxi   = _md;         // pointer to Array of all Subdomains
 fetiInfo = _fetiInfo;  // Feti solver information
 subToSub = _subToSub;  // subdomain to subdomain connectivity

 // Initialize pointers and variables
 subToEdge      = 0;    // subdomain to edge connectivity
 edgeToSub      = 0;    // edge to subdomain connectivity
 edgeToEdge     = 0;    // edge to edge connectivity 
 
 // Define FETI tolerance and maximum number of iterations
 double fetiTolerance = fetiInfo->tol;
 epsilon2 = fetiTolerance*fetiTolerance;
 maxiter  = fetiInfo->maxit;

 oSet = 0;

 // Initialize MPC data
 MPCSize = _MPCsize;
 CoarseMPC = 0;
 MPCSolver = 0;

 // Store the number of Fourier modes 
 numModes = mdAxi[0]->locBCs->numModes;

 computeLengths(halfSize, interface, internal);

 // Copy the sign of the subdomains
 int *SubSign = (int *) dbg_alloca(sizeof(int)*nsub);
 int i;
 for (i=0; i<nsub; ++i)
     SubSign[i] = mdAxi[i]->InterfSign;

 // Print some information of the local solvers
 switch (fetiInfo->solvertype) {
   default:
   case 0:
     fprintf(stderr," ... Complex Skyline local solvers  ... \n");
     break;
   case 1:
     fprintf(stderr," ... Complex Sparse local solvers   ... \n");
     break;
 }

 // Switch between the parallelization on Fourier modes for 1 subdomain
 // and the parallelization on the subdomains
 int nP = (nsub>numModes) ? nsub : numModes+1;

 // Allocate Subdomain matrices
 startTimerMemory(times.constructMatrices, times.memorySubMatrices);
 int **glBoundMap = (int **) dbg_alloca(sizeof(int*)*nsub);
 int **glInternalMap = (int **) dbg_alloca(sizeof(int*)*nsub);
 execParal(nsub, this, &FetiHAxiSolver::defineMatrices, glBoundMap, 
           glInternalMap);

 // Construct and Assemble Sub-domain matrices i.e. K, Kii, Kib, Kbb, Kuc
 fprintf(stderr," ... Assemble Subdomain Matrices    ... \n");
 timedParal(times.assembleMat, nP, this, 
            &FetiHAxiSolver::constructAndAssembleMat, SubSign, glBoundMap,
            glInternalMap);
 execParal(nsub,this,&FetiHAxiSolver::deleteGlMap,glBoundMap,glInternalMap);
 stopTimerMemory(times.constructMatrices, times.memorySubMatrices);

 // Factor matrices: K and Kii (if dirichlet preconditioner))
 fprintf(stderr," ... Factor Subdomain Matrices      ... \n");
 timedParal(times.factorMat, nP, this, &FetiHAxiSolver::factorMatrices);

 // Allocate interface buffer
 execParal(nsub, this, &FetiHAxiSolver::allocateBuffer);

 makeCoarse();

 times.memoryFETI += threadManager->memoryUsed();

}


void
FetiHAxiSolver::computeLengths(int &halfsize, 
                DistrInfo &interface, DistrInfo &internal) {

 // Compute total interface length, total internal length
 // and total half interface length
 int tInterfLen    = 0;
 int tLocalLen     = 0;
 int halfInterfLen = 0;
 int iSub;
 int totalFourier  = 2*numModes+1;

 for(iSub = 0; iSub < nsub; ++iSub) {
   interface.domLen[iSub] = totalFourier*mdAxi[iSub]->interfLen();
   internal.domLen[iSub]  = totalFourier*mdAxi[iSub]->localLen();

   mdAxi[iSub]->halfOffset = halfInterfLen;
   halfInterfLen += totalFourier*mdAxi[iSub]->halfInterfLen();
   tInterfLen    += totalFourier*mdAxi[iSub]->interfLen();
   tLocalLen     += totalFourier*mdAxi[iSub]->localLen();
 }

 interface.len = tInterfLen;
 internal.len  = tLocalLen;
 halfSize      = halfInterfLen;

}


void
FetiHAxiSolver::defineMatrices(int iS, int **glBoundMap, int **glInternalMap) {

 mdAxi[iS]->defineKs(glBoundMap, glInternalMap);

}


void 
FetiHAxiSolver::constructAndAssembleMat(int iP, int *SubSign, 
                int **glBoundMap, int **glInternalMap) {

 if (nsub>numModes) {
   mdAxi[iP]->makeKs(glBoundMap[iP], glInternalMap[iP]);
   mdAxi[iP]->Assemble();
   mdAxi[iP]->addInterface(SubSign);
 } 
 else 
   for (int iSub=0; iSub<nsub; ++iSub) {
     mdAxi[iSub]->makeKs(glBoundMap[iSub], glInternalMap[iSub], iP, iP+1);
     mdAxi[iSub]->Assemble(iP, iP+1);
     mdAxi[iSub]->addInterface(SubSign, iP, iP+1);
   }

}


void
FetiHAxiSolver::deleteGlMap(int iSub, int **glBoundMap, int **glInternalMap) {

 mdAxi[iSub]->deleteGlMap(glBoundMap, glInternalMap);

}


void
FetiHAxiSolver::allocateBuffer(int iSub) {

  mdAxi[iSub]->interfBuff = new DComplex[mdAxi[iSub]->interfLen()];

}


void
FetiHAxiSolver::factorMatrices(int iP) {

 if (nsub>numModes) {
   mdAxi[iP]->factorKC();
   mdAxi[iP]->factorKiiC();
 } 
 else 
   for (int iSub=0; iSub<nsub; ++iSub) {
     mdAxi[iSub]->factorKC(iP, iP+1);
     mdAxi[iSub]->factorKiiC(iP, iP+1);
   }

}


void
FetiHAxiSolver::countEdges(int iSub, int *edges) {

  // Get my global subdomain number
  int myNum = mdAxi[iSub]->localSubNum();

  // Count a subdomain's edges
  edges[myNum] = 0;
  int j;
  SComm *sc = mdAxi[iSub]->getSComm();
  int numNeighbors = sc->numNeighb;
  for (j=0; j<numNeighbors; ++j) {
      int neighbNum = (sc->subNums)[j];
      if (myNum<neighbNum)
          edges[myNum] += 1;
  }

}


void
FetiHAxiSolver::numberEdges(int iSub, int *eP, int* edges) {

 // first compute each subdomain's starting edge number
 int myNum = mdAxi[iSub]->localSubNum();
 int startI = eP[myNum];
 int fP = subToSub->offset(myNum)-myNum;
 eP[myNum] = fP;
 
 // number all subdomain's edges
 int jSub;
 SComm *sc = mdAxi[iSub]->getSComm();
 int numNeighbors = sc->numNeighb;
 for (jSub=0; jSub<numNeighbors; ++jSub) {
    int neighbNum = (sc->subNums)[jSub];
    if (myNum<neighbNum) {
       edges[fP] = startI;
       startI++;
       }
    else 
       edges[fP] = -1;
    // Send the numbered edges to all neighbors 
    // to complete the edge vector
    mdAxi[iSub]->setSendData(jSub, edges+fP);
    fP++;
 }

}


void
FetiHAxiSolver::receiveNeighbEdgeNums(int iSub, int *eP, int* edges) {

 // first compute each subdomain's starting edge number
 int myNum = mdAxi[iSub]->localSubNum();
 // number all subdomain's edges
 int jSub;
 SComm *sc = mdAxi[iSub]->getSComm();
 int numNeighbors = sc->numNeighb;
 for (jSub=0; jSub<numNeighbors; ++jSub) {
    int *en = (int *) sc->getExchangeData(jSub);
    if (*en >= 0) edges[eP[myNum]+jSub] = *en;
 }

}


void 
FetiHAxiSolver::makeEdgeConnectivity() {

 // First count number of edges per subdomain
 int *cx = new int[nsub+1];

 int i;
 for(i=0; i<nsub; ++i)
   cx[i] = 0;

 // count edges in parallel
 execParal(nsub, this, &FetiHAxiSolver::countEdges, cx);

 // Sum each subdomain's edge count to compute total number of edges
 // We modify 'edges' at the same time so that it has the first edge number
 // for each subdomain
 int numEdges = 0;
 for(i=0; i<nsub; ++i) {
   int tmp = numEdges;
   numEdges += cx[i];
   cx[i] = tmp;
 }

 // Find total number of edges including duplicated ones.
 int totalNumEdges = 0;
 for (i=0; i<nsub; ++i)
    totalNumEdges += mdAxi[i]->getSComm()->numNeighb; 

 int *connect = new int[totalNumEdges];
 cx[nsub] = totalNumEdges;

 execParal(nsub, this, &FetiHAxiSolver::numberEdges, cx, connect);

 execParal(nsub, this, &FetiHAxiSolver::receiveNeighbEdgeNums, cx, connect);

 // construct subdomain to edge connectivity
 subToEdge = new Connectivity(nsub, cx, connect);

 // create the edge to subdomain connectivity
 edgeToSub = subToEdge->reverse();

 // create the edge to edge connectivity
 edgeToEdge = edgeToSub->transcon(subToEdge);

}
 

void
FetiHAxiSolver::makeCoarse() {

 startTimerMemory(times.coarse1, times.memoryGtG);

 int totalFourier;

 fprintf(stderr," ... Make Coarse Matrices           ... \n");

 Connectivity *coarseConnectivity = 0;
 
 makeEdgeConnectivity();
 coarseConnectivity = edgeToEdge;

 // Renumber the coarse problem equations
 // 0 - no renumbering
 // 1 = sloan renumbering
 // 2 = rcm renumbering

 compStruct renumber = coarseConnectivity->renumByComponent(1);

 coarseEqs = new DofSetArray(coarseConnectivity->csize(), renumber.renum);

 int i,j;
 for (i=0; i<nsub; ++i) {
    int Neighbors = mdAxi[i]->getSComm()->numNeighb;
    for (j=0; j<Neighbors; ++j)
       coarseEqs->setWeight((*subToEdge)[i][j], mdAxi[i]->edgeDofSize[j]);
 }

 coarseEqs->finish();

 CoarseSolver = 0;

 switch (fetiInfo->nonLocalQ) {
   default:
   case 0:
     totalFourier = numModes+1;
     break;
   case 1:
     totalFourier = 2*numModes+1;
     break;
 }

 CoarseSolver = new ComplexSolver *[totalFourier];
 CoarseVector = new ComplexVector *[2*numModes+1];

 if (MPCSize)
   CoarseMPC = new FullMC *[2*numModes+1];
 
 execParal(totalFourier, this, &FetiHAxiSolver::allocateCoarseSolver,
           coarseConnectivity);

 execParal(2*numModes+1, this, &FetiHAxiSolver::allocateCoarseVectors);

 times.numCRNs = coarseEqs->size();

 execParal(nsub, this, &FetiHAxiSolver::makeCoarseGridDofs);
 timedParal(times.consMatrix, totalFourier, this, 
            &FetiHAxiSolver::assembleCoarse);

 stopTimerMemory(times.coarse1, times.memoryGtG);

 fprintf(stderr," ... Done Coarse Matrices           ... \n");

 // Factor Coarse matrices
 fprintf(stderr," ... Factor Coarse Matrices         ... \n");
 startTimerMemory(times.pfactor, times.memoryFactor);
 execParal(totalFourier, this, &FetiHAxiSolver::factorCoarseMatrices);
 stopTimerMemory(times.pfactor, times.memoryFactor);

 if (CoarseSolver[0]) {
   times.numRBMs = CoarseSolver[0]->numRBM();
   for (int Fourier=1; Fourier<totalFourier; ++Fourier)
     if (CoarseSolver[Fourier]->numRBM() > times.numRBMs)
        times.numRBMs = CoarseSolver[Fourier]->numRBM();
 }
 
}


void
FetiHAxiSolver::allocateCoarseSolver(int Fourier, 
                Connectivity *coarseConnectivity) {

 double tolerance = fetiInfo->grbm_tol;

 // CoarseSolver zero valued in the constructor
 if (coarseEqs->size()==0)
   CoarseSolver[Fourier] = 0; 
 else {
   switch (fetiInfo->gtgSolver) {
     default:
     case 0:
       CoarseSolver[Fourier] = new SkyMatrixC(coarseConnectivity,coarseEqs,
                               tolerance);
       break;
     case 1:
       CoarseSolver[Fourier] = new BLKSparseMatrixC(coarseConnectivity,
                               coarseEqs, tolerance);
       break;
   }
 }
}


void
FetiHAxiSolver::allocateCoarseVectors(int Fourier) {

 CoarseVector[Fourier]=new ComplexVector(coarseEqs->size(),DComplex(0.0,0.0));

 if (MPCSize) {

   //
   // Rows    of CoarseMPC[.] = coarseEqs->size()
   // Columns of CoarseMPC[.] = MPCSize
   //

   CoarseMPC[Fourier] = new FullMC(coarseEqs->size(),MPCSize); 
   (*CoarseMPC[Fourier]) = DComplex(0.0,0.0);

 }

}


void
FetiHAxiSolver::makeCoarseGridDofs(int iSub) {

 mdAxi[iSub]->makeCoarseGridDofs(coarseEqs, subToEdge);

}


void
FetiHAxiSolver::assembleCoarse(int Fourier) {

 int iSub;

 ComplexSparseMatrix *spm;

 if (CoarseSolver[Fourier]==0) {
   spm = 0;
 }
 else {
   switch (fetiInfo->gtgSolver) {
     default :
     case 0 :
       spm = (SkyMatrixC*) CoarseSolver[Fourier];
       break;
     case 1 :
       spm = (BLKSparseMatrixC*) CoarseSolver[Fourier];
       break;
   }  
 }

 for (iSub=0; iSub<nsub; ++iSub) {
   mdAxi[iSub]->assembleCoarse(Fourier, spm, CoarseMPC);
 }

 if (Fourier%10==0) {
   fprintf(stderr," ... Make Coarse : Fourier %d\n", Fourier);
 }

}


void
FetiHAxiSolver::factorCoarseMatrices(int Fourier) {

 if (CoarseSolver[Fourier])
   CoarseSolver[Fourier]->factor();

}


void
FetiHAxiSolver::makeStaticLoad(DistrComplexVector &f) {

 fprintf(stderr," ... Building the Force             ...\n");
 f.zero();
 int nP = (nsub>numModes) ? nsub : 2*numModes+1; 
 timedParal1R(times.buildRhs,nP,this,
               &FetiHAxiSolver::makeLocalStaticLoad,f);

}


void
FetiHAxiSolver::makeLocalStaticLoad(int iP, DistrComplexVector &f) {

 if (nsub>numModes) {
   mdAxi[iP]->buildRHS(f.subData(mdAxi[iP]->localSubNum()));
 }
 else {
   for (int iSub=0; iSub<nsub; ++iSub)
     mdAxi[iSub]->buildRHS(f.subData(mdAxi[iSub]->localSubNum()), iP, iP+1);
 }

}


void 
FetiHAxiSolver::solve(DistrComplexVector &f, DistrComplexVector &u) {

 fprintf(stderr," ... Complex GCR solver selected    ...\n");

 switch (fetiInfo->precno) {
   case (FetiInfo::noPrec) : 
      fprintf(stderr," ... with no preconditioner         ...\n");
      break;
   case (FetiInfo::lumped) : 
      fprintf(stderr," ... with lumped preconditioner     ...\n");
      break;
   case (FetiInfo::dirichlet) : {
      fprintf(stderr,"\n\n");
      fprintf(stderr," ... DIRICHLET NOT AVAILABLE YET    ...\n");
      fprintf(stderr," ... USE LUMPED PRECONDITIONER      ...\n");
      fprintf(stderr,"\n\n");
      }
   case (FetiInfo::identity) : {
      fprintf(stderr,"\n\n");
      fprintf(stderr," ... IDENTITY NOT AVAILABLE YET    ...\n");
      fprintf(stderr," ... USE LUMPED PRECONDITIONER      ...\n");
      fprintf(stderr,"\n\n");
      }
      break; 
 }
   
 fprintf(stderr," --------------------------------------\n");
 fprintf(stderr," |     Iterations loop monitoring     |\n");
 fprintf(stderr," --------------------------------------\n");

 // Set timers to evaluate the loop
 double temps;
 double tLoopMin = 1.0e6;
 double tLoopMax = 0.0;
 double tLoopAvg = 0.0;

 startTimerMemory(times.solve, times.memorySolve);

 times.memoryDV -= threadManager->memoryUsed();
 DistrComplexVector r(interface);
 DistrComplexVector z(interface);
 DistrComplexVector p(interface); 
 DistrComplexVector Fz(interface);
 DistrComplexVector Fp(interface);
 DistrComplexVector lambda(interface);
 DistrComplexVector ft(internal);
 deltaU = new DistrComplexVector(internal);
 deltaU->zero();
 deltaF = new DistrComplexVector(internal);
 deltaF->zero();
 times.memoryDV += threadManager->memoryUsed();

 ft = f;

 times.memoryOSet -= threadManager->memoryUsed();
 oSet = new GCRC(halfSize, fetiInfo->maxortho);
 times.memoryOSet += threadManager->memoryUsed();

 // Compute the square of a pseudo-norm of fVec
 double pseudoFNormSq = real(f^f);

 // Set some initial vectors at zero
 lambda.zero();

 // Compute r0 = P^T d = P^T B K^-1 f
 startTimerMemory(times.project, times.memoryProject1);
 computePtd(r,f);
 stopTimerMemory(times.project, times.memoryProject1);

 double r0Norm2 = real(r^r);

 double errorEstimator;
 int iter;
 double lastz2 = 0.0;
 double z2;

 for (iter=0; iter<maxiter; ++iter) {

   temps = - getTime();
   
   // Precondition z = M^-1 r
   errorEstimator = preCondition(r, z);

   if ((iter % fetiInfo->numPrint() == 0) && (fetiInfo->numPrint() > 0)) {
     double rError = (r0Norm2==0.0) ? sqrt(real(r^r)) : 
                                      sqrt(real(r^r)/r0Norm2);
     fprintf(stderr," %3d Relative Primal Error = %e Relative Dual Error"
     " = %e\n",iter, sqrt(errorEstimator/pseudoFNormSq), rError);
   }

   // Test for convergence or maximum number of iterations, this
   // way, on exit of the loop we are guaranteed that alpha is 
   // correct
   if (((z2 = errorEstimator) < epsilon2*pseudoFNormSq) ||  
      (iter == maxiter - 1)) {
      if (z2< epsilon2*pseudoFNormSq) {
        times.converged = 1;
      }
      break;
   }
   
   // Check for stagnation
   if (fabs(z2-lastz2) < 1.0e-6*lastz2) {
     fprintf(stderr,"\n STAGNATION : \n\n Relative Primal Error Reached = "
             " %e %e %e %e \n",sqrt(z2/pseudoFNormSq), fabs(z2-lastz2),
             z2, lastz2);
     times.converged = 0;
     break;
   } 

   lastz2 = z2;
  
   // Apply F_I P
   startTimerMemory(times.project, times.memoryProject1);
   applyFP(z,Fz);
   stopTimerMemory(times.project, times.memoryProject1);

   // Orthogonalization
   orthogonalize(z,Fz,p,Fp);

   // Update
   DComplex FpFp = Fp^Fp;
   DComplex Fpr  = Fp^r;
   DComplex zeta = Fpr/FpFp;
   lambda.linAdd(zeta,p);
   r.linAdd(-zeta,Fp);

   times.reOrtho -= getTime();
   orthoAdd(p, Fp, FpFp); 
   times.reOrtho += getTime();

   temps += getTime();
   tLoopAvg += temps;
   if (temps<tLoopMin)
     tLoopMin = temps;
   if (temps>tLoopMax)
     tLoopMax = temps;

 }

 times.numIter = iter;
 times.iterations[0].numFetiIter = iter;
 times.iterations[0].finalPrimal = sqrt(errorEstimator/pseudoFNormSq);
 double rError = (r0Norm2==0.0) ? sqrt(real(r^r)) :
                                  sqrt(real(r^r)/r0Norm2);
 times.iterations[0].finalDual   = rError;

 f = ft;
 recoverU(f,lambda,u);

 stopTimerMemory(times.solve, times.memorySolve);

 if (iter) {
   fprintf(stderr," --------------------------------------\n");
   fprintf(stderr," ... Min time for one iteration : %12.4e s ...\n",
           tLoopMin/1000.0);
   fprintf(stderr," ... Avg time for one iteration : %12.4e s ...\n",
           tLoopAvg/(1000.0*iter));
   fprintf(stderr," ... Max time for one iteration : %12.4e s ...\n",
           tLoopMax/1000.0);
   fprintf(stderr," --------------------------------------\n");
 }

}


void
FetiHAxiSolver::computePtd(DistrComplexVector &r, DistrComplexVector&f) {

 int totalFourier = 2*numModes+1;
 int Fourier;

 int nP = (nsub>numModes) ? nsub : 2*numModes+1;

 execParal1R(nP, this, &FetiHAxiSolver::multQtBK, f);

 execParal(totalFourier, this, &FetiHAxiSolver::assembleCoarseVector);

 execParal(totalFourier, this, &FetiHAxiSolver::solveCoarsePb);

 execParal1R(nP, this, &FetiHAxiSolver::subtractQw,f);

 startTimerMemory(times.sAndJ, times.memorySAndJ);
 timedParal2R(times.solveAndJump, nP, this, &FetiHAxiSolver::multBK, r, f);
 stopTimerMemory(times.sAndJ, times.memorySAndJ);

 for (Fourier=0; Fourier<totalFourier; ++Fourier) {
   execParal2R(nsub, this, &FetiHAxiSolver::sendInterf, r, Fourier);
   execParal2R(nsub, this, &FetiHAxiSolver::InterfDiff, r, Fourier);
 }

}


void
FetiHAxiSolver::zeroLocalCVec(int iP) {
 
 if (nsub>numModes)
   mdAxi[iP]->zeroLocalCVec();
 else
   for (int iSub=0; iSub<nsub; ++iSub)
     mdAxi[iSub]->zeroLocalCVec(iP, iP+1);

}


void 
FetiHAxiSolver::multQtBK(int iP, DistrComplexVector &f) {

// Compute and store locally for each mode : 
// Q^T B K^-1 ( - f )

 if (nsub>numModes) {
   mdAxi[iP]->multLocalQtBK(f.subData(mdAxi[iP]->localSubNum()));
 }
 else {
   for (int iSub=0; iSub<nsub; ++iSub)
     mdAxi[iSub]->multLocalQtBK(f.subData(mdAxi[iSub]->localSubNum()),iP,iP+1); 
 }

}


void
FetiHAxiSolver::assembleCoarseVector(int Fourier) {

// Assemble all the local contribution for each Fourier 

 int i, iSub; 

 CoarseVector[Fourier]->zero();
 for (iSub=0; iSub<nsub; ++iSub) {
   int numCDofs = mdAxi[iSub]->numCoarseDofs();
   int *dofs = mdAxi[iSub]->coarseEqNums;
   DComplex *localCVec = mdAxi[iSub]->getLocalCVec(Fourier);
   for (i=0; i<numCDofs ; ++i)
     (*CoarseVector[Fourier])[dofs[i]] += localCVec[i];
 }

}
               

void
FetiHAxiSolver::solveCoarsePb(int Fourier) {

// Compute (CoarseSolver)^-1 CoarseVector
// Overwrite the solution in CoarseVector

 int numF;
 switch (fetiInfo->nonLocalQ) {
   default:
   case 0:
     numF = (Fourier%2==0) ? Fourier/2 : (Fourier+1)/2; 
     break;
   case 1:
     numF = Fourier;
     break;
 }
   
 if (CoarseSolver[numF]==0)
    return;

 CoarseSolver[numF]->reSolve(*CoarseVector[Fourier]);
 
}


void
FetiHAxiSolver::subtractQw(int iP, DistrComplexVector &f) {

// Compute f -= Q . CoarseVector

 if (nsub>numModes) {
   mdAxi[iP]->subtractQw(f.subData(mdAxi[iP]->localSubNum()), 
                           CoarseVector);
 } 
 else {
   for (int iSub=0; iSub<nsub; ++iSub) 
     mdAxi[iSub]->subtractQw(f.subData(mdAxi[iSub]->localSubNum()),
                           CoarseVector, 0, iP, iP+1);
 }

}


void
FetiHAxiSolver::multBK(int iP, DistrComplexVector &r, DistrComplexVector &f) {

// Compute for each subdomain 
//    r = B K^-1 f

 if (nsub>numModes) {
   int sn = mdAxi[iP]->localSubNum();  
   DComplex *buff = new DComplex[(2*numModes+1)*mdAxi[iP]->localLen()];
   mdAxi[iP]->multKf(buff,f.subData(sn));
   mdAxi[iP]->computeBw(r.subData(sn),buff);
   delete[] buff;
 }
 else {
   for (int iSub=0; iSub<nsub; ++iSub) {
     int sn = mdAxi[iSub]->localSubNum();  
     DComplex *buff = (DComplex *) dbg_alloca(sizeof(DComplex)*
                       mdAxi[iSub]->localLen());
     mdAxi[iSub]->multKf(buff,f.subData(sn), iP, 0, iP, iP+1);
     mdAxi[iSub]->computeBw(r.subData(sn),buff, iP, iP+1);
   }
 }

}


void
FetiHAxiSolver::sendInterf(int iSub, DistrComplexVector &r, int &Fourier) {

 int sn = mdAxi[iSub]->localSubNum();
 mdAxi[iSub]->sendInterf(r.subData(sn)+Fourier*mdAxi[iSub]->interfLen());

}


void
FetiHAxiSolver::InterfDiff(int iSub, DistrComplexVector &r, int &Fourier) {

 int sn = mdAxi[iSub]->localSubNum();
 mdAxi[iSub]->interfaceJump(r.subData(sn)+Fourier*mdAxi[iSub]->interfLen());

}


double
FetiHAxiSolver::preCondition(DistrComplexVector &r, 
                DistrComplexVector &Mr) {

 startTimerMemory( times.precond, times.memoryPrecond );

// Multiply the vector r by the preconditioner M
// Compute the error estimator

 int nP = (nsub>numModes) ? nsub : 2*numModes+1;
 timedParal4R(times.preconditioner, nP, this, &FetiHAxiSolver::multKbb, 
              r, Mr, *deltaU, *deltaF);

 int Fourier;
 int totalFourier = 2*numModes+1;

 for (Fourier=0; Fourier<totalFourier; ++Fourier) {
    timedParal2R(times.preconditioner,nsub,this,&FetiHAxiSolver::sendInterf,
                 Mr, Fourier);
    timedParal2R(times.preconditioner,nsub,this,&FetiHAxiSolver::InterfDiff, 
                 Mr, Fourier);
 }

 DComplex primalResidual = 0.0;
 DComplex *subDots = (DComplex *) dbg_alloca(sizeof(DComplex)*nsub);

 // Compute each partial true norm
 timedParal(times.preconditioner, nsub, this, &FetiHAxiSolver::normDeltaF, 
            subDots, deltaF);

 // Now sum each partial norm over all subdomains
 int i;
 for(i = 0; i < nsub; ++i) 
   primalResidual += subDots[i];

 // No preconditioner case
 if (fetiInfo->precno == FetiInfo::noPrec) {
    Mr = r;
 }

 stopTimerMemory( times.precond, times.memoryPrecond );

 // return true primal residual
 return std::abs(primalResidual);

}


void
FetiHAxiSolver::multKbb(int iP,DistrComplexVector &r,
                DistrComplexVector &Mr, DistrComplexVector &dU,
                DistrComplexVector &dF) {

 if (nsub>numModes) {
   int sn = mdAxi[iP]->localSubNum();
   mdAxi[iP]->multKbbC(r.subData(sn),Mr.subData(sn),
                         dU.subData(sn),dF.subData(sn));  
 }
 else {
   for (int iSub=0; iSub<nsub; ++iSub) {
     int sn = mdAxi[iSub]->localSubNum();
     mdAxi[iSub]->multKbbC(r.subData(sn),Mr.subData(sn),
                           dU.subData(sn),dF.subData(sn), iP, iP+1);  
   }
 }

}


void
FetiHAxiSolver::normDeltaF(int iSub, DComplex *subDots, 
                DistrComplexVector *dF) {

 int Fourier;
 int totalFourier = 2*numModes+1;
 int localSize = mdAxi[iSub]->localLen();
 int offSet = 0;

 DComplex *buff = dF->subData(iSub);
 subDots[iSub] = DComplex(0.0,0.0);

 for (Fourier=0; Fourier<totalFourier; ++Fourier) {
   mdAxi[iSub]->sendDeltaF(buff + offSet);
   subDots[iSub] += mdAxi[iSub]->collectAndDotDeltaF(buff + offSet);
   offSet += localSize;
 }
 
}


void
FetiHAxiSolver::applyFP(DistrComplexVector &z,DistrComplexVector &Fz) {

 int Fourier;
 int totalFourier = 2*numModes+1;
 
 int nP = (nsub>numModes) ? nsub : 2*numModes+1;

 execParal1R(nP, this, &FetiHAxiSolver::multTransposeBKQ, z);

 execParal(totalFourier, this, &FetiHAxiSolver::assembleCoarseVector);

 execParal(totalFourier, this, &FetiHAxiSolver::solveCoarsePb);

 startTimerMemory(times.sAndJ, times.memorySAndJ);
 timedParal2R(times.solveAndJump, nP, this, &FetiHAxiSolver::finishFPz,
              z, Fz);
 stopTimerMemory(times.sAndJ, times.memorySAndJ);

 for (Fourier=0; Fourier<totalFourier; ++Fourier) {
    execParal2R(nsub, this, &FetiHAxiSolver::sendInterf, Fz, Fourier);
    execParal2R(nsub, this, &FetiHAxiSolver::InterfDiff, Fz, Fourier);
 }

}


void
FetiHAxiSolver::multTransposeBKQ(int iP, DistrComplexVector &z) {

// Compute and store locally for each mode :
// Q^T B K^-1 B^T ( - z )

 if (nsub>numModes) {
   int sn = mdAxi[iP]->localSubNum();
   mdAxi[iP]->multTransposeBKQ(z.subData(sn));
 }
 else {
   for (int iSub=0; iSub<nsub; ++iSub) {
     int sn = mdAxi[iSub]->localSubNum();
     mdAxi[iSub]->multTransposeBKQ(z.subData(sn), iP, iP+1);
   }
 }
 
}


void
FetiHAxiSolver::finishFPz(int iP, DistrComplexVector &z,
                          DistrComplexVector &FPz) {

 if (nsub>numModes) {
   int sn = mdAxi[iP]->localSubNum();
   mdAxi[iP]->finishFPz(CoarseVector,z.subData(sn),FPz.subData(sn));
 }
 else {
   for (int iSub=0; iSub<nsub; ++iSub) {
     int sn = mdAxi[iSub]->localSubNum();
     mdAxi[iSub]->finishFPz(CoarseVector,z.subData(sn),FPz.subData(sn),iP,iP+1);
   }
 }

}


void
FetiHAxiSolver::orthogonalize(DistrComplexVector &r, DistrComplexVector &Fr,
               DistrComplexVector &p, DistrComplexVector &Fp) {

 times.reOrtho -= getTime();

 if (oSet->numDir() == 0) {
   p = r;
   Fp = Fr;
 }  
 else {

   DComplex *hp = (DComplex *) dbg_alloca(halfSize*sizeof(DComplex));
   DComplex *hOp= (DComplex *) dbg_alloca(halfSize*sizeof(DComplex));

   DComplex *hFp = (DComplex *) dbg_alloca(halfSize*sizeof(DComplex));
   DComplex *hOFp= (DComplex *) dbg_alloca(halfSize*sizeof(DComplex));

   DistrComplexVector *vec2 = 0;
   timedParal(times.orthogonalize,nsub,this,
              &FetiHAxiSolver::gatherHalfInterface,&r,vec2,hp,hOp);

   vec2 = 0;
   timedParal(times.orthogonalize,nsub,this,
              &FetiHAxiSolver::gatherHalfInterface,&Fr,vec2,hFp,hOFp);

   oSet->orthogonalizeTimed(times.orthogonalize,hp,hFp,hOp,hOFp);

   int Fourier;
   int totalFourier = 2*numModes+1;

   for (Fourier=0; Fourier<totalFourier; ++Fourier) {
     timedParal(times.orthogonalize,nsub,this,
               &FetiHAxiSolver::scatterHalfInterface, hOp, &p, &Fourier);
     timedParal2R(times.orthogonalize,nsub, this, 
                 &FetiHAxiSolver::rebuildInterface, p, Fourier);
   }

   for (Fourier=0; Fourier<totalFourier; ++Fourier) {
     timedParal(times.orthogonalize, nsub, this,
               &FetiHAxiSolver::scatterHalfInterface, hOFp, &Fp, &Fourier);
     timedParal2R(times.orthogonalize, nsub, this, 
                 &FetiHAxiSolver::rebuildInterface, Fp, Fourier);
   }
 
 }
  
 times.reOrtho += getTime();

}


void
FetiHAxiSolver::gatherHalfInterface(int iSub, DistrComplexVector *v1, 
                DistrComplexVector *v2, DComplex *v3, DComplex *v4) {

 DComplex *vec1 =  v1->subData(mdAxi[iSub]->localSubNum());
 DComplex *vec2 =  (v2) ? v2->subData(mdAxi[iSub]->localSubNum()) : 0;

 DComplex *vec3 = v3 + mdAxi[iSub]->halfOffset;
 DComplex *vec4 = v4 + mdAxi[iSub]->halfOffset;

 int Fourier;
 int totalFourier = 2*numModes+1;
 int interfSize = mdAxi[iSub]->interfLen();
 int interfHalfSize = mdAxi[iSub]->halfInterfLen();

 if (vec2) 
   for (Fourier=0; Fourier<totalFourier; ++Fourier) 
      mdAxi[iSub]->getHalfInterf(vec1+Fourier*interfSize,
                   vec3+Fourier*interfHalfSize,vec2+Fourier*interfSize,
                   vec4+Fourier*interfHalfSize);
 else
   for (Fourier=0; Fourier<totalFourier; ++Fourier)
       mdAxi[iSub]->getHalfInterf(vec1+Fourier*interfSize, 
                    vec3+Fourier*interfHalfSize);

}


void
FetiHAxiSolver::scatterHalfInterface(int iSub, DComplex *v1, 
                DistrComplexVector* v2, int *Fourier) {

 DComplex *vec1 = v1 + mdAxi[iSub]->halfOffset;
 DComplex *vec2 = v2->subData( mdAxi[iSub]->localSubNum() );

 int interfHalfSize = mdAxi[iSub]->halfInterfLen();
 int interfSize = mdAxi[iSub]->interfLen();

 mdAxi[iSub]->scatterHalfInterf( vec1+(*Fourier)*interfHalfSize, 
                                 vec2+(*Fourier)*interfSize );
 mdAxi[iSub]->sendInterf(vec2+(*Fourier)*interfSize);

}


void
FetiHAxiSolver::rebuildInterface(int iSub, DistrComplexVector &v,
                int &Fourier) {

 mdAxi[iSub]->rebuildInterf(v.subData(mdAxi[iSub]->localSubNum())+
                            Fourier*mdAxi[iSub]->interfLen());

}


void
FetiHAxiSolver::orthoAdd(DistrComplexVector &p, DistrComplexVector &Fp, 
                DComplex pFp) {

 DComplex *hp  = (DComplex *) dbg_alloca(halfSize*sizeof(DComplex));
 DComplex *hFp = (DComplex *) dbg_alloca(halfSize*sizeof(DComplex));

 timedParal(times.orthogonalize, nsub, this,
            &FetiHAxiSolver::gatherHalfInterface, &p, &Fp, hp, hFp);

 times.memoryOSet -= threadManager->memoryUsed();  
 oSet->orthoAddTimed( times.orthogonalize, hp, hFp, pFp*DComplex(0.5,0.0) );
 times.memoryOSet += threadManager->memoryUsed();

}


void
FetiHAxiSolver::recoverU(DistrComplexVector &f, 
                DistrComplexVector &z, DistrComplexVector &u) {

 int totalFourier = 2*numModes+1;
 int nP = (nsub>numModes) ? nsub : 2*numModes+1;

 execParal2R(nP, this, &FetiHAxiSolver::subtractBt, z, f);

 execParal1R(nP, this, &FetiHAxiSolver::multQtBK, f);

 execParal(totalFourier, this, &FetiHAxiSolver::assembleCoarseVector);

 execParal(totalFourier, this, &FetiHAxiSolver::solveCoarsePb);

 execParal1R(nP, this, &FetiHAxiSolver::subtractQw,f);

 execParal2R(nP, this, &FetiHAxiSolver::localSolution, f, u);

}


void
FetiHAxiSolver::subtractBt(int iP, DistrComplexVector &z,
                           DistrComplexVector &f) {

// Compute for each subdomain
//    f -= B^T z 

 if (nsub>numModes) {
   int sn = mdAxi[iP]->localSubNum();
   mdAxi[iP]->subtractBt(z.subData(sn),f.subData(sn));
 }
 else {
   for (int iSub=0; iSub<nsub; ++iSub) {
     int sn = mdAxi[iSub]->localSubNum();
     mdAxi[iSub]->subtractBt(z.subData(sn), f.subData(sn), iP, iP+1);
   }
 }

}


void
FetiHAxiSolver::localSolution(int iP, DistrComplexVector &f,
                              DistrComplexVector &u) {

// Compute for each subdomain
//    u = K^-1  f 

 if (nsub>numModes) {
   int sn = mdAxi[iP]->localSubNum();
   mdAxi[iP]->multKf(u.subData(sn),f.subData(sn));
 }
 else {
   for (int iSub=0; iSub<nsub; ++iSub) {
     int sn = mdAxi[iSub]->localSubNum();
     mdAxi[iSub]->multKf(u.subData(sn),f.subData(sn), 0, 0, iP, iP+1);
   }
 }

}


Timings&
FetiHAxiSolver::getTimers() {

 return times;

}


double
FetiHAxiSolver::getSolutionTime() {

 return times.solve;

}


