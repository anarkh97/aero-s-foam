#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/SubDomain.h>
#include <Driver.d/PolygonSet.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <HelmAxi.d/AxiHElem.h>
#include <HelmAxi.d/FourierHelmBCs.h>
#include <Feti.d/FetiInfo.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <HelmAxi.d/MDAxiData.h>
#include <Threads.d/PHelper.h>
#include <Feti.d/DistrVector.h>
#include <HelmAxi.d/FetiHAxi.d/DistrComplexVector.h>
#include <HelmAxi.d/FetiHAxi.d/GCRC.h>
#include <HelmAxi.d/FetiHAxi.d/FetiHAxi.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>


void 
FetiHAxiSolver::solveMPC(DistrComplexVector &f, ComplexVector &g, 
                DistrComplexVector &u) {

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

 buildMPCSolver();

 startTimerMemory(times.solve, times.memorySolve);

 times.memoryDV -= threadManager->memoryUsed();
 DistrComplexVector rl(interface);
 DistrComplexVector zl(interface);
 DistrComplexVector pl(interface);
 DistrComplexVector Fzl(interface);
 DistrComplexVector Fpl(interface);
 DistrComplexVector lambda(interface);
 DistrComplexVector ft(internal);
 deltaU = new DistrComplexVector(internal);
 deltaU->zero();
 deltaF = new DistrComplexVector(internal);
 deltaF->zero();
 times.memoryDV += threadManager->memoryUsed();

 ft = f;

 int Fourier;
 int totalFourier = 2*numModes+1;
 ComplexVector **lambdaVec = new ComplexVector *[totalFourier];
 for (Fourier=0; Fourier<totalFourier; ++Fourier)
   lambdaVec[Fourier] = new ComplexVector(coarseEqs->size());

 ComplexVector mu(g);
 mu.zero();

 times.memoryOSet -= threadManager->memoryUsed();
 oSet = new GCRC(halfSize, fetiInfo->maxortho);
 times.memoryOSet += threadManager->memoryUsed();

 // Compute the square of a pseudo-norm of fVec
 double pseudoFNormSq = real(f^f);

 // Set some initial vectors at zero
 lambda.zero();

 // Compute r0 = P^T d 

 startTimerMemory(times.project, times.memoryProject1);
 computePtd(rl,f,g,mu,lambdaVec);
 stopTimerMemory(times.project, times.memoryProject1);

 double r0Norm2 = real(rl^rl);

 double errorEstimator;
 int iter;
 double lastz2 = 0.0;
 double z2;

 for (iter=0; iter<maxiter; ++iter) {

   temps = - getTime();

   // At each iteration, the MPC is exactly satisfied.
   // Hence, this gives no contribution for the error.
  
   // Precondition z = M^-1 r with identity for rm
   errorEstimator = preCondition(rl, zl);

   if ((iter % fetiInfo->numPrint() == 0) && (fetiInfo->numPrint() > 0)) {
     double rError = (r0Norm2==0.0) ? sqrt(real(rl^rl)) :
                                      sqrt(real(rl^rl)/r0Norm2);
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

   // Apply F_I P for the MPC case
   startTimerMemory(times.project, times.memoryProject1);
   applyFPMPC(zl,Fzl,lambdaVec);
   stopTimerMemory(times.project, times.memoryProject1);

   // Orthogonalization
   orthogonalize(zl,Fzl,pl,Fpl);

   // Update
   DComplex FpFp = Fpl^Fpl;
   DComplex Fpr  = Fpl^rl;
   DComplex zeta = Fpr/FpFp;
   lambda.linAdd(zeta,pl);
   rl.linAdd(-zeta,Fpl);

   times.reOrtho -= getTime();
   orthoAdd(pl, Fpl, FpFp);
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
 double rError = (r0Norm2==0.0) ? sqrt(real(rl^rl)) :
                                  sqrt(real(rl^rl)/r0Norm2);
 times.iterations[0].finalDual   = rError;

 f = ft;
 recoverUMPC(f,lambda,g,u,mu,lambdaVec);

 for (Fourier=0; Fourier<totalFourier; ++Fourier)
   delete lambdaVec[Fourier];
 delete[] lambdaVec;

 stopTimerMemory(times.solve, times.memorySolve);

// FILE *fid = fopen("VectorMu","w");
// mu.print();

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
FetiHAxiSolver::buildMPCSolver() {

 // For a matrix
 //
 //   |  A_1     0   B_1 |
 //   |   0     A_2  B_2 |
 //   | B_1^T  B_2^T  C  |
 //
 // MPCSolver = C - B_1^T A_1^-1 B_1 - B_2^T A_2^-1 B_2

 int totalFourier = 2*numModes+1;
 int nT = threadManager->numThr();

 startTimerMemory(times.coarse2, times.memoryPCtFPC); 

 execParal(totalFourier, this, &FetiHAxiSolver::startSchur);

 MPCSolver = new SkyMatrixC(MPCSize, 0.0);

 execParal(nT, this, &FetiHAxiSolver::finishSchur);

 stopTimerMemory(times.coarse2, times.memoryPCtFPC);

 fprintf(stderr," >>>> Assembly of MPCSolver : %e s\n", times.coarse2/1000.0);

 // Factorize the matrix solver on Lagrange multipliers for MPC

 times.pfactor2 -= getTime();
 MPCSolver->parallelFactor();
 times.pfactor2 += getTime();

 fprintf(stderr," >>>> Factorization of MPC Solver : %e s \n", 
         times.pfactor2/1000.0); 

}


void 
FetiHAxiSolver::startSchur(int Fourier) {

// CoarseMPC stores after the routine 
//       { - (Q^T B K^-1 B^T Q)^-1 } * { - Q^T B K^-1 C^T }

 int lsize;
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

 lsize = (CoarseSolver[numF]==0) ? 0 : CoarseSolver[numF]->dim();

 if (lsize)
   CoarseSolver[numF]->reSolve(CoarseMPC[Fourier]); 
 else  
   CoarseMPC[Fourier]->setNewSize(0);

}


void
FetiHAxiSolver::finishSchur(int iT) {

 int i, j, iSub;
 int locMPCSize;
 int nT = threadManager->numThr();
 int Fourier;
 int totalFourier = 2*numModes+1;

 DComplex buffer(0.0, 0.0);
 DComplex *pointer = (DComplex *) dbg_alloca(sizeof(DComplex)*4);
 pointer[2] = DComplex(0.0, 0.0);
 pointer[3] = DComplex(0.0, 0.0);
 //int *dofs = (int *) dbg_alloca(sizeof(int)*2);

 for (iSub=0; iSub<nsub; ++iSub) {
   locMPCSize = mdAxi[iSub]->getMPCSize();
   if (locMPCSize==0)
      continue;
   int numCDofs = mdAxi[iSub]->numCoarseDofs();
   int *mapping = mdAxi[iSub]->globalToLocalMPC;
   if (numCDofs) { 
     for (i=iT; i<MPCSize; i+=nT) {

       if (mapping[i]==-1)
           continue;

       for (j=i; j<MPCSize; ++j) {

         if (mapping[j]==-1)
           buffer = DComplex(0.0, 0.0); 
         else
           buffer = -(*mdAxi[iSub]->CKCt)[mapping[i]][mapping[j]-mapping[i]];

         for (Fourier=0; Fourier<totalFourier; ++Fourier) {
           for (int k=0; k<numCDofs; ++k) {
             buffer += mdAxi[iSub]->CKQ[Fourier][k*locMPCSize+mapping[i]]*
                       (*CoarseMPC[Fourier])[mdAxi[iSub]->coarseEqNums[k]][j];
           }
         }

         MPCSolver->addPoint(buffer, i, j);

       }
     }
   }
   else {
     for (i=iT; i<MPCSize; i+=nT) {

       if (mapping[i]==-1)
           continue;

       for (j=i; j<MPCSize; ++j) {
         if (mapping[j]==-1)
           continue;
         MPCSolver->addPoint(-(*mdAxi[iSub]->CKCt)[mapping[i]][mapping[j]-
                             mapping[i]], i, j);
       }

     }
   }
 }

}


void 
FetiHAxiSolver::computePtd(DistrComplexVector &rl, DistrComplexVector &f,
                           ComplexVector &g, ComplexVector &muVector,
                           ComplexVector **lambdaVec) {

 int totalFourier = 2*numModes+1;
 int Fourier;

 int nP = (nsub>numModes) ? nsub : 2*numModes+1;
 int nT = threadManager->numThr();

 execParal1R(nP, this, &FetiHAxiSolver::multQtBK, f);

 execParal(totalFourier, this, &FetiHAxiSolver::assembleCoarseVector);

 ComplexVector rm(g);
 double coeff = 1.0;

 // rm =  g - \sum_{subd,mode} C_mode^{s^T} K_mode^{s^-1} f_mode^s

 execParal(nT, this, &FetiHAxiSolver::assembleMPCVector, rm.data(), &coeff);

 // Build the RHS for the Schur equation

 execParal(nT, this, &FetiHAxiSolver::copyMu, rm.data(), muVector.data());

 execParal(nT, this, &FetiHAxiSolver::buildRHSSchur,muVector.data());

 // Solve the Schur equation
 MPCSolver->reSolve(muVector);
 
 // Recover the eliminated variables
 execParal(totalFourier, this, &FetiHAxiSolver::recoverSchur, lambdaVec,
           &muVector); 

 // Final treatment for the residual on mu

 execParal(nP, this, &FetiHAxiSolver::zeroLocalCVec);

 execParal(nP, this, &FetiHAxiSolver::addCKQ, lambdaVec, &muVector);

 coeff = 1.0;

 execParal(nT, this, &FetiHAxiSolver::assembleMPCVector, rm.data(), &coeff);

 fprintf(stderr," --------------------------------------------\n");
 fprintf(stderr," Norm of residual on mu = %e \n",sqrt(rm.squareNorm()));
 fprintf(stderr," --------------------------------------------\n");

 // Final treatment for the residual on lambda

 // The next routine overwrites f
 execParal(nP, this, &FetiHAxiSolver::subtractQwMPC, &f, lambdaVec,
           &muVector);

 startTimerMemory(times.sAndJ, times.memorySAndJ);
 timedParal2R(times.solveAndJump, nP, this, &FetiHAxiSolver::multBK, rl, f);
 stopTimerMemory(times.sAndJ, times.memorySAndJ);

 for (Fourier=0; Fourier<totalFourier; ++Fourier) {
   execParal2R(nsub, this, &FetiHAxiSolver::sendInterf, rl, Fourier);
   execParal2R(nsub, this, &FetiHAxiSolver::InterfDiff, rl, Fourier);
 }

}


void
FetiHAxiSolver::assembleMPCVector(int iT, DComplex *v, double *alpha) {

// Assemble all the local contribution for each Fourier
//     \sum_{subd} C K^-1 (-f)

 int i, iSub;
 int locMPCSize;
 int nT = threadManager->numThr();
 int Fourier;
 int totalFourier = 2*numModes+1;

 for (iSub=0; iSub<nsub; ++iSub) {
   locMPCSize = mdAxi[iSub]->getMPCSize();
   if (locMPCSize==0)
     continue;
   int *mapping = mdAxi[iSub]->globalToLocalMPC;
   int numCDofs = mdAxi[iSub]->numCoarseDofs();
   for (i=iT; i<MPCSize; i+=nT) {
     if (mapping[i]==-1)
       continue;
     for (Fourier=0; Fourier<totalFourier; ++Fourier) {
       DComplex *localVec = mdAxi[iSub]->getLocalCVec(Fourier);
       v[i] += (*alpha)*localVec[mapping[i]+numCDofs];
     }
   }
 }

}


void 
FetiHAxiSolver::copyMu(int iT, DComplex *rm, DComplex *muVector) {

 int i;
 int nT = threadManager->numThr();  

 for (i=iT; i<MPCSize; i+=nT) {
   muVector[i] = - rm[i];
 }

}


void
FetiHAxiSolver::buildRHSSchur(int iT, DComplex *buffer) {

// Compute (CoarseMPC[Fourier])^T * CoarseVector[Fourier]
// ~ or ~ C K^-1 B^T Q (Q^T B K^-1 B^T Q)^-1 ( CoarseVector[Fourier] )
// ~ when called from computePtd ~
//        C K^-1 B^T Q (Q^T B K^-1 B^T Q)^-1 ( Q^T B K^-1 (-f) )
// ~ when called from applyFP ~
//        C K^-1 B^T Q (Q^T B K^-1 B^T Q)^-1 ( Q^T B K^-1 B^T (-z) )
//

 int i,k;
 int size;
 int nT = threadManager->numThr();  
 int Fourier;
 int totalFourier = 2*numModes+1;

 size = CoarseVector[0]->size();

 for (i=iT; i<MPCSize; i+=nT) {
   for (Fourier=0; Fourier<totalFourier; ++Fourier) {
     for (k=0; k<size; ++k) {
       buffer[i] += (*CoarseMPC[Fourier])[k][i]*(*CoarseVector[Fourier])[k];
     }
   }
 }

}


void
FetiHAxiSolver::recoverSchur(int Fourier, ComplexVector **lambda,
                ComplexVector *mu) {

 int i,k;
 int lsize;
 int numF;

 lambda[Fourier]->copy(*CoarseVector[Fourier]);

 switch (fetiInfo->nonLocalQ) {
   default:
   case 0:
     numF = (Fourier%2==0) ? Fourier/2 : (Fourier+1)/2;
     break;
   case 1:
     numF = Fourier;
     break;
 }

 lsize = (CoarseSolver[numF]==0) ? 0 : CoarseSolver[numF]->dim();
 if (lsize == 0)
   return;

 CoarseSolver[numF]->reSolve((*lambda[Fourier])); 

 for (i=0; i<lsize; ++i) {
   for (k=0; k<MPCSize; ++k) {
     (*lambda[Fourier])[i] -= (*CoarseMPC[Fourier])[i][k]*(*mu)[k];
   }
 }

}


void
FetiHAxiSolver::subtractQwMPC(int iP, DistrComplexVector *f,
                ComplexVector **lambdaVec, ComplexVector *mu) {

 if (nsub>numModes) { 
   mdAxi[iP]->subtractQw(f->subData(mdAxi[iP]->localSubNum()),
                           lambdaVec, mu);
 }
 else {
   for (int iSub=0; iSub<nsub; ++iSub)
     mdAxi[iSub]->subtractQw(f->subData(mdAxi[iSub]->localSubNum()),
                         lambdaVec, mu, iP, iP+1);
 }

}


void 
FetiHAxiSolver::applyFPMPC(DistrComplexVector &zl, DistrComplexVector &Fzl,
                           ComplexVector **lambdaVec) { 

 int Fourier;
 int totalFourier = 2*numModes+1;

 int nP = (nsub>numModes) ? nsub : 2*numModes+1;
 int nT = threadManager->numThr();

 execParal1R(nP, this, &FetiHAxiSolver::multTransposeBKQ, zl);

 execParal(totalFourier, this, &FetiHAxiSolver::assembleCoarseVector);

 // Build the RHS for the Schur equation

 ComplexVector muVector(MPCSize, DComplex(0.0,0.0));
 double coeff = -1.0;

 execParal(nT,this,&FetiHAxiSolver::assembleMPCVector,muVector.data(),&coeff);

 execParal(nT, this, &FetiHAxiSolver::buildRHSSchur,muVector.data());

 // Solve the Schur equation
 MPCSolver->reSolve(muVector);

 // Recover the eliminated variables
 execParal(totalFourier, this, &FetiHAxiSolver::recoverSchur, lambdaVec,
           &muVector);

 // Final treatment for the residual on mu

 execParal(nP, this, &FetiHAxiSolver::addCKQ, lambdaVec, &muVector);

 ComplexVector Fzm(MPCSize, DComplex(0.0,0.0));

 execParal(nT, this, &FetiHAxiSolver::assembleMPCVector, Fzm.data(), &coeff);

 fprintf(stderr," --------------------------------------------- \n");
 fprintf(stderr," Norm of residual on mu = %e \n",sqrt(Fzm.squareNorm()));
 fprintf(stderr," --------------------------------------------- \n");

 // Final treatment for the residual on lambda

 startTimerMemory(times.sAndJ, times.memorySAndJ);
 timedParal4R(times.solveAndJump, nP, this, &FetiHAxiSolver::finishFzl,
              zl, muVector, lambdaVec, Fzl); 
 stopTimerMemory(times.sAndJ, times.memorySAndJ);

 for (Fourier=0; Fourier<totalFourier; ++Fourier) {
   execParal2R(nsub, this, &FetiHAxiSolver::sendInterf, Fzl, Fourier);
   execParal2R(nsub, this, &FetiHAxiSolver::InterfDiff, Fzl, Fourier);
 }

}


void
FetiHAxiSolver::addCKQ(int iP, ComplexVector **lambdaVec,
                      ComplexVector *mu) {

 if (nsub>numModes) 
   mdAxi[iP]->addCKQw(lambdaVec, mu);
 else
   for (int iSub=0; iSub<nsub; iSub++)
     mdAxi[iSub]->addCKQw(lambdaVec, mu, iP, iP+1);

}


void
FetiHAxiSolver::finishFzl(int iP, DistrComplexVector &zl, 
              ComplexVector &mu, ComplexVector** &lambda, 
              DistrComplexVector &Fzl) {

 if (nsub>numModes) {
   int sn = mdAxi[iP]->localSubNum();
   mdAxi[iP]->finishFzl(zl.subData(sn), mu.data(), lambda, Fzl.subData(sn));
 }
 else {
   for (int iSub=0; iSub<nsub; ++iSub) {
     int sn = mdAxi[iSub]->localSubNum();
     mdAxi[iSub]->finishFzl(zl.subData(sn), mu.data(), lambda, Fzl.subData(sn),
                  iP, iP+1);
   }
 }

}


void
FetiHAxiSolver::recoverUMPC(DistrComplexVector &f, DistrComplexVector &z,
                         ComplexVector &g, DistrComplexVector &u, 
                         ComplexVector &muVector, ComplexVector **lambdaVec) {

 int nP = (nsub>numModes) ? nsub : 2*numModes+1;
 int nT = threadManager->numThr();
 int totalFourier = 2*numModes+1;

 execParal2R(nP, this, &FetiHAxiSolver::subtractBt, z, f);

 execParal1R(nP, this, &FetiHAxiSolver::multQtBK, f);

 execParal(totalFourier, this, &FetiHAxiSolver::assembleCoarseVector);

 // rm =  g - \sum_{subd,mode} C_mode^{s^T} K_mode^{s^-1} f_mode^s - B^{s^T} z

 ComplexVector rm(g);
 double coeff = 1.0;

 execParal(nT, this, &FetiHAxiSolver::assembleMPCVector, rm.data(), &coeff);

 // Build the RHS for the Schur equation

 execParal(nT, this, &FetiHAxiSolver::copyMu, rm.data(), muVector.data()); 

 execParal(nT, this, &FetiHAxiSolver::buildRHSSchur,muVector.data());

 // Solve the Schur equation
 MPCSolver->reSolve(muVector);

 // Recover the eliminated variables
 execParal(totalFourier, this, &FetiHAxiSolver::recoverSchur, lambdaVec,
           &muVector);

 execParal(nP, this, &FetiHAxiSolver::subtractQwMPC, &f,lambdaVec,&muVector);

 execParal2R(nP, this, &FetiHAxiSolver::localSolution, f, u);

}


