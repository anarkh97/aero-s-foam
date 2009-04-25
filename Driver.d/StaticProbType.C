#include <Timers.d/StaticTimers.h>
#include <cassert>
template < class Scalar,
           class OpSolver, 
           class VecType, 
	   class PostProcessor, 
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType, 
	      PostProcessor, ProblemDescriptor, ComplexVecType >
  ::solve()
{
 probDesc->preProcess();

 rhs = new VecType(probDesc->solVecInfo());
 sol = new VecType(probDesc->solVecInfo());
 sol->zero(); 
 if((domain->solInfo().inpc || domain->solInfo().noninpc) && domain->solInfo().type != 2) // type != 2 means not FETI solver
   sol->setn(domain->numUncon());

 //solver=0; postProcessor=0;  
 allOps=0; postProcessor=0;

 //solver = probDesc->getSolver();
 allOps = probDesc->getAllOps();
 postProcessor = probDesc->getPostProcessor();

 int padeN = domain->solInfo().padeN;  // 1 for 1-point pade or taylor, 2 for 2-point pade, etc

 if(domain->solInfo().doFreqSweep) { // multiple RHS solve
   filePrint(stderr, " ... Frequency Sweep Analysis       ... \n");
   filePrint(stderr, " ... Number of coarse freqencies = %3d     ... \n", domain->coarse_frequencies->size());
   if(padeN > 1)
     filePrint(stderr, " ... Number of extrapolated freqencies = %3d     ... \n", domain->frequencies->size());
   else 
     filePrint(stderr, " ... Number of interpolated freqencies = %3d     ... \n", domain->frequencies->size());
   int nRHS = domain->solInfo().nFreqSweepRHS;
   if(nRHS==0) filePrint(stderr, " *** ERROR: nRHS = 0 \n");

   // some initialization
   probDesc->setIWaveDir(0); 
   VecType **sol_prev = new VecType * [(nRHS+1)*padeN];
   for(int i = 0; i < (nRHS+1)*padeN; ++i)
     sol_prev[i] = new VecType(probDesc->solVecInfo());
   double *h = new double[padeN];

   //--- UH --- Data for Pade Lanczos
   double *PadeLanczos_VtKV = 0, *PadeLanczos_Vtb = 0;
   VecType **PadeLanczos_solprev = 0;
   if (domain->solInfo().freqSweepMethod == SolverInfo::PadeLanczos) {
     filePrint(stderr, " ... Describe implementation PadeLanczos ... \n");
     //--- Allocate new set of pointers to have them without any gap
     PadeLanczos_solprev = new VecType * [nRHS*padeN];
     for (int ii = 0; ii < padeN; ++ii) {
       for (int iRHS = 0; iRHS < nRHS; ++iRHS)
         PadeLanczos_solprev[iRHS + ii*nRHS] = sol_prev[iRHS + ii*(nRHS+1)];
     }
   } // if (domain->solInfo().freqSweepMethod == SolverInfo::PadeLanczos)
   //--- UH --- Data for Pade Lanczos
   
   // loop over coarse grid
   bool first_time = true;
   double w0;
   int count = 0;
   while(domain->coarse_frequencies->size() > 0) {
     int offset = count*(nRHS+1);
     sol_prev[offset]->zero();
     domain->isCoarseGridSolve = true;
     double wc = domain->coarse_frequencies->front();
     if(count == 0) w0 = wc;

     // if this isn't the first coarse solve then rebuild solver with new frequency
     if(!first_time) rebuildSolver(wc); 
     else first_time = false;

     //----- UH ------
     if (domain->solInfo().freqSweepMethod == SolverInfo::PadeLanczos) {
       filePrint(stderr,"\n ... Build Subspace M-orthonormal Basis      ... \n");
       /*probDesc->*/ PadeLanczos_BuildSubspace(nRHS, PadeLanczos_solprev + count * nRHS, 
                                           PadeLanczos_solprev, count * nRHS);
     }
     else {
       for(int iRHS = 0; iRHS < nRHS; iRHS++) {
         filePrint(stderr,"\n ... Solving RHS #%3d               ... \n",iRHS);
         if(iRHS > 0) { 
           // PJSA: modified 5-29-04 now passing u and all derivatives to getFreqSweepRHS (needed for higher order sommerfeld)
           // note: sol_prev[0] = 0, sol_prev[1] = u, sol_prev[2] = u^(1), etc.
           *sol_prev[offset+iRHS] = *sol;
           probDesc->getFreqSweepRHS(rhs, sol_prev+offset, iRHS);
         } 
         else probDesc->getRHS(*rhs);
         sol->zero();
         allOps->sysSolver->solve(*rhs, *sol);
         rhs->zero();
       } // for (int iRHS = 0; iRHS < nRHS; iRHS++)
       *sol_prev[offset+nRHS] = *sol;
     }
     //----- UH ------
     // PJSA 9-22-06: fix for coupled multi point pade
     if(domain->solInfo().isCoupled) for(int i=1; i<(nRHS+1); ++i) scaleDisp(*sol_prev[offset+i]);
     bool printTimers = ((domain->coarse_frequencies->size()+domain->frequencies->size()) > 1) ? false : true;
     
     //----- UH ------
     if (domain->solInfo().freqSweepMethod == SolverInfo::PadeLanczos) {
       //--- We destroy the arrays PadeLanczos_VtKV and PadeLanczos_VtB when V is incomplete.
       //--- This is not optimal as we recompute several times some coefficients.
       //--- Outputting should be done only when the basis V is complete.
       if(PadeLanczos_VtKV) { delete[] PadeLanczos_VtKV; PadeLanczos_VtKV = 0; }
       if(PadeLanczos_Vtb) { delete[] PadeLanczos_Vtb; PadeLanczos_Vtb = 0; }

       //--- sol_prev[offset] stores the normalized solution vector !
       //--- So we need to recompute the solution.
       //--- Note that wc is referred only on the first call to this routine.
       //--- It is used to get K (as K is actually storing K - wc*wc*M)
       /*probDesc->*/ PadeLanczos_Evaluate((count+1)*nRHS, PadeLanczos_solprev, 
                                      PadeLanczos_VtKV, PadeLanczos_Vtb, wc, sol, wc*wc);
       postProcessor->staticOutput(*sol, *rhs, printTimers);
/* PJSA 10-09-08 MOVED THIS ABOVE PadeLanczos_Evaluate
       //--- We destroy the arrays PadeLanczos_VtKV and PadeLanczos_VtB when V is incomplete.
       //--- This is not optimal as we recompute several times some coefficients.
       //--- Outputting should be done only when the basis V is complete.
       if (count + 1 < padeN) {
         delete[] PadeLanczos_VtKV;
         PadeLanczos_VtKV = 0;
         delete[] PadeLanczos_Vtb;
         PadeLanczos_Vtb = 0;
       } // if (count + 1 < padeN)
*/
     }
     else {
       postProcessor->staticOutput(*sol_prev[offset+1], *rhs, printTimers);
     }
     //----- UH ------

     if(domain->solInfo().freqSweepMethod == SolverInfo::Pade ||
        domain->solInfo().freqSweepMethod == SolverInfo::Fourier) {
       double hi = wc-w0;
       h[count] = hi;
     }
     count++;

     if(!(domain->solInfo().freqSweepMethod == SolverInfo::PadeLanczos && count == padeN && domain->coarse_frequencies->size() > 1)) // PJSA 10-09-08 temporary fix to get sweep working
       domain->coarse_frequencies->pop_front();

     if(count == padeN) {  // this means there are enough coarse solves to solve for pade polynomial coefs & reconstruct
       filePrint(stderr, " --------------------------------------\n");
       count = 0;
       double max_freq_before_rebuild = wc;
       int ncoarse = domain->coarse_frequencies->size();
       if((ncoarse > 0) && (padeN==1)) 
         max_freq_before_rebuild = wc + (domain->coarse_frequencies->front() - wc)/2.0;

       // now loop over all the freqencies & reconstruct solutions for each
       // starting with 2nd (1st has already been solved & output)
       StaticTimers *times = probDesc->getStaticTimers();
       while((domain->frequencies->size() > 0) && ((domain->frequencies->front() < max_freq_before_rebuild) 
             || (domain->coarse_frequencies->size() == 0))) {
         startTimerMemory(times->timeFreqSweep, times->memoryFreqSweep);
         domain->isCoarseGridSolve = false;
         double w = domain->frequencies->front();
         double deltaw = w - w0;
         if(domain->solInfo().isAcousticHelm()) {
           filePrint(stderr, " ... Reconstructing solution for k = %f...\n", w/domain->fluidCelerity);
         }
         else {
           filePrint(stderr, " ... Reconstructing solution for f = %f ...\n", w/(2.0*PI));
         }
         //--------
         switch(domain->solInfo().freqSweepMethod) {
           case SolverInfo::Taylor : 
           {
             double _IsConverged = 1.0e-12;
             *sol = *sol_prev[1];
             double prevSolSqNorm, solSqNorm = sol->sqNorm();
             double relError = 0.0;
             //cerr << "  i = 0, sol->sqNorm() = " << sol->sqNorm() << endl;
             for(int i=1; i<nRHS; ++i) {
               prevSolSqNorm = solSqNorm;
               sol->linAdd(1.0/DFactorial(i)*pow(deltaw,i), *sol_prev[i+1]);
               solSqNorm = sol->sqNorm();
               relError = fabs((solSqNorm-prevSolSqNorm)/prevSolSqNorm);
               //cerr << "  i = " << i << ", solSqNorm = " << solSqNorm
               //     << ", factor = " << 1.0/double(DFactorial(i))*pow(deltaw,i) 
               //     << ", sol_prev[i+1]->sqNorm() = " << sol_prev[i+1]->sqNorm() 
               //     << ", relError = " << relError << endl;
               if(relError < _IsConverged) break;
             }
           } break;
           case SolverInfo::Pade :
             if(threadManager->numThr() > 1 && domain->solInfo().type == 2)
               probDesc->pade(sol, sol_prev, h, deltaw);
             else 
               pade(sol, sol_prev, h, deltaw);
             break;
           case SolverInfo::Pade1 :
             pade1(sol, sol_prev, deltaw);
             break;
           //--- UH --- Evaluate the Pade approximation
           case SolverInfo::PadeLanczos :
             /*probDesc->*/ PadeLanczos_Evaluate(padeN * nRHS, PadeLanczos_solprev, 
                                            PadeLanczos_VtKV, PadeLanczos_Vtb, w, sol);
             break;
           //--- UH ---
           case SolverInfo::Fourier :
             fourier(sol, sol_prev, h, deltaw);
             break;
         } // switch (domain->solInfo().freqSweepMethod)
         //--------
         geoSource->setImpe(w); // for output file and/or restart
         bool printTimers = ((ncoarse+domain->frequencies->size()) > 1) ? false : true;
         stopTimerMemory(times->timeFreqSweep, times->memoryFreqSweep);
         postProcessor->staticOutput(*sol, *rhs, printTimers); 
         domain->frequencies->pop_front();
       } // while((domain->frequencies->size() > 0) ... )

       if(domain->solInfo().freqSweepMethod == SolverInfo::PadeLanczos) { // PJSA 10-09-08 temporary fix to get sweep working
                                                                          // not optimal, can we reuse part of the basis? 
         first_time = true; // so the solver isn't rebuild
       }
       else {
         if((padeN > 1) && (ncoarse > 0)) {
           count = (padeN-ncoarse > 1) ? padeN-ncoarse : 1;
           int coffset = (padeN-count)*(nRHS+1);
           for(int i = 0; i < (nRHS+1)*count; ++i) *sol_prev[i] = *sol_prev[coffset+i];
           double deltah = h[padeN-count];
           w0 += deltah;
           for(int i=0; i<count; ++i) h[i] = h[padeN-count+i] - deltah;
         } // if((padeN > 1) && (ncoarse > 0))
       }

     } // if (count == padeN)

   } // while(domain->coarse_frequencies->size() > 0)
   for(int i = 0; i < (nRHS+1)*padeN; ++i)
     delete sol_prev[i];
   delete [] sol_prev;
   delete [] h;
   //--- UH --- Data for Pade Lanczos
   if (PadeLanczos_VtKV)
     delete[] PadeLanczos_VtKV;
   if (PadeLanczos_Vtb)
     delete[] PadeLanczos_Vtb;
   if (PadeLanczos_solprev)
     delete[] PadeLanczos_solprev;
   //--- UH --- Data for Pade Lanczos
 } // if(domain->solInfo().doFreqSweep)
 else if(domain->solInfo().probType == SolverInfo::HelmholtzDirSweep) { // PJSA: never tested this in FEM/DPH
   int iRHS;
   int nDir = domain->getNumWaveDirections();
   if (nDir==0) nDir =1;
   for (iRHS = 0; iRHS < nDir; iRHS++) {
     probDesc->setIWaveDir(iRHS);
     if(nDir > 1) filePrint(stderr,"\n Running RHS #%d\n",iRHS);
     probDesc->getRHS(*rhs);
     allOps->sysSolver->solve(*rhs, *sol);
     postProcessor->staticOutput(*sol, *rhs);
     rhs->zero();
     sol->zero();
   }
 }
 else if(domain->solInfo().noninpc) {
   // initialization
   probDesc->getRHS(*rhs);
   VecType *psi_u = new VecType(probDesc->solVecInfo(sfem->getP()));
   if(domain->solInfo().type != 2) psi_u->setn(domain->numUncon());

   psi_u->zero();
   for(int i=0; i<domain->solInfo().nsample; ++i) {
     cout << endl << "Realization number  " << i+1 << endl;
     if(i>0) {
       sfem->genXiPsi(i); // seed=i 
       probDesc->assignRandMat();
       probDesc->rebuildSolver(); 
       //solver = probDesc->getSolver();
     }
     allOps->sysSolver->solve(*rhs,*sol);
     sfem_noninpc->update_Psi_u(sol,psi_u);
     sol->zero();
   }
   delete sol;
   sol = 0;
   sfem_noninpc->compute_Coefs(psi_u); 
   VecType *sol1 = new VecType(probDesc->solVecInfo(1));
   if(domain->solInfo().type != 2) sol1->setn(domain->numUncon());

   sfem_noninpc->computeMean(psi_u,sol1);
   postProcessor->staticOutput(*sol1, *rhs, false, 1);
   delete sol1;
   VecType *sol2 = new VecType(probDesc->solVecInfo(1));
   if(domain->solInfo().type != 2) sol2->setn(domain->numUncon());

   sfem_noninpc->computeStdDev(psi_u,sol2);
   postProcessor->staticOutput(*sol2, *rhs, false, 2);
   delete sol2;

   VecType *sol3 = new VecType(probDesc->solVecInfo(1));
   if(domain->solInfo().type != 2) sol3->setn(domain->numUncon());
   int nsample_output = sfem->getnsamp_out();
   for(int i=0; i< nsample_output; ++i) {
     sfem_noninpc->computePdf(i,psi_u,sol3);
     postProcessor->staticOutput(*sol3, *rhs, true, 3);
   }
   delete sol3;
   stochStress(psi_u); 
 }
 else if(domain->solInfo().inpc) {
   probDesc->getRHSinpc(*rhs);
   allOps->sysSolver->solve(*rhs,*sol);
//   postProcessor->staticOutput(*sol, *rhs, true, 0); /

   VecType *sol1 = new VecType(probDesc->solVecInfo(1));
   if(domain->solInfo().type != 2) sol1->setn(domain->numUncon());
   sfem_inpc->computeMean(sol,sol1);
   postProcessor->staticOutput(*sol1, *rhs, true, 1);
   delete sol1;
 
   VecType *sol2 = new VecType(probDesc->solVecInfo(1)); 
   if(domain->solInfo().type != 2) sol2->setn(domain->numUncon());
   sfem_inpc->computeStdDev(sol,sol2);
   postProcessor->staticOutput(*sol2, *rhs, false, 2);
   delete sol2;

   VecType *sol3 = new VecType(probDesc->solVecInfo(1));
   if(domain->solInfo().type != 2) sol3->setn(domain->numUncon());
   int nsample_output = sfem->getnsamp_out();
   for(int i=0; i<nsample_output; ++i) {
     sfem_inpc->computePdf(i,sol,sol3); 
     postProcessor->staticOutput(*sol3, *rhs, false, 3);
   }
   delete sol3;

   probDesc->retrieveElemset(); 
   delete allOps->sysSolver; 
   allOps->sysSolver = 0;
   stochStress(sol); // YYY Commented out temporarily
 }
 else {
   probDesc->getRHS(*rhs);
   allOps->sysSolver->solve(*rhs,*sol);
   if(domain->solInfo().isCoupled) scaleDisp(*sol); // PJSA 9-22-06
   postProcessor->staticOutput(*sol, *rhs);
 }
 geoSource->closeOutputFiles(); 
}

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::stochStress(VecType *sol)
{
 
//  int n = domain->numUncon();
//  int P = sfem->getP();
//  int n_node = geoSource->numNode();
//  int npsize=P*n_node;

  MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>* mltintg = new MultInteg<Scalar, VecType, PostProcessor, ProblemDescriptor>();


//  int numNodes = geoSource->numNode();  // PJSA 8-26-04 don't want to print displacements for internal nodes
  Scalar *globVal = 0;
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  int qmd=4;
  int dof;
//  int iNode;
  for(int i = 0; i < numOutInfo; ++i)  { 
    if(oinfo[i].ndtype == 0) continue; // YYY fix ndflag issue
    if(oinfo[i].interval == 1) {
      dof = -1;

/*      if(domain->solInfo().noninpc)
        mltintg->computeStressStat(qmd,sol,i,oinfo[i].type,postProcessor,probDesc,oinfo[i].ndtype);
      else
        mltintg->computeStressStat(qmd,sol,i,oinfo[i].type,postProcessor,probDesc,oinfo[i].ndtype);
*/
      switch(oinfo[i].type)  {
     
        default:
        //  fprintf(stderr, " *** WARNING: output %d is not supported \n", i);
          break;

        case OutputInfo::StressXX:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SXX,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SXX,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StressYY:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SYY,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SYY,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StressZZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SZZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SZZ,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StressXY:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SXY,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SXY,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StressYZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SYZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SYZ,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StressXZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,SXZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,SXZ,postProcessor,probDesc,oinfo[i].ndtype);
          break; 
        case OutputInfo::StressVM:
            if(domain->solInfo().noninpc)
             mltintg->computeStressStat(qmd,sol,i,VON,postProcessor,probDesc,oinfo[i].ndtype);
            else
             mltintg->computeStressStat(qmd,sol,i,VON,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainXX:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EXX,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EXX,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainYY:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EYY,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EYY,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainZZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EZZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EZZ,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainXY:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EXY,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EXY,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainYZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EYZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EYZ,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainXZ:
           if(domain->solInfo().noninpc)
            mltintg->computeStressStat(qmd,sol,i,EXZ,postProcessor,probDesc,oinfo[i].ndtype);
           else
            mltintg->computeStressStat(qmd,sol,i,EXZ,postProcessor,probDesc,oinfo[i].ndtype);
          break;
        case OutputInfo::StrainVM:
            if(domain->solInfo().noninpc)
             mltintg->computeStressStat(qmd,sol,i,STRAINVON,postProcessor,probDesc,oinfo[i].ndtype);
            else
             mltintg->computeStressStat(qmd,sol,i,STRAINVON,postProcessor,probDesc,oinfo[i].ndtype);
          break;

      }
    } 
    if(globVal)
      { delete [] globVal; globVal = 0; }
  } 
  delete mltintg;
}

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::rebuildSolver(double w1)
{
  dbg_alloca(0);
  geoSource->setOmega(w1);
  if(domain->solInfo().isAcousticHelm()) {  // set new wave number for acoustic helmholtz
    filePrint(stderr,"\n ... Rebuilding LHS with wave number = %f ... \n", w1/domain->fluidCelerity);
    SPropContainer& sProps = geoSource->getStructProps();
    for(int iProp=0;iProp<geoSource->getNumProps();iProp++) {
       if(sProps[iProp].kappaHelm!=0.0 || sProps[iProp].kappaHelmImag!=0.0) { // Fluid
         complex<double> k1 = w1/sProps[iProp].soundSpeed;
         sProps[iProp].kappaHelm = real(k1);
         sProps[iProp].kappaHelmImag = imag(k1);
       } 
    }
  }
  else {
    filePrint(stderr,"\n ... Rebuilding LHS with frequency = %f ... \n", w1/(2.0*PI));
  }
  probDesc->rebuildSolver();
  //solver = probDesc->getSolver();

  if(a) {
    for(int i = 0; i < ia; ++i) delete a[i];
    delete [] a;
    a = 0;
  }
  if(b) {
    for(int i = 0; i < ib; ++i) delete b[i];
    delete [] b;
    b = 0;
  }
  if(P) { delete P; P = 0; }
  if(Q) { delete Q; Q = 0; }

  if(ca) {
    for(int i = 0; i < ia; ++i) delete ca[i];
    delete [] ca;
    ca = 0;
  }
  if(cb) {
    for(int i = 0; i < ib; ++i) delete cb[i];
    delete [] cb;
    cb = 0;
  }
  if(cP) { delete cP; cP = 0; }
  if(cQ) { delete cQ; cQ = 0; }

}

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::scaleDisp(VecType &u)
{
  probDesc->scaleDisp(u);
}


//------------------------------------------------------------------------------

