#include <Utils.d/dbg_alloca.h>
#include <Threads.d/Paral.h>
#include <Threads.d/PHelper.h>
#include <Math.d/mathUtility.h>
#include <Driver.d/Domain.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/DBSparseMatrix.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/Memory.h>
#include <Driver.d/GeoSource.h>

#ifdef DISTRIBUTED
#include <Dist.d/DistDom.h>
#endif

template<class Scalar>
GenMultiDomainStatic<Scalar>::GenMultiDomainStatic(Domain *d)
{
 domain = d;
/*
 switch(domain->solInfo().fetiInfo.version) {
  default:
  case FetiInfo::feti1:
    filePrint(stderr," ... FETI-1 is Selected             ...\n");
    break;
  case FetiInfo::feti2:
    filePrint(stderr," ... FETI-2 is Selected             ...\n");
    break;
  case FetiInfo::fetidp:
    if (!(domain->solInfo().fetiInfo.dph_flag)) 
      filePrint(stderr," ... FETI-Dual/Primal is Selected   ...\n");
    else 
      filePrint(stderr," ... FETI-DPH is Selected           ...\n");
    break;
  }
*/
#ifdef DISTRIBUTED
 decDomain = new GenDistrDomain<Scalar>(domain);
#else
 decDomain = new GenDecDomain<Scalar>(domain);
#endif

 times  = new StaticTimers;

 // debug: check static timers are initialized
 //MatrixTimers &mt = domain->getTimers();
}

template<class Scalar>
DistrInfo &
GenMultiDomainStatic<Scalar>::solVecInfo()
{
 return solver->localInfo();
}

template<class Scalar>
DistrInfo &
GenMultiDomainStatic<Scalar>::solVecInfo(int i) 
{
 cerr << "Warning : GenMultiDomainStatic<Scalar>::solVecInfo(int i) should not be called" << endl;
 return solver->localInfo();
}


template<class Scalar>
void
GenMultiDomainStatic<Scalar>::getRHS(GenDistrVector<Scalar> &rhs)
{
 times->formRhs -= getTime();
 execParal1R(decDomain->getNumSub(), this, &GenMultiDomainStatic<Scalar>::subGetRHS, rhs);

 // eigen mode projector 
 if(domain->solInfo().modeFilterFlag)
   project(rhs);

 times->formRhs += getTime();
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::subGetRHS(int isub, GenDistrVector<Scalar>& rhs)
{
 GenSubDomain<Scalar> *sd = decDomain->getSubDomain(isub);
 GenStackVector<Scalar> subrhs(rhs.subData(isub), rhs.subLen(isub));
 sd->buildRHSForce(subrhs, (*allOps.Kuc)[isub]);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::preProcess()
{
 // Makes renumbering, connectivities and dofsets
 startTimerMemory(times->preProcess, times->memoryPreProcess);
 decDomain->preProcess();
 stopTimerMemory(times->preProcess, times->memoryPreProcess);

 // Construct FETI solver and any other matrices required
 times->getFetiSolverTime -= getTime();

 if(domain->solInfo().freqSweepMethod == SolverInfo::PadeLanczos || 
    domain->solInfo().freqSweepMethod == SolverInfo::GalProjection) {
   for(int i=0; i<decDomain->getNumSub(); ++i) decDomain->getSubDomain(i)->makeAllDOFs();
   decDomain->buildOps(allOps, 0.0, 0.0, 1.0);
   solver = (GenFetiSolver<Scalar> *) allOps.sysSolver;
 }
 else {
   allOps.sysSolver = solver = decDomain->getFetiSolver();
   allOps.Kuc = new GenSubDOp<Scalar>(decDomain->getNumSub());
   for(int i=0; i<decDomain->getNumSub(); ++i) (*allOps.Kuc)[i] = decDomain->getSubDomain(i)->Kuc;
 }

 times->getFetiSolverTime += getTime();

 int useModeFilter = domain->solInfo().modeFilterFlag;
 if(useModeFilter) {
   filePrint(stderr, " ... MODEfilter requested          ...\n");
   eigmode_projector_prep();
 }

}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::rebuildSolver()
{
  times->getFetiSolverTime -= getTime(); // PJSA 3-30-06
  if(domain->solInfo().freqSweepMethod == SolverInfo::PadeLanczos ||
     domain->solInfo().freqSweepMethod == SolverInfo::GalProjection) {
    GenMDDynamMat<Scalar> ops;
    ops.sysSolver = allOps.sysSolver;
    ops.K = allOps.K;
    ops.Kuc = allOps.Kuc;
    decDomain->rebuildOps(ops, 0.0, 0.0, 1.0); // just rebuild solver and K
  }
  else {
    solver->rebuildLHSfreq();
  }
  times->getFetiSolverTime += getTime(); // PJSA 3-30-06
}


template<class Scalar>
void
GenMultiDomainStatic<Scalar>::scaleDisp(GenDistrVector<Scalar> &u)
{
 decDomain->scaleDisp(u);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::clean()
{
}

template<class Scalar>
GenFetiSolver<Scalar> *
GenMultiDomainStatic<Scalar>::getSolver()
{
 return solver;
}

template<class Scalar>
GenMultiDomainPostProcessor<Scalar> *
GenMultiDomainStatic<Scalar>::getPostProcessor()
{
 return new GenMultiDomainPostProcessor<Scalar>(decDomain, solver, times);
}

template<class Scalar>
void
GenMultiDomainPostProcessor<Scalar>::staticOutput(GenDistrVector<Scalar> &sol, GenDistrVector<Scalar> &force,
                                                  bool printTimers, int ndflag)
{
 startTimerMemory(times->output, times->memoryOutput);
 decDomain->postProcessing(sol, force);
 stopTimerMemory(times->output, times->memoryOutput);

 times->numSubdomain = decDomain->getNumSub();

 // --- Memory computations should be done so that other problems using
 // --- FETI can use the same routines.

 // --- NOTE: the average, min, and max should be done per processor,
 // ---       not per subdomain.

 int i;
 long (*memory)=(long *) dbg_alloca(sizeof(long)*decDomain->getNumSub());
 long totMemPrec = 0, totMemK = 0;

 // Kii memory calculations
 for(i=0; i<decDomain->getNumSub(); ++i) memory[i] = 0;
 execParal(decDomain->getNumSub(), this,  
           &GenMultiDomainPostProcessor<Scalar>::getMemoryPrec, memory);
 for(i=0; i<decDomain->getNumSub(); ++i) totMemPrec += memory[i];

 // K memory calculations
 for(i=0; i<decDomain->getNumSub(); ++i) memory[i] = 0;
 execParal(decDomain->getNumSub(), this, 
           &GenMultiDomainPostProcessor<Scalar>::getMemoryK, memory);
 for(i=0; i<decDomain->getNumSub(); ++i) totMemK += memory[i];

 Timings &fetiTimers = solver->getTimers();
 fetiTimers.preconditioner.addOverAll(totMemPrec, 0.0);
 fetiTimers.kMem.addOverAll(totMemK, 0.0);

#ifdef DISTRIBUTED
 double mem1 = (double) totMemPrec;
 if(structCom) mem1 = structCom->globalSum(mem1);
 totMemPrec = (long) mem1;

 mem1 = (double) totMemK;
 if(structCom) mem1 = structCom->globalSum(mem1);
 totMemK = (long) mem1;
#endif

 times->memoryPrecond = totMemPrec;
 times->memoryK       = totMemK;

 if(printTimers) {
   filePrint(stderr," ... Print Timers                   ... \n");
   switch(domain->solInfo().fetiInfo.version) {
     default:
     case FetiInfo::feti1:
     case FetiInfo::feti2:
       times->printStaticTimers(domain->getTimers(), 
                                solver->getSolutionTime(), 
                                domain->solInfo() , 
                                solver->getTimers(),
                                geoSource->getCheckFileInfo()[0],
                                domain);
       break;

     case FetiInfo::fetidp:
       times->printFetiDPtimers(domain->getTimers(),
                                solver->getSolutionTime(),
                                domain->solInfo() ,
                                solver->getTimers(),
                                geoSource->getCheckFileInfo()[0],
                                domain);
       break;
   }
//   geoSource->closeOutputFiles();
 }

 filePrint(stderr," --------------------------------------\n");
}

template<class Scalar>
void
GenMultiDomainPostProcessor<Scalar>::getMemoryK(int iSub, long *memory)
{
 memory[iSub] = decDomain->getSubDomain(iSub)->getMemoryK();
}

template<class Scalar>
void
GenMultiDomainPostProcessor<Scalar>::getMemoryPrec(int iSub, long *memory)
{
 memory[iSub] = decDomain->getSubDomain(iSub)->getMemoryPrec();
}

// FETI-H
template<class Scalar>
void
GenMultiDomainStatic<Scalar>::setIWaveDir(int _i)
{
 int i;
 for(i=0;i<decDomain->getNumSub();i++)
   decDomain->getSubDomain(i)->iWaveDir = _i;
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::getFreqSweepRHS(GenDistrVector<Scalar> *rhs, 
                                              GenDistrVector<Scalar> **sol_prev, int iRHS)
{
 if(domain->solInfo().freqSweepMethod == SolverInfo::PadeLanczos || 
    domain->solInfo().freqSweepMethod == SolverInfo::GalProjection) {
   for(int i=0;i<decDomain->getNumSub();i++) {
     decDomain->getSubDomain(i)->M = (*allOps.M)[i];
     decDomain->getSubDomain(i)->Muc = (GenCuCSparse<Scalar> *)(*allOps.Muc)[i];
     if (allOps.C) decDomain->getSubDomain(i)->C = (*allOps.C)[i];
     if (allOps.C_deriv) {
       decDomain->getSubDomain(i)->C_deriv =
          new GenSparseMatrix<Scalar>*[iRHS+1];
       (decDomain->getSubDomain(i)->C_deriv)[0] = ((*allOps.C_deriv)[0])[i];
       for(int j=1;j<iRHS+1;j++) (decDomain->getSubDomain(i)->C_deriv)[j] = 0;
     }
     if (allOps.Cuc_deriv) {
       decDomain->getSubDomain(i)->Cuc_deriv =
          new GenSparseMatrix<Scalar>*[iRHS+1];
       (decDomain->getSubDomain(i)->Cuc_deriv)[0] = ((*allOps.Cuc_deriv)[0])[i];
       for(int j=1;j<iRHS+1;j++) (decDomain->getSubDomain(i)->Cuc_deriv)[j] = 0;
     }
   }
 }
 solver->getFreqSweepRHS(rhs, sol_prev, iRHS);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::getRHS(GenDistrVector<Scalar> &rhs, double omega,
                                     double deltaomega)
{
 solver->makeStaticLoad(rhs,omega,deltaomega);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::pade(GenDistrVector<Scalar> *sol, GenDistrVector<Scalar> **sol_prev, double *h, double x)
{
  solver->pade(sol, sol_prev, h, x);
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::eigmode_projector_prep()
{
  //if(Rmem) return; // already done it or requested some other filter

  // Read computed eigenvectors from file EIGENMODES
  // ======================================
#ifdef DISTRIBUTED
  char *filename = new char[40];
  sprintf(filename,"EIGENMODES%d",structCom->myID());
  BinFileHandler modefile(filename, "r");
#else
  BinFileHandler modefile("EIGENMODES", "r");
#endif
  if(modefile.get_fileid() <= 0) { fprintf(stderr, " *** Error: Failed to open EIGENMODES file ***\n"); exit(-1); }

  modefile.read(&numR, 1);
  filePrint(stderr," ... Reading %d modes from EIGENMODES file ...\n", numR);

  int eigsize;
  modefile.read(&eigsize, 1);
  if(eigsize != solVecInfo().totLen()) { fprintf(stderr, " *** Error: Bad data in EIGENMODES file %d %d ***\n", eigsize, solVecInfo().totLen()); exit(-1); }

  Rmem = new GenDistrVector<Scalar> * [numR];
  double *data = new double[eigsize];
  for(int i = 0; i < numR; ++i) {
    Rmem[i] = new GenDistrVector<Scalar>(solVecInfo());
    modefile.read(data, eigsize);
    for(int j=0; j<eigsize; ++j) Rmem[i]->data()[j] = data[j];
  }
  delete [] data;

  // Build (U_c^T*U_c)^{-1} for projector
  // ======================================
  GenFSFullMatrix<Scalar> RtR(numR, numR);
  for(int i=0; i<numR; ++i) {
    Scalar rii = (*(Rmem[i]))*(*(Rmem[i]));
    RtR[i][i] = rii;
    for(int j=0; j<i; ++j) {
      Scalar rij = (*(Rmem[i]))*(*(Rmem[j]));
      RtR[i][j] = rij;
      RtR[j][i] = rij;
    }
  }
  RtRinverse = new GenFSFullMatrix<Scalar>(numR, numR);
  (*RtRinverse) = RtR.invert();
}

template<class Scalar>
void
GenMultiDomainStatic<Scalar>::project(GenDistrVector<Scalar> &b)
{
  Scalar *Rtb = new Scalar[numR];
  for(int i=0; i<numR; ++i) Rtb[i] = (*Rmem[i]) * b;
  Scalar *x = new Scalar[numR];
  RtRinverse->mult(Rtb, x);
  for(int i=0; i<b.size(); ++i) 
    for(int j=0; j<numR; ++j) b.data()[i] -= Rmem[j]->data()[i]*x[j];

  ///*check:*/ for(int i=0; i<numR; ++i) cerr << (*Rmem[i]) * b << " "; cerr << endl;
}

