#include <stdio.h>

#include <Threads.d/Paral.h>
#include <Driver.d/Domain.h>
#include <Paral.d/MDNLStatic.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Driver.d/SubDomain.h>
#include <Threads.d/PHelper.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Math.d/Vector.h>
#include <Math.d/mathUtility.h>
#include <Utils.d/DistHelper.h>
#include <Paral.d/MDStatic.h>
#ifdef DISTRIBUTED
#include <Dist.d/DistDom.h>
#endif

//#define PRINT_NLTIMERS

void
MDNLStatic::getSubStiffAndForce(int isub, DistrGeomState &geomState, 
                                DistrVector &res, DistrVector &elemIntForce)
{
 SubDomain *sd = decDomain->getSubDomain(isub);

 StackVector residual(res.subData(isub), res.subLen(isub));

 // eIF = element internal force
 StackVector eIF(elemIntForce.subData(isub), elemIntForce.subLen(isub));

 sd->getStiffAndForce(*geomState[isub], eIF, allCorot[isub], kelArray[isub],
                      residual);

 // PJSA: start LMPC code
 // filePrint(stderr, " ... Processing LMPCs for non-linear FETI ...\n");
 sd->updateMpcRhs(*geomState[isub], decDomain->getMpcToSub());
 double *mpcForces = new double[sd->numMPCs()]; // don't delete  
 solver->getLocalMpcForces(isub, mpcForces);  // mpcForces set to incremental mpc lagrange multipliers
 sd->addMpcForceIncrement(mpcForces);  // mpcForces set to total mpc lagrange multipliers
 // cerr << "mpcForces = "; for(int i=0; i<sd->numMPCs(); ++i) cerr << mpcForces[i] << " "; cerr << endl;
 sd->constraintProductTmp(mpcForces, residual); // C^T*lambda added to force residual
 // PJSA: end LMPC code
}

void
MDNLStatic::makeSubKelArrays(int isub)
{
 SubDomain *sd = decDomain->getSubDomain(isub);
 sd->createKelArray(kelArray[isub]);
}

void
MDNLStatic::makeSubCorotators(int isub)
{
 SubDomain *sd  = decDomain->getSubDomain(isub);
 int numele     = sd->numElements();
 allCorot[isub] = new Corotator*[numele];
 sd->createCorotators(allCorot[isub]);
}

// Constructor
MDNLStatic::MDNLStatic(Domain *d)
{
 domain = d;
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
#ifdef DISTRIBUTED
 decDomain = new GenDistrDomain<double>(domain);
#else
 decDomain = new GenDecDomain<double>(domain);
#endif
 numSystems = 0;
}

DistrInfo&
MDNLStatic::solVecInfo()
{
 return solver->localInfo();
}

DistrInfo&
MDNLStatic::sysVecInfo()
{
 return solver->localInfo();
}


DistrInfo&
MDNLStatic::elemVecInfo()
{
 return *decDomain->elementVectorInfo();
}

int
MDNLStatic::checkConvergence(int iter, double normDv, double normRes)
{
 times->timeCheck -= getTime();

 if(iter == 0) { 
   firstDv  = normDv;
   firstRes = normRes;
   lastRes  = 0.0;
 }

 double relativeDv  = normDv/firstDv;
 double relativeRes = normRes/firstRes;
 lastRes = normRes;

 filePrint(stderr,"----------------------------------------------------\n");
 filePrint(stderr,"Nonlinear Iter #%d\tcurrent dv   = % e\n \t\t\t"
                  "first dv     = % e\n \t\t\trelative dv  = % e\n",
                   iter+1, normDv, firstDv, relativeDv);
 filePrint(stderr,"                \tcurrent Res  = % e\n \t\t\t"
                  "first Res    = % e\n \t\t\trelative Res = % e\n",
                   normRes, firstRes, relativeRes);

 filePrint(stderr,"----------------------------------------------------\n");
 
 int converged = 0;

 // KHP: Charbel requested convergence be monitored based on residual only.
 // Check incremental displacement
 //if(relativeDv <= domain->solInfo().getNLInfo().tolRes)
 //  converged = 1;

 // Check to see if residual has converged
 if(relativeRes <= domain->solInfo().getNLInfo().tolRes)
  converged = 1;

 // Divergence check
 // if( normDv > 1000.0*firstDv || normRes > 1000.0*firstRes)
 if( normDv > 1000.0*firstDv)
   converged = -1;

 // Store residual norm and dv norm for output
 times->norms[numSystems].normDv      = normDv;
 times->norms[numSystems].relativeDv  = relativeDv;
 times->norms[numSystems].normRes     = normRes;
 times->norms[numSystems].relativeRes = relativeRes;

 numSystems += 1;

 times->timeCheck += getTime();

 return converged;
}

double
MDNLStatic::getStiffAndForce(DistrGeomState& geomState, 
                             DistrVector& residual, DistrVector& elementInternalForce,
                             DistrVector&)
{
 times->buildStiffAndForce -= getTime();

 execParal3R(decDomain->getNumSub(),this,&MDNLStatic::getSubStiffAndForce,geomState,
            residual,elementInternalForce);

 times->buildStiffAndForce += getTime();
 
 return sqrt(solver->getFNormSq(residual));
}


DistrGeomState*
MDNLStatic::createGeomState()
{
 times->timeGeom -= getTime();
 
 DistrGeomState* geomState = new DistrGeomState(decDomain);
 
 times->timeGeom += getTime();
#ifdef PRINT_NLTIMERS 
 filePrint(stderr," ... Time to Make Subdomains Geometry States %18e\n",
           times->timeGeom/1000.0);
#endif
 return geomState;
}

void
MDNLStatic::updatePrescribedDisplacement(DistrGeomState *geomState, double)
{
 times->timePresc -= getTime();

 execParal1R(decDomain->getNumSub(),this,&MDNLStatic::updatePrescribedDisp,*geomState);

 times->timePresc += getTime();
#ifdef PRINT_NLTIMERS
 filePrint(stderr," ... Time to Update Prescribed Displacements %18e\n",
           times->timePresc/1000.0);
#endif
}

void
MDNLStatic::updatePrescribedDisp(int isub, DistrGeomState& geomState)
{
 SubDomain *sd = decDomain->getSubDomain(isub);
 sd->updatePrescribedDisp(geomState[isub], deltaLambda);
} 

int
MDNLStatic::reBuild(int iteration, int step, DistrGeomState& geomState )
{
/*
 times->rebuild -= getTime();

 int rebuildFlag = 0;
 // Always Rebuild at the beginning of a step, then rebuild every N iterations
 // within a step

 // rebuild every certain number of newton iterations or 
 // at the beginning of a load step
 if( iteration % domain->solInfo().getNLInfo().updateK == 0 || 
     iteration == 0)
 // rebuild until a certain iteration
 //if( iteration < domain->solInfo().getNLInfo()updateK  || 
 //     iteration == 0 ||
 //    numSystems == 1) 
 {
   filePrint(stderr,"===> REBUILDING TANGENT STIFFNESS MATRIX\n");
   times->norms[numSystems].rebuildTang = 1;
   solver->reBuild(kelArray, geomState, iteration, step);
   rebuildFlag = 1;
 } else
   times->norms[numSystems].rebuildTang = 0;

// This part is necessary if you decide not to rebuild the preconditioner
// we found rebuilding the preconditioner is benefitial in cpu time so we
// decided to always rebuild the preconditioner if the tangent stiffness
// matrix was being rebuilt.
// } else {
//   fprintf(stderr,"===> REBUILDING ERROR ESTIMATOR\n");
//   solver->reBuildErrorEstimator(kelArray);
// }
 times->rebuild += getTime();
 
 return rebuildFlag;
*/
 times->rebuild -= getTime();

 int rebuildFlag = 0;

 if (iteration % domain->solInfo().getNLInfo().updateK == 0) {
   GenMDDynamMat<double> allOps;
   allOps.sysSolver = solver;
   decDomain->rebuildOps(allOps, 0.0, 0.0, 1.0, kelArray);
   rebuildFlag = 1;
 }

 times->rebuild += getTime();

 return rebuildFlag;

}

void
MDNLStatic::makeSubDofs(int isub)
{
 SubDomain *sd = decDomain->getSubDomain(isub);
 sd->makeAllDOFs();
}

void
MDNLStatic::preProcess()
{
 // Structure to store timers
 times = new StaticTimers;

 times->memoryPreProcess -= threadManager->memoryUsed();
 
 // Constructs renumbering, connectivities and dofsets
 times->preProcess -= getTime();
 decDomain->preProcess();
 times->preProcess += getTime();

 // Get number of subdomains
 int numSub = decDomain->getNumSub();

 // Make subdomain's degrees of freedom
 times->makeDOFs -= getTime();
 execParal(numSub, this, &MDNLStatic::makeSubDofs);
 times->makeDOFs += getTime();
#ifdef PRINT_NLTIMERS
 filePrint(stderr," ... Time to Make Subdomain Dofs %30e\n",
         times->makeDOFs/1000.0);
#endif 
 // Make subdomain's corotators
 times->corotatorTime -= getTime();
 allCorot = new Corotator**[numSub]; 
 execParal(numSub, this, &MDNLStatic::makeSubCorotators);
 times->corotatorTime += getTime();
#ifdef PRINT_NLTIMERS
 filePrint(stderr," ... Time to Make Subdomain Corotators %24e\n", 
         times->corotatorTime/1000.0);
#endif

 times->memoryPreProcess += threadManager->memoryUsed();

 // NOTE: count FETI memory separately from pre-process memory in
 // timing file

 // Construct FETI Solver, build GtG, subdomain matrices, etc.
 times->getFetiSolverTime -= getTime();
 //solver = decDomain->getFetiSolver();
 GenMDDynamMat<double> allOps;

 double Kcoef = 1.0;
 double Mcoef = 0.0;
 double Ccoef = 0.0;

 filePrint(stderr," ... Building Solver                ...\n");
/*
 Rbm *rigidBodyModes = 0;
 if(domain->solInfo().rbmflg) rigidBodyModes = domain->constructRbm(); // new policy is to construct rbms if GRBM is requested in input file
                                                                       // but only use them when it is appropriate to do so. In nonlinear statics it is not
                                                                       // since the nullity of the tangent stiffness matrix may be less than the nullity
                                                                       // of the number of rigid body modes
*/
 decDomain->buildOps(allOps, Mcoef, Ccoef, Kcoef, (Rbm **) 0);
 solver = (GenFetiSolver<double> *) allOps.sysSolver;
 times->getFetiSolverTime += getTime();
#ifdef PRINT_NLTIMERS
 filePrint(stderr," ... Time to Construct FETI Solver %28e\n", times->getFetiSolverTime/1000.0);
#endif

 // Make subdomain's array of stiffness matrices
 times->memoryPreProcess -= threadManager->memoryUsed();
 times->kelArrayTime -= getTime();
 kelArray = new FullSquareMatrix*[numSub];
 execParal(numSub, this, &MDNLStatic::makeSubKelArrays);
 times->kelArrayTime += getTime();
 times->memoryPreProcess += threadManager->memoryUsed();
#ifdef PRINT_NLTIMERS
 filePrint(stderr," ... Time to Make Subdomain Stiffness Arrays %18e\n", 
        times->kelArrayTime/1000.0);
#endif
 tolerance = domain->solInfo().getNLInfo().tolRes;
}

int
MDNLStatic::getMaxit()
{
 return domain->solInfo().getNLInfo().maxiter;
}

// Just for defining a minimum and maximum delta Lambda
double
MDNLStatic::getScaleFactor()
{
 return domain->solInfo().getNLInfo().lfactor;
}

double
MDNLStatic::getDeltaLambda0()
{
 deltaLambda = domain->solInfo().getNLInfo().dlambda;
 return deltaLambda;
}

double
MDNLStatic::getMaxLambda()
{
 return domain->solInfo().getNLInfo().maxLambda;
}

void
MDNLStatic::getRHS(DistrVector& rhs, DistrGeomState *gs)
{
 // ... BUILD THE RHS FORCE (external + gravity + nonhomogeneous)
 times->formRhs -= getTime();

 solver->makeStaticLoad(rhs,gs); //HB: add DistrGeomState for computing 
                                 //    follower forces (i.e. pressure)
 times->formRhs += getTime(); 
#ifdef PRINT_NLTIMERS 
 filePrint(stderr," ... Time to Assemble External Forces %25e\n",
           times->formRhs/1000.0);
#endif
}

FetiSolver *
MDNLStatic::getSolver()
{
  return solver;
}

MultiDomainPostProcessor *
MDNLStatic::getPostProcessor()
{
 return new MultiDomainPostProcessor(decDomain, solver, times);
}

void
MDNLStatic::staticOutput(DistrGeomState *geomState, double lambda,
                         DistrVector &Force, DistrVector &)
{
 startTimerMemory(times->output, times->memoryOutput);

 decDomain->postProcessing(geomState, allCorot, lambda);

 stopTimerMemory(times->output, times->memoryOutput);
#ifdef PRINT_NLTIMERS 
 filePrint(stderr," ... Time to Postprocess %35e\n", times->output/1000.0);
#endif
}

void
MDNLStatic::printTimers()  
{
  int i;
  long (*memory)=(long *) dbg_alloca(sizeof(long)*decDomain->getNumSub());
  for (i = 0; i < decDomain->getNumSub(); ++i)
    memory[i] = 0;

  MultiDomainPostProcessor *mdpp = getPostProcessor();

  execParal(decDomain->getNumSub(), mdpp,
           &MultiDomainPostProcessor::getMemoryPrec, memory);
  long totMemPrec = 0;
  for(i = 0; i < decDomain->getNumSub(); ++i)
    totMemPrec += memory[i];

  for(i = 0; i < decDomain->getNumSub(); ++i)
    memory[i] = 0;

  execParal(decDomain->getNumSub(), mdpp,
           &MultiDomainPostProcessor::getMemoryK, memory);

  long totMemK = 0;
  for(i=0; i<decDomain->getNumSub(); ++i)
    totMemK += memory[i];

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

  times->memoryK = totMemK;
  times->memoryPrecond = totMemPrec;

  times->timeTimers -= getTime();
 
  times->numSubdomain = decDomain->getNumSub();

  times->printTimers(domain, solver->getTimers(),
                     solver->getSolutionTime());
		    
  times->timeTimers += getTime();
#ifdef PRINT_NLTIMERS 
  filePrint(stderr," ... Time to Print Timers %e\n", times->timeTimers/1000.0);
#endif
}

/*
template<class Scalar>
void
GenMultiDomainPostProcessor<Scalar>::staticOutput(DistrGeomState *geomState, double lambda)
{
 times->output -= getTime();
 //decDomain->postProcessing(geomState, deformations, lambda);
 times->output += getTime();
#ifdef PRINT_NLTIMERS
 filePrint(stderr," ... Time to Postprocess %e\n", times->output/1000.0);
#endif
}
*/
