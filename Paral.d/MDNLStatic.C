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

void
MDNLStatic::getSubStiffAndForce(int isub, DistrGeomState &geomState, 
                                DistrVector &res, DistrVector &elemIntForce, double lambda)
{
 SubDomain *sd = decDomain->getSubDomain(isub);

 StackVector residual(res.subData(isub), res.subLen(isub));

 // eIF = element internal force
 StackVector eIF(elemIntForce.subData(isub), elemIntForce.subLen(isub));

 sd->getStiffAndForce(*geomState[isub], eIF, allCorot[isub], kelArray[isub],
                      residual, lambda);
}

double
MDNLStatic::norm(DistrVector &vec)
{
 return sqrt(solver->getFNormSq(vec));
}

void
MDNLStatic::addMpcForces(DistrVector &vec)
{
  execParal1R(decDomain->getNumSub(), this, &MDNLStatic::subAddMpcForces, vec);
}

void
MDNLStatic::subAddMpcForces(int isub, DistrVector &vec)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  double *mpcForces = new double[sd->numMPCs()]; // don't delete  
  solver->getLocalMpcForces(isub, mpcForces);  // mpcForces set to incremental mpc lagrange multipliers
  StackVector localvec(vec.subData(isub), vec.subLen(isub));
  sd->constraintProductTmp(mpcForces, localvec); // C^T*lambda added to force residual
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
                             DistrVector&, double lambda)
{
 times->buildStiffAndForce -= getTime();

 execParal4R(decDomain->getNumSub(), this, &MDNLStatic::getSubStiffAndForce, geomState,
             residual, elementInternalForce, lambda);

 times->buildStiffAndForce += getTime();

 return sqrt(solver->getFNormSq(residual)); // XXXX this should include the active constraint functions' values
}

DistrGeomState*
MDNLStatic::createGeomState()
{
 times->timeGeom -= getTime();
 DistrGeomState* geomState = new DistrGeomState(decDomain);
 times->timeGeom += getTime();
 return geomState;
}

void
MDNLStatic::updatePrescribedDisplacement(DistrGeomState *geomState, double)
{
 times->timePresc -= getTime();
 execParal1R(decDomain->getNumSub(),this,&MDNLStatic::updatePrescribedDisp,*geomState);
 times->timePresc += getTime();
}

void
MDNLStatic::updatePrescribedDisp(int isub, DistrGeomState& geomState)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  sd->updatePrescribedDisp(geomState[isub], deltaLambda);
} 

int
MDNLStatic::reBuild(int iteration, int step, DistrGeomState& geomState)
{
 if(step == 1 && iteration == 0) return 0; // XXXX
 times->rebuild -= getTime();
 int rebuildFlag = 0;

 if (iteration % domain->solInfo().getNLInfo().updateK == 0) {

   GenMDDynamMat<double> allOps;
   allOps.sysSolver = solver;
   decDomain->rebuildOps(allOps, 0.0, 0.0, 1.0, kelArray);

/*
   GenMDDynamMat<double> allOps;
   decDomain->buildOps(allOps, 0.0, 0.0, 1.0, (Rbm **)0, kelArray);
   solver = (GenFetiSolver<double> *) allOps.sysSolver;
*/
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
 
 // Make subdomain's corotators
 times->corotatorTime -= getTime();
 allCorot = new Corotator**[numSub]; 
 execParal(numSub, this, &MDNLStatic::makeSubCorotators);
 times->corotatorTime += getTime();

 times->memoryPreProcess += threadManager->memoryUsed();

 // NOTE: count FETI memory separately from pre-process memory in
 // timing file

 // Construct FETI Solver, build GtG, subdomain matrices, etc.
 times->getFetiSolverTime -= getTime();
 GenMDDynamMat<double> allOps;
 decDomain->buildOps(allOps, 0.0, 0.0, 1.0);
 solver = (GenFetiSolver<double> *) allOps.sysSolver;
 times->getFetiSolverTime += getTime();

 // Make subdomain's array of stiffness matrices
 times->memoryPreProcess -= threadManager->memoryUsed();
 times->kelArrayTime -= getTime();
 kelArray = new FullSquareMatrix*[numSub];
 execParal(numSub, this, &MDNLStatic::makeSubKelArrays);
 times->kelArrayTime += getTime();
 times->memoryPreProcess += threadManager->memoryUsed();

 tolerance = domain->solInfo().getNLInfo().tolRes;

 domain->InitializeStaticContactSearch(decDomain->getNumSub(), decDomain->getAllSubDomains()); // YYYY
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

bool
MDNLStatic::linesearch()
{
  return domain->solInfo().getNLInfo().linesearch;
}

void
MDNLStatic::getRHS(DistrVector& rhs)
{
  // ... BUILD THE RHS FORCE (not including follower forces and internal force)
  times->formRhs -= getTime();
  execParal1R(decDomain->getNumSub(), this, &MDNLStatic::subGetRHS, rhs);
  times->formRhs += getTime(); 
}

void
MDNLStatic::subGetRHS(int isub, DistrVector& rhs)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector subrhs(rhs.subData(isub), rhs.subLen(isub));
  sd->computeConstantForce<double>(subrhs);
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
}

void
MDNLStatic::updateMpcRhs(DistrGeomState &geomState)
{
  decDomain->setContactGap(&geomState, solver);
}

void
MDNLStatic::updateContactConditions(DistrGeomState* geomState)
{

  domain->UpdateSurfaces(geomState, 1, decDomain->getAllSubDomains());
  domain->PerformStaticContactSearch();
  domain->deleteLMPCs();
  domain->ExpComputeMortarLMPC();
  //domain->printLMPC();
  domain->CreateMortarToMPC();
  ((GenFetiDPSolver<double> *) solver)->deleteG();
  decDomain->reProcessMPCs();
  ((GenFetiDPSolver<double> *) solver)->reconstructMPCs(decDomain->mpcToSub_dual, decDomain->mpcToMpc, decDomain->mpcToCpu);
  //solver = decDomain->getFetiSolver();

}

