#include <cstdio>

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
                                DistrVector &res, DistrVector &elemIntForce, double lambda,
                                DistrGeomState *refState)
{
 SubDomain *sd = decDomain->getSubDomain(isub);

 StackVector residual(res.subData(isub), res.subLen(isub));

 // eIF = element internal force
 StackVector eIF(elemIntForce.subData(isub), elemIntForce.subLen(isub));

 GeomState *subRefState = (refState) ? (*refState)[isub] : 0;
 sd->getStiffAndForce(*geomState[isub], eIF, allCorot[isub], kelArray[isub],
                      residual, lambda, 0, subRefState);
}

double
MDNLStatic::getResidualNorm(DistrVector &r, DistrGeomState &geomState)
{
 // XXXX should also include constraint error i.e. pos_part<gap>
 DistrVector w(r);
 execParal1R(decDomain->getNumSub(), this, &MDNLStatic::addConstraintForces, w); // w = r + C^T*lambda
                  // note C = grad(gap) has already been updated in getStiffAndForce.
                  // XXXX need to make sure lambda_i is correctly mapped to C_i. I think this is done
                  // correctly only for the case of one contactsurfaces pair
 return sqrt(solver->getFNormSq(w));
}

void
MDNLStatic::addConstraintForces(int isub, DistrVector &vec)
{
  // I need to treat the contact forces from CONTACTSURFACES separately due to search,
  // the ith lagrange multiplier at iteration n may not correspond to the ith constraint
  // after updating the contact surfaces
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector localvec(vec.subData(isub), vec.subLen(isub));
  sd->addConstraintForces(mu[isub], lambda[isub], localvec);      // C^T*lambda added to vec
}

void
MDNLStatic::makeSubKelArrays(int isub)
{
 SubDomain *sd = decDomain->getSubDomain(isub);
 sd->createKelArray(kelArray[isub]);
}

void
MDNLStatic::deleteSubKelArrays(int isub)
{
 if(kelArray[isub]) delete [] kelArray[isub];
}

void
MDNLStatic::makeSubCorotators(int isub)
{
 SubDomain *sd  = decDomain->getSubDomain(isub);
 int numele     = sd->numElements();
 allCorot[isub] = new Corotator*[numele];
 sd->createCorotators(allCorot[isub]);
}

void
MDNLStatic::deleteSubCorotators(int isub)
{
 SubDomain *sd = decDomain->getSubDomain(isub);
 if(allCorot[isub]) {
   for (int iElem = 0; iElem < sd->numElements(); ++iElem) {
     if(allCorot[isub][iElem] && (allCorot[isub][iElem] != dynamic_cast<Corotator*>(sd->getElementSet()[iElem])))
       delete allCorot[isub][iElem];
   }
   delete [] allCorot[isub];
 }
}


MDNLStatic::MDNLStatic(Domain *d)
{
 domain = d;

#ifdef DISTRIBUTED
 decDomain = new GenDistrDomain<double>(domain);
#else
 decDomain = new GenDecDomain<double>(domain);
#endif
 numSystems = 0;
 mu = 0; lambda = 0;
 solver = 0;
 kelArray = 0;
 allCorot = 0;
 times = 0;
}

MDNLStatic::~MDNLStatic()
{
  if(mu) delete [] mu;
  if(lambda) delete [] lambda;
  if(solver) delete solver;
  if(times) delete times;
  if(allCorot) {
    execParal(decDomain->getNumSub(), this, &MDNLStatic::deleteSubCorotators);
    delete [] allCorot;
  }
  if(kelArray) {
    execParal(decDomain->getNumSub(), this, &MDNLStatic::deleteSubKelArrays);
    delete [] kelArray;
  }
  if(decDomain) delete decDomain;
}

DistrInfo&
MDNLStatic::solVecInfo()
{
 return decDomain->solVecInfo();
}

DistrInfo&
MDNLStatic::sysVecInfo()
{
 return decDomain->sysVecInfo();
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

 if(verboseFlag) {
   filePrint(stderr,"----------------------------------------------------\n");
   filePrint(stderr,"Nonlinear Iter #%d\tcurrent dv   = % e\n \t\t\t"
                    "first dv     = % e\n \t\t\trelative dv  = % e\n",
                     iter+1, normDv, firstDv, relativeDv);
   filePrint(stderr,"                \tcurrent Res  = % e\n \t\t\t"
                    "first Res    = % e\n \t\t\trelative Res = % e\n",
                     normRes, firstRes, relativeRes);

   filePrint(stderr,"----------------------------------------------------\n");
 }
 
 int converged = 0;

 // KHP: Charbel requested convergence be monitored based on residual only.
 // Check incremental displacement
 //if(relativeDv <= domain->solInfo().getNLInfo().tolRes)
 //  converged = 1;

 // Check to see if residual has converged
 if(relativeRes <= domain->solInfo().getNLInfo().tolRes && relativeDv <= domain->solInfo().getNLInfo().tolInc)
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
                             DistrVector&, double _lambda, DistrGeomState *refState)
{
 times->buildStiffAndForce -= getTime();

 updateConstraintTerms(&geomState, _lambda);

 execParal5R(decDomain->getNumSub(), this, &MDNLStatic::getSubStiffAndForce, geomState,
             residual, elementInternalForce, _lambda, refState);

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
 solver = allOps.sysSolver;
 if(allOps.K) delete allOps.K;
 if(allOps.Kuc) delete allOps.Kuc;
 if(allOps.M) delete allOps.M;
 if(allOps.Muc) delete allOps.Muc;
 if(allOps.C) delete allOps.C;
 if(allOps.Cuc) delete allOps.Cuc;
 times->getFetiSolverTime += getTime();

 // Make subdomain's array of stiffness matrices
 times->memoryPreProcess -= threadManager->memoryUsed();
 times->kelArrayTime -= getTime();
 kelArray = new FullSquareMatrix*[numSub];
 execParal(numSub, this, &MDNLStatic::makeSubKelArrays);
 times->kelArrayTime += getTime();
 times->memoryPreProcess += threadManager->memoryUsed();

 tolerance = domain->solInfo().getNLInfo().tolRes;

 domain->InitializeStaticContactSearch(MortarHandler::CTC, decDomain->getNumSub(), decDomain->getAllSubDomains());

 mu = new std::map<std::pair<int,int>,double>[decDomain->getNumSub()];
 lambda = new std::vector<double>[decDomain->getNumSub()];
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

ParallelSolver *
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
                         DistrVector &Force, DistrVector &, DistrGeomState *refState)
{
  startTimerMemory(times->output, times->memoryOutput);
  decDomain->postProcessing(geomState, allCorot, lambda, (SysState<GenDistrVector<double> > *) 0,
                            (GenDistrVector<double> *) 0, refState);
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

  delete mdpp;

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
MDNLStatic::updateConstraintTerms(DistrGeomState* geomState, double _lambda)
{
  GenFetiDPSolver<double> *fetiSolver = dynamic_cast<GenFetiDPSolver<double> *>(solver);
  if(fetiSolver) {
    execParal(decDomain->getNumSub(), this, &MDNLStatic::getConstraintMultipliers);
    if(domain->GetnContactSurfacePairs()) {
      domain->UpdateSurfaces(MortarHandler::CTC, geomState, decDomain->getAllSubDomains());
      domain->PerformStaticContactSearch(MortarHandler::CTC);
      domain->deleteSomeLMPCs(mpc::ContactSurfaces);
      domain->ExpComputeMortarLMPC(MortarHandler::CTC);
      domain->CreateMortarToMPC();
      decDomain->reProcessMPCs();
      fetiSolver->reconstructMPCs(decDomain->mpcToSub_dual, decDomain->mpcToMpc, decDomain->mpcToCpu);
    }
    // set the gap for the linear constraints
    decDomain->setConstraintGap(geomState, fetiSolver, _lambda);
  }
}

void
MDNLStatic::getConstraintMultipliers(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  mu[isub].clear();
  lambda[isub].clear();
  sd->getConstraintMultipliers(mu[isub], lambda[isub]);
}

void
MDNLStatic::updateStates(DistrGeomState *refState, DistrGeomState& geomState)
{
  execParal2R(decDomain->getNumSub(), this, &MDNLStatic::subUpdateStates, refState, &geomState);
}

void
MDNLStatic::subUpdateStates(int isub, DistrGeomState *refState, DistrGeomState *geomState)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  GeomState *subRefState = (refState) ? (*refState)[isub] : 0;
  sd->updateStates(subRefState, *(*geomState)[isub], allCorot[isub]);
}

