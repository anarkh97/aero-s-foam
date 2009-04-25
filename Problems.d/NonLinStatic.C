#include <Utils.d/dbg_alloca.h>
#include <stdio.h>
#include <Utils.d/Memory.h>

#include <Driver.d/Domain.h>
#include <Problems.d/NonLinStatic.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/TemperatureState.h>
#include <Solvers.d/Solver.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/FullSquareMatrix.h>
#include <Timers.d/GetTime.h>

extern int verboseFlag;

int
NonLinStatic::solVecInfo()
{
 return domain->numUncon();
}

int
NonLinStatic::sysVecInfo()
{
 return 0;
}

double
NonLinStatic::getStiffAndForce(GeomState& geomState, 
                         Vector& residual,Vector& elementInternalForce, 
                         Vector &)
{
 times->buildStiffAndForce -= getTime();


 domain->getStiffAndForce(geomState, elementInternalForce, allCorot, 
                          kelArray, residual);

 times->buildStiffAndForce += getTime();

 return sqrt(residual*residual);
}

void
NonLinStatic::updatePrescribedDisplacement(GeomState *geomState, double)
{
 // Measure the time necessary to update the prescribed displacments

 times->timePresc -= getTime();

 int numDirichlet = domain->nDirichlet();

 BCond *dbc = domain->getDBC();

 double delta = domain->solInfo().getNLInfo().dlambda;

 geomState->updatePrescribedDisplacement(dbc, numDirichlet, delta);

 times->timePresc += getTime();

}

int
NonLinStatic::checkConvergence(int iter, double normDv, double normRes)
{
 // Measure time necessary to check for convergence of Newton algorithm

 times->timeCheck -= getTime();

 if(iter == 0) {
   firstDv  = normDv;
   firstRes = normRes;
 }

 double relativeDv  = normDv/firstDv;
 double relativeRes = normRes/firstRes;

 if(verboseFlag)
  {
    fprintf(stderr,"----------------------------------------------------\n");
    fprintf(stderr,"Newton Iter    #%d\tcurrent dv   = % e\n \t\t\t"
                   "first dv     = % e\n \t\t\trelative dv  = % e\n",
                    iter+1, normDv, firstDv, relativeDv);
    fprintf(stderr,"                \tcurrent Res  = % e\n \t\t\t"
                   "first Res    = % e\n \t\t\trelative Res = % e\n",
                    normRes, firstRes, relativeRes);

    fprintf(stderr,"----------------------------------------------------\n");
  }

 int converged = 0;

 // Check relative convergence criteria
 if(normRes <= tolerance*firstRes)
   converged = 1;

 // Check Divergence
 if(normRes > 1000*firstRes) converged = -1;

 // Store residual norm and dv norm for output
 times->norms[iter].normDv      = normDv;
 times->norms[iter].relativeDv  = relativeDv;
 times->norms[iter].normRes     = normRes;
 times->norms[iter].relativeRes = relativeRes;

 times->timeCheck += getTime();

 return converged;

}

GeomState*
NonLinStatic::createGeomState()
{
 times->timeGeom -= getTime();

 GeomState *geomState;
 if(domain->solInfo().soltyp == 2) 
   geomState = (GeomState *) new TemperatureState( *domain->getDSA(),*domain->getCDSA(),domain->getNodes());
 else
   geomState = new GeomState( *domain->getDSA(),*domain->getCDSA(),domain->getNodes()); 

 times->timeGeom += getTime();

 //fprintf(stderr," ... Time to build GeomState %e\n",times->timeGeom/1000.0);
 //fflush(stderr);

 return geomState;
}

int
NonLinStatic::reBuild(int iteration, int step, GeomState&)
{

 // NOTE: rebuild also includes factoring

 times->rebuild -= getTime();
  
 int rebuildFlag = 0;

 // KHP: MODIFICATION
 if( iteration % domain->solInfo().getNLInfo().updateK == 0 )  {
   cerr << "REBUILDING SOLVER\n";
   solver->reBuild(kelArray);
   rebuildFlag = 1;
 }
 times->rebuild += getTime();

 return rebuildFlag;
}

int
NonLinStatic::elemVecInfo()
{
  return domain->maxNumDOF();
}

int
NonLinStatic::getMaxit()
{
  return domain->solInfo().getNLInfo().maxiter;
}

double
NonLinStatic::getScaleFactor()
{
 return domain->solInfo().getNLInfo().lfactor;
}

double
NonLinStatic::getDeltaLambda0()
{
 return domain->solInfo().getNLInfo().dlambda;
}

double
NonLinStatic::getMaxLambda()
{
 return domain->solInfo().getNLInfo().maxLambda;
}

void
NonLinStatic::getRHS(Vector& rhs, GeomState *gs)
{
 // ... BUILD THE RHS FORCE (external + gravity + nonhomogeneous)
 times->formRhs -= getTime();
 domain->buildRHSForce<double>(rhs, 0, gs);
 times->formRhs += getTime();
}

void
NonLinStatic::preProcess()
{
 // Allocate space for the Static Timers
 times = new StaticTimers;

 //times->preProcess -= getTime();
 startTimerMemory(times->preProcess, times->memoryPreProcess);

 // Makes renumbering, connectivities and dofsets
 fprintf(stderr," ... Constructing Connectivities    ...\n");
 fflush(stderr);

 times->timePre -= getTime();
 domain->preProcessing();
 times->timePre += getTime();

 //fprintf(stderr," ... Time to construct Renumbering & Connectivities %e\n", times->timePre/1000.0);
 //fflush(stderr);

 int numdof = domain->numdof();

 fprintf(stderr," ... Making Boundary Conditions     ...\n");
 fflush(stderr);

 times->makeBCs -= getTime();

 int *bc = (int *) dbg_alloca(sizeof(int)*numdof);
 bcx = new double[numdof];

 // Make the boundary conditions info
 domain->make_bc(bc, bcx);

 times->makeBCs += getTime();

 //fprintf(stderr," ... Time to make Boundary Conditions %e\n", times->makeBCs/1000.0);
 //fflush(stderr);

 // Now, call make_constrainedDSA(bc) to 
 // built c_dsa that will incorporate all 
 // the boundary conditions info

 fprintf(stderr," ... Making Degrees of Freedom      ...\n");
 fflush(stderr);

 times->makeDOFs -= getTime();
 domain->make_constrainedDSA(bc);
 domain->makeAllDOFs();
 times->makeDOFs += getTime();

 //fprintf(stderr," ... Time to make Degrees of Freedom %e\n", times->makeDOFs/1000.0);
 //fflush(stderr);

 stopTimerMemory(times->preProcess, times->memoryPreProcess);
 AllOps<double> allOps;

 double Kcoef = 1.0;
 double Mcoef = 0.0;
 double Ccoef = 0.0;

 fprintf(stderr," ... Building Solver                ...\n");
 fflush(stderr);

 long buildMem = -memoryUsed();
 times->timeBuild -= getTime();

 Rbm *rigidBodyModes = 0;
 if(domain->solInfo().rbmflg) 
   rigidBodyModes = domain->constructRbm(domain->probType());
 
 domain->buildOps<double>(allOps, Kcoef, Mcoef, Ccoef, rigidBodyModes);
 times->timeBuild += getTime();
 buildMem += memoryUsed();

 //fprintf(stderr," ... Time to build solver %e\n",times->timeBuild/1000.0);
 //fprintf(stderr," ... Memory to build solver %12.4f Mb\n", (double)buildMem/(1024*1024));
 //fflush(stderr);

 solver = allOps.sysSolver;

 fprintf(stderr," ... Creating Element Corotators    ...\n");
 fflush(stderr);

 // ... ALLOCATE MEMORY FOR THE ARRAY OF COROTATORS
 startTimerMemory(times->preProcess, times->memoryPreProcess);
 times->corotatorTime -= getTime();
 allCorot = new Corotator *[domain->numElements()];

 // ... CREATE THE ARRAY OF POINTERS TO COROTATORS
 domain->createCorotators(allCorot);
 times->corotatorTime += getTime();

 //fprintf(stderr," ... Time to create Element Corotators %e\n", times->corotatorTime/1000.0);
 fflush(stderr);

 // ... CREATE THE ARRAY OF ELEMENT STIFFNESS MATRICES
 fprintf(stderr," ... Creating Elem. Stiff. Array    ...\n");
 fflush(stderr);

 times->kelArrayTime -= getTime();
 domain->createKelArray(kelArray);
 times->kelArrayTime += getTime();
 stopTimerMemory(times->preProcess, times->memoryPreProcess);

 //fprintf(stderr," ... Time to create Element Stiffness Array %e\n", times->kelArrayTime/1000.0);
 //fflush(stderr);

 // Set the nonlinear tolerance used for convergence
 tolerance = domain->solInfo().getNLInfo().tolRes;

 //times->preProcess += getTime();
}

Solver *
NonLinStatic::getSolver()
{
  return solver;
}

SingleDomainPostProcessor<double, Vector, Solver> *
NonLinStatic::getPostProcessor()
{
 return new SingleDomainPostProcessor<double,Vector,Solver>(domain,bcx,times,solver);
}

void
NonLinStatic::printTimers()
{
 times->timeTimers -= getTime();

 long memoryUsed = solver->size();
 double solveTime  = solver->getSolutionTime();

 times->printStaticTimers(solveTime, memoryUsed, domain);

 times->timeTimers += getTime();


}

void
NonLinStatic::staticOutput(GeomState *geomState, double lambda, Vector& force,
                Vector &)
{

  times->output -= getTime();
  Vector dummyForce(domain->numUncon(), 0.0);
  domain->postProcessing(geomState, force, dummyForce, lambda, 1, 0, 0, allCorot);
  times->output += getTime();
}
