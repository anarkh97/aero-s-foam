#include <Utils.d/dbg_alloca.h>
#include <cstdio>
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

NonLinStatic::NonLinStatic(Domain *d)
{
  domain = d;
  kelArray = 0;
  allCorot = 0;
  bcx = 0;
  solver = solver;
  prec = 0;

  if(domain->GetnContactSurfacePairs())
     domain->InitializeStaticContactSearch(MortarHandler::CTC);
}

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
NonLinStatic::getStiffAndForce(GeomState& geomState, Vector& residual, Vector& elementInternalForce, 
                               Vector &, double lambda, GeomState *refState)
{
  times->buildStiffAndForce -= getTime();

  if(domain->GetnContactSurfacePairs()) {
    domain->UpdateSurfaces(MortarHandler::CTC, &geomState);
    domain->PerformStaticContactSearch(MortarHandler::CTC);
    domain->deleteSomeLMPCs(mpc::ContactSurfaces);
    domain->ExpComputeMortarLMPC(MortarHandler::CTC);
    domain->UpdateContactSurfaceElements();

    if(solver) delete solver;
    if(prec) delete prec;
    if(allCorot) delete [] allCorot; allCorot = 0;  // memory leak?
    if(kelArray) delete [] kelArray; kelArray = 0;
    preProcess(false); // TODO consider case domain->solInfo().getNLInfo().updateK > 1
    elementInternalForce.initialize(domain->maxNumDOF());
  }

  domain->getStiffAndForce(geomState, elementInternalForce, allCorot, 
                           kelArray, residual, lambda, 0, refState);

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
    filePrint(stderr,"----------------------------------------------------\n");
    filePrint(stderr,"Newton Iter    #%d\tcurrent dv   = % e\n \t\t\t"
                   "first dv     = % e\n \t\t\trelative dv  = % e\n",
                    iter+1, normDv, firstDv, relativeDv);
    filePrint(stderr,"                \tcurrent Res  = % e\n \t\t\t"
                   "first Res    = % e\n \t\t\trelative Res = % e\n",
                    normRes, firstRes, relativeRes);

    filePrint(stderr,"----------------------------------------------------\n");
  }

 int converged = 0;

 // Check relative convergence criteria
 if(normRes <= tolerance*firstRes)
   converged = 1;

 // Check Divergence
 if(normRes > 10000*firstRes) converged = -1;

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
   geomState = new GeomState( *domain->getDSA(),*domain->getCDSA(),domain->getNodes(),&domain->getElementSet()); 

 times->timeGeom += getTime();

 return geomState;
}

int
NonLinStatic::reBuild(int iteration, int step, GeomState&)
{
 times->rebuild -= getTime();

 int rebuildFlag = 0;

 if (iteration % domain->solInfo().getNLInfo().updateK == 0) {
   //PJSA 11/5/09: new way to rebuild solver (including preconditioner), now works for any solver
   spm->zeroAll();
   AllOps<double> ops;
   if(spp) { // rebuild preconditioner as well as the solver
     spp->zeroAll();
     ops.spp = spp;
   }
   domain->makeSparseOps<double>(ops, 1.0, 0.0, 0.0, spm, kelArray);
   solver->factor();
   if(prec) prec->factor();
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

bool
NonLinStatic::linesearch()
{
 return domain->solInfo().getNLInfo().linesearch;
}

void
NonLinStatic::getRHS(Vector& rhs)
{
 // ... BUILD THE RHS FORCE (not including follower or internal forces)
 times->formRhs -= getTime();
 domain->computeConstantForce<double>(rhs);
 times->formRhs += getTime();
}

void
NonLinStatic::preProcess(bool factor)
{
 // Allocate space for the Static Timers
 times = new StaticTimers;

 startTimerMemory(times->preProcess, times->memoryPreProcess);

 times->timePre -= getTime();
 domain->preProcessing();
 times->timePre += getTime();

 int numdof = domain->numdof();

 times->makeBCs -= getTime();
 int *bc = (int *) dbg_alloca(sizeof(int)*numdof);
 bcx = new double[numdof];

 // Make the boundary conditions info
 domain->make_bc(bc, bcx);

 times->makeBCs += getTime();

 // Now, call make_constrainedDSA(bc) to 
 // built c_dsa that will incorporate all 
 // the boundary conditions info

 times->makeDOFs -= getTime();
 domain->make_constrainedDSA(bc);
 domain->makeAllDOFs();
 times->makeDOFs += getTime();

 stopTimerMemory(times->preProcess, times->memoryPreProcess);
 AllOps<double> allOps;

 long buildMem = -memoryUsed();
 times->timeBuild -= getTime();

 Rbm *rigidBodyModes = 0;
 if(domain->solInfo().rbmflg) rigidBodyModes = domain->constructRbm(); // new policy is to construct rbms if GRBM is requested in input file
                                                                       // but only use them when it is appropriate to do so. In nonlinear statics it is not
                                                                       // since the nullity of the tangent stiffness matrix may be less than the nullity
                                                                       // of the number of rigid body modes
 
 domain->buildOps<double>(allOps, 1.0, 0.0, 0.0, (Rbm *) 0, kelArray, factor);
 times->timeBuild += getTime();
 buildMem += memoryUsed();

 solver = allOps.sysSolver;
 spm = allOps.spm;
 prec = allOps.prec;
 spp = allOps.spp;

 // ... ALLOCATE MEMORY FOR THE ARRAY OF COROTATORS
 startTimerMemory(times->preProcess, times->memoryPreProcess);
 times->corotatorTime -= getTime();
 allCorot = new Corotator *[domain->numElements()];

 // ... CREATE THE ARRAY OF POINTERS TO COROTATORS
 domain->createCorotators(allCorot);
 times->corotatorTime += getTime();

 // ... CREATE THE ARRAY OF ELEMENT STIFFNESS MATRICES
 times->kelArrayTime -= getTime();
 domain->createKelArray(kelArray);
 times->kelArrayTime += getTime();
 stopTimerMemory(times->preProcess, times->memoryPreProcess);

 // Set the nonlinear tolerance used for convergence
 tolerance = domain->solInfo().getNLInfo().tolRes;
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
                           Vector &, GeomState *refState)
{
  times->output -= getTime();
  Vector dummyForce(domain->numUncon(), 0.0);
  domain->postProcessing(geomState, force, dummyForce, lambda, 1, 0, 0, allCorot,
                         (FullSquareMatrix *) 0, (double *) 0, (double *) 0, refState);
  times->output += getTime();
}

double
NonLinStatic::getEnergy(double lambda, Vector& force, GeomState* geomState)
{
  Vector sol(domain->numUncon(), 0.0);
  for(int i=0; i<domain->numNode(); ++i) {
    Node *node_i = domain->getNodes()[i];
    int xloc  = domain->getCDSA()->locate(i, DofSet::Xdisp);
    if(xloc >= 0 && node_i) {
      sol[xloc]  = ( (*geomState)[i].x - node_i->x);
    }
    int yloc  = domain->getCDSA()->locate(i, DofSet::Ydisp);
    if(yloc >= 0 && node_i)
      sol[yloc]  = ( (*geomState)[i].y - node_i->y);
    int zloc  = domain->getCDSA()->locate(i, DofSet::Zdisp);
    if(zloc >= 0 && node_i)
      sol[zloc]  = ( (*geomState)[i].z - node_i->z);
    double rot[3];
    mat_to_vec((*geomState)[i].R,rot);
    int xrot  = domain->getCDSA()->locate(i, DofSet::Xrot);
    if(xrot >= 0 && node_i)
      sol[xrot]  = rot[0];
    int yrot  = domain->getCDSA()->locate(i, DofSet::Yrot);
    if(yrot >= 0 && node_i)
      sol[yrot]  = rot[1];
    int zrot  = domain->getCDSA()->locate(i, DofSet::Zrot);
    if(zrot >= 0 && node_i)
      sol[zrot]  = rot[2];
  }

  // compute external energy not including follower forces
  double Wext = -lambda*force*sol;

  // compute internal energy and external energy due to follower forces
  // XXXX note: need to include pressure, gravity and thermal forces in getElementEnergy for this to work
  double Wela = 0.0;
  for(int i = 0; i < domain->numElements(); ++i)
     Wela += allCorot[i]->getElementEnergy(*geomState, domain->getNodes());
  //cerr << "Wext = " << Wext << ", Wela = " << Wela << endl;

  // Total Energy = Wext + Wela
  return Wext + Wela;
}

double
NonLinStatic::getResidualNorm(Vector &res)
{
  CoordinateMap *m = dynamic_cast<CoordinateMap *>(solver);
  if(m) return m->norm(res);
  else return res.norm();
}
