#include <stdio.h>
#include <stdlib.h>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/Domain.h>
#include <Driver.d/Dynam.h>
#include <Problems.d/EigenDescr.h>

#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/mathUtility.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/Vector.h>
#include <Math.d/VectorSet.h>

#include <Solvers.d/Solver.h>
#include <Solvers.d/Rbm.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Corotational.d/GeomState.h>

int
SingleDomainEigen::solVecSize()
{
 // ... returns number of unconstrained dof
 return domain->numUncon();
}

int
SingleDomainEigen::solVecInfo()
{
 // ... returns number of unconstrained dof
 return domain->numUncon();
}

void
SingleDomainEigen::preProcess()
{

 // Allocate space for the Static Timers
 times = new StaticTimers;

 times->preProcess -= getTime();
 // Makes renumbering, connectivities and dofsets
 domain->preProcessing();

 // Total number of dof
 int numdof  = domain->numdof();

 times->makeBCs -= getTime();
 //int *bc  = (int *)    dbg_alloca(sizeof(int)*numdof);
 //bcx = (double *) dbg_alloca(sizeof(double)*numdof);
 int *bc  = new int[numdof]; //HB
 bcx = new double[numdof]; //HB

 // ... make boundary conditions
 domain->make_bc(bc,bcx);
 times->makeBCs += getTime();

 // ... make constrained dof set array
 times->makeDOFs -= getTime();
 domain->make_constrainedDSA();

 // ... construct all dofs
 domain->makeAllDOFs();
 times->makeDOFs += getTime();

 //ADDED FOR HEV PROBLEM, EC, 20070820
 if (domain->solInfo().HEV)  {
   domain->makeAllDOFsFluid();
 }
 
 // if we have initial displacements, we have to consider
 // the nonlinear tangent stiffness matrix instead of the
 // linear stiffness matrix. This is due to the prestress.

 geomKelArray = 0;
 if( domain->numInitDisp6() > 0 && domain->solInfo().gepsFlg == 1) {

   if(domain->solInfo().buckling) {
     fprintf(stderr," ... Buckling Analysis              ...\n");
     times->kelArrayTime -= getTime();
     domain->createKelArray(geomKelArray);
     times->kelArrayTime += getTime();
   }

   domain->computeGeometricPreStress(allCorot, geomState, kelArray, times,
                                     geomKelArray);

   if(domain->solInfo().buckling == 1) kelArray = 0; //HB: is there any memory leak here ? It seems to me that
                                                     //kelArray was allocated in domain->computeGeometricPreStress
 }

 times->preProcess += getTime();
}

SDEigenPostProcessor *
SingleDomainEigen::getPostProcessor()
{
 return new SDEigenPostProcessor(domain, times, bcx);
}

void
SingleDomainEigen::buildEigOps( DynamMat &dMat )
{
 AllOps<double> allOps;

 //if(addedMass) allOps.M = domain->constructDBSparseMatrix<double>((DofSetArray *)0, nodeToNode_added);
 //else 
 allOps.M = domain->constructDBSparseMatrix<double>();
 // Used for printing out K during debugging.
 allOps.K = domain->constructDBSparseMatrix<double>();

 // construct geometric rigid body modes if necessary
 //Rbm *rigidBodyModes = 0;
 if(domain->solInfo().rbmflg) { 
   //filePrint(stderr, " ... Constructing Geometric RBMs    ... \n");
   dMat.rigidBodyModes = domain->constructRbm();
 }

 if(domain->solInfo().hzemFlag) {
   //filePrint(stderr, " ... Constructing HZEMs             ... \n");
   dMat.rigidBodyModes = domain->constructHzem();
 }

 // construct rigid body modes for sloshing problems if necessary
 if(domain->solInfo().slzemFlag) { 
   //filePrint(stderr, " ... Constructing Sloshing RBMs in EigenDescr.C  ... \n");
   dMat.rigidBodyModes = domain->constructSlzem();
 }

 // build stiffness and mass matrices
 domain->buildOps<double>(allOps, 1.0, 0.0, 0.0, dMat.rigidBodyModes, kelArray);
 dMat.dynMat  = allOps.sysSolver;
 dMat.M       = allOps.M;

 if(domain->solInfo().explicitK)
   dMat.refK    = allOps.K;

 // If we are doing a pre-stress computation.
 if(geomKelArray && domain->solInfo().gepsFlg == 1) {
   fprintf(stderr," ... Assembling Geometric Stiffness ...\n");
   Connectivity *allDOFs = domain->getAllDOFs();
   dMat.M->zeroAll();
   int iele;
   for(iele=0; iele<domain->numElements(); ++iele) {
     dMat.M->add(geomKelArray[iele],(*allDOFs)[iele]);
   }
 }
}

void
SingleDomainEigen::reBuild( DynamMat &dMat )
{
 AllOps<double> allOps;
 dMat.M->zeroAll();
 allOps.M    = dMat.M ;
 allOps.sysSolver = dMat.dynMat;

 // rebuild stiffness and mass matrices
 // watch: no rigid body modes assumed

 domain->rebuildOps<double>(allOps, 1.0, 0.0, 0.0);
}

int SingleDomainEigen::getNumEigen()
{
  return domain->solInfo().nEig;
}

int SingleDomainEigen::getEigenSolverType()
{
  return domain->solInfo().eigenSolverType;
}

int SingleDomainEigen::getEigenSolverSubType()
{
  return domain->solInfo().eigenSolverSubType;
}

void
SingleDomainEigen::getSubSpaceInfo(int& subspacesize, int& maxIter, 
                                   double& tolEig, double& tolJac, bool &explicitK )
{
 subspacesize = domain->solInfo().subspaceSize;
 tolEig       = domain->solInfo().tolEig;
 tolJac       = domain->solInfo().tolJac;

 int numIter1 = domain->solInfo().maxitEig;
 int numIter2 = 10*subspacesize;

 maxIter = myMax(numIter1, numIter2);

 if (numIter1 > 0) maxIter = numIter1;

 explicitK = domain->solInfo().explicitK;
}

void
SingleDomainEigen::error(int subspacesize, int numRbm)
{
 fprintf(stderr,"*** subspace size must be bigger than number of rbm   ***\n");
 fprintf(stderr,"subspace size = %d\n",subspacesize);
 fprintf(stderr,"number of rbm = %d\n",numRbm);
}

#ifdef WINDOWS
#define drand48 rand
#endif

void
SingleDomainEigen::initQ(Vector* Q, int size)
{
  int numdof = domain->numUncon();
  int i, j;
  if(size == 0) {
    fprintf(stderr," *** ERROR: Division by zero within initQ routine. \n");
    return;
  }

  Vector t(numdof);

  // This has been changed to uses random vectors for the initial subspace.
  // The old way (commented out) could lead to a linear dependence between
  // RBMs and initial vectors.
  // WARNING: does drand48() exist on all machines?
  for(i=0; i<size; ++i) {
    t.zero();
    for(j=0; j<numdof; ++j)
      t[j] = drand48();
    Q[i] = t;
  }
}

void
SingleDomainEigen::printTimers(Solver *solver)
{
 long memoryUsed = solver->size();
 double solveTime  = solver->getSolutionTime();

 times->printStaticTimers(solveTime, memoryUsed, domain);
}

void
SDEigenPostProcessor::eigenOutput(Vector& eigenValues, VectorSet& eigenVectors, int convEig)
{
 times->output -= getTime(); 
 domain->eigenOutput(eigenValues,eigenVectors,bcx,convEig);
 times->output += getTime(); 
}

