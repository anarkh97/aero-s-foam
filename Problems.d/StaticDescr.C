#include <Utils.d/dbg_alloca.h>
#include <iostream>
#include <cstdio>
#include <cmath>
// #include <Problems.d/StaticDescr.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Corotational.d/GeomState.h>
#include <Solvers.d/Rbm.h>

#include <Math.d/FullMatrix.h>
#include <Sfem.d/Sfem.h>

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#endif

extern Sfem *sfem;
typedef FSFullMatrix FullMatrix;
template <class Scalar> class GenVector;
typedef GenVector<DComplex> ComplexVector;
//#define DEBUG_RBM_FILTER
template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::projector_prep(Rbm *rbms)
{
 numR = rbms->numRBM(); 
 if (!numR) return;

 // KHP: store this pointer to the RBMs to use in the actual
 //      projection step 
 int useProjector=domain->solInfo().filterFlags;
 if(useProjector) {
   filePrint(stderr," ... Building the RBM Projector     ...\n");
 }

 int useHzemFilter = domain->solInfo().hzemFilterFlag;
 if(useHzemFilter) {
   filePrint(stderr," ... Building the HZEM Projector    ...\n");
 }
 int ndof = solVecInfo();
 Rmem = new double[numR*ndof];
 if (useProjector) rbms->getRBMs(Rmem);
 if (useHzemFilter) {
   for(int n=0; n<ndof; ++n) Rmem[n] = 1.;
 }

 //ADDED FOR SLOSHING PROBLEM, EC, 20070723
 int useSlzemFilter = domain->solInfo().slzemFilterFlag;
 if(useSlzemFilter) {
   filePrint(stderr," ... Building the SLZEM Projector    ...\n");
   filePrint(stderr," ... THIS PART OF THE CODE HAS NOT BEEN DEBUGGED!!!    ...\n");
 }

 if (useSlzemFilter) {
   for(int n=0; n<ndof; ++n) Rmem[n] = 1.;
 }

 StackFSFullMatrix Rt(numR, ndof, Rmem);

/*
 int i,j;
 for(i=0; i<numR; ++i) {
   double sum = 0.0;
   for(j=0; j<ndof; ++j)
     sum +=  Rt[i][j]*Rt[i][j];
   sum = sqrt(sum);
   for(j=0; j<ndof; ++j)
     Rt[i][j] *= (1.0/sum);
 }
*/
 FullMatrix R = Rt.transpose();
#ifdef DEBUG_RBM_FILTER
 R.print("R","R");
#endif

 FullMatrix RtR(numR,numR);
 Rt.mult(R,RtR);
#ifdef DEBUG_RBM_FILTER
 RtR.print("RtIR","RtIR");
#endif

 FullMatrix RtRinverse = RtR.invert();
#ifdef DEBUG_RBM_FILTER
 RtRinverse.print("RtRinverse");
#endif

 X = new FullMatrix(ndof,numR);
 R.mult(RtRinverse,(*X));
#ifdef DEBUG_RBM_FILTER
 X->print("RRtRinverse");
#endif
}

inline void applyMult(FSFullMatrix &M, GenVector<double> &x, GenVector<double> &y)
{
 M.mult(x,y);
}
inline void applyMult(FSFullMatrix &M, ComplexVector &x, ComplexVector &y)
{
 fprintf(stderr, "Programming error: Project was called on complex vector\n");
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::project(VectorType &f)
{
 if (!numR) return;

 int ndof = f.size();

 typedef typename VectorType::DataType Scalar;
 Scalar *yMem = (Scalar *) dbg_alloca(ndof*sizeof(Scalar));
 Scalar *zMem = (Scalar *) dbg_alloca(ndof*sizeof(Scalar));

 GenStackVector<Scalar> y(ndof,yMem);
 GenStackVector<Scalar> z(ndof,zMem);

 StackFSFullMatrix Rt(numR, ndof, Rmem);

 // y = Rt*f
 y.zero(); applyMult(Rt,f,y);

#ifdef DEBUG_RBM_FILTER
 y.print("Y");
#endif

 // z = X*y
 z.zero(); applyMult(*X,y,z);

 // f = f - z;
 f.linC(1.0, f, -1.0, z);

 // check R^T*f = 0
 //y.zero(); applyMult(Rt,f,y); cerr << "5. y = "; print_debug(y);
}


template<class T, class VectorType, class SolverType>
int
SingleDomainStatic<T, VectorType, SolverType>::solVecInfo()
{
 int ret = domain->numUncon();
 if(domain->solInfo().inpc) ret *= sfem->getP();
 return ret;
}


template<class T, class VectorType, class SolverType>
int
SingleDomainStatic<T, VectorType, SolverType>::solVecInfo(int i)
{
 int ret = (domain->numUncon());
 ret *= i;
 return ret;
}


template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::clean()
{
// RT
  solver->clean_up();
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::getRHS(VectorType &rhs)
{
 filePrint(stderr," ... Building the Force             ...\n");

 // ... BUILD THE RHS FORCE (external + gravity + nonhomogeneous)
 startTimerMemory(times->formRhs, times->memoryRhs);
 domain->template buildRHSForce<T>(rhs, kuc);
 
 // rigid body mode projector (or eigen mode projector)
 bool useProjector = (domain->solInfo().filterFlags || domain->solInfo().modeFilterFlag);
 if(useProjector) 
   project(rhs);

 stopTimerMemory(times->formRhs, times->memoryRhs);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::getRHSinpc(VectorType &rhs)
{
 filePrint(stderr," ... Building the Force   (inpc)          ...\n");

 // ... BUILD THE RHS FORCE (external + gravity + nonhomogeneous)
 startTimerMemory(times->formRhs, times->memoryRhs);
 rhs.zero();
 for (int i=0; i<(*allOps.rhs_inpc).size(); ++i) rhs[i] = (*allOps.rhs_inpc)[i]; 
 // print forces   
/* char forcefile[20];
 FILE *fileforce;
 fileforce=fopen("forcefile","w");
 for (int j=0;j<rhs.size();j++) fprintf(fileforce,"%g\n",rhs[j]);
 fclose(fileforce);
*/
 int useProjector=domain->solInfo().filterFlags;
 if(useProjector)
   project(rhs);

 stopTimerMemory(times->formRhs, times->memoryRhs);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::preProcess()
{
// filePrint(stderr," ... SingleDomainStatic<T, VectorType, SolverType>::preProcess() called         ...\n");
 // Allocate space for the Static Timers
 times = new StaticTimers;

 startTimerMemory(times->preProcess, times->memoryPreProcess);

 // Makes renumbering, connectivities and dofsets
 domain->preProcessing();

 int numdof = domain->numdof();

 times->makeBCs -= getTime();
 int *bc = (int *) dbg_alloca(sizeof(int)*numdof);
 bcx     = new T[numdof];

 // Make boundary conditions info
 //if(domain->solInfo().isAcousticHelm()) { // PJSA 1-15-08
 if(domain->getImplicitFlag() || domain->nCDirichlet()) {
   if(domain->getImplicitFlag()) bcxC = new DComplex[numdof * domain->getNumWaveDirections()];
   else bcxC = new DComplex[numdof];
   ((HData *)domain)->make_bc(domain, bc, bcxC);
   for(int i=0; i<numdof; ++i) ScalarTypes::copy(bcx[i],bcxC[i]); // temp fix, needs to be done for every direction before post processing
 }
 else domain->make_bc(bc,bcx);
 times->makeBCs += getTime();

 // Now, call make_constrainedDSA(bc) to  built c_dsa 
 // that will incorporate all the boundary conditions info
 times->makeDOFs -= getTime();
 domain->make_constrainedDSA(bc);
 domain->makeAllDOFs();
 times->makeDOFs += getTime();

 // if we have initial displacements, we have to consider
 // the nonlinear tangent stiffness matrix instead of the
 // linear stiffness matrix. This is due to the prestress.

 kelArray  = 0;
 geomState = 0;
 allCorot  = 0;

 if(domain->numInitDisp6() > 0 && domain->solInfo().gepsFlg == 1) {
   filePrint(stderr," ... Static Problem with Initial Displacement %d\n",
             domain->numInitDisp6());
   FullSquareMatrix *geomKelArray=0, *dummy=0;
   domain->computeGeometricPreStress(allCorot, geomState, kelArray, times,
                                     geomKelArray, dummy);
 }

 stopTimerMemory(times->preProcess, times->memoryPreProcess);

 domain->template getSolverAndKuc<T>(allOps, kelArray);  // PJSA 10-5-04 for freq sweep compatibility
 solver = allOps.sysSolver;                              // also need M, Muc, C and Cuc (in allOps)
 kuc = allOps.Kuc;
 kcc = allOps.Kcc;

 Rbm *rigidBodyModes = 0;
 int useProjector=domain->solInfo().filterFlags;
 if(useProjector) {
   std::cout << " ... RBMfilter Requested            ..." << std::endl;
   rigidBodyModes = domain->constructRbm();
   projector_prep(rigidBodyModes);
 }
 
 int useHzemFilter = domain->solInfo().hzemFilterFlag;
 if(useHzemFilter) {
   std::cout << " ... HZEMfilter Requested           ..." << std::endl;
   rigidBodyModes = domain->constructHzem();
   projector_prep(rigidBodyModes);
 }

 //ADDED FOR SLOSHING PROBLEM, EC, 20070723
 int useSlzemFilter = domain->solInfo().slzemFilterFlag;
 if(useSlzemFilter) {
   std::cout << " ... SLZEMfilter Requested           ..." << std::endl;
   std::cout << " ... THIS PART OF THE CODE HAS NOT BEEN DEBUGGED!! ..." << std::endl;
   rigidBodyModes = domain->constructSlzem();
   projector_prep(rigidBodyModes);
 }

 int useModeFilter = domain->solInfo().modeFilterFlag;
 if(useModeFilter) {
   std::cout << " ... MODEfilter requested          ..." << std::endl;
   eigmode_projector_prep();
 }

}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::scaleDisp(VectorType &sol)
{
  domain->scaleDisp(sol.data());
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::scaleInvDisp(VectorType &sol)
{
  domain->scaleInvDisp(sol.data());
}


template<class T, class VectorType, class SolverType>
SolverType *
SingleDomainStatic<T, VectorType, SolverType>::getSolver()
{
  return solver;
}

template<class T, class VectorType, class SolverType>
SingleDomainPostProcessor<T, VectorType, SolverType> *
SingleDomainStatic<T, VectorType, SolverType>::getPostProcessor()
{
 return new SingleDomainPostProcessor<T,VectorType,SolverType>(domain,bcx,times,solver,kuc,kcc);
}

// ... The next function is the only one needed for Statics and FAcoustics
template<class T, class VectorType, class SolverType>
void
SingleDomainPostProcessor<T, VectorType, SolverType>::staticOutput(VectorType &sol, 
                                                                   VectorType &force, 
                                                                   bool printTimers, int ndflag)
{
#ifdef DISTRIBUTED
 if(structCom->myID() != 0) return; // used for parallel mumps, only one process should write the output file
#endif
 startTimerMemory(times->output, times->memoryOutput);
 domain->template postProcessing<T>(sol,bcx,force,ndflag,0,0,0,kuc,kcc);
 stopTimerMemory(times->output, times->memoryOutput);

 long memoryUsed = 0;
 double solveTime  = 0.0;

 memoryUsed = solver->size();
 solveTime  = solver->getSolutionTime();
 if(printTimers) {
   times->printStaticTimers(solveTime, memoryUsed, domain);
   //geoSource->closeOutputFiles();
 }

 filePrint(stderr," --------------------------------------\n");
}

template<class T, class VectorType, class SolverType>
void
SingleDomainPostProcessor<T, VectorType, SolverType>::staticOutput(GeomState &geomState, 
                                                                   double x)
{
  times->output -= getTime();
  domain->template postProcessing<T>(geomState,x);
  times->output += getTime();

  MatrixTimers mt   = domain->getTimers();
  double solveTime  = solver->getSolutionTime();
  double memoryUsed = solver->getMemory();
}

template<class T, class VectorType, class SolverType>
void
SingleDomainPostProcessor<T, VectorType, SolverType>::staticOutput(GeomState *geomState, 
                                                                   double x)
{
  startTimerMemory(times->output, times->memoryOutput);
  domain->template postProcessing<T>(*geomState,x);
  stopTimerMemory(times->output, times->memoryOutput);

  MatrixTimers mt  = domain->getTimers();
  double solveTime = solver->getSolutionTime();
}


//------------------------------------------------------------------------------

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::getFreqSweepRHS(VectorType *rhs, VectorType **u, int k)
{ 
  startTimerMemory(times->formRhs, times->memoryRhs);
// RT: 11/17/2008
//  double omega2 = (domain->isFluidElement(0) && !domain->solInfo().isCoupled) 
//                   ? domain->getElementSet()[0]->helmCoef() : geoSource->shiftVal();
  double omega2 = geoSource->shiftVal();

  double omega = sqrt(omega2);

  VectorType *vec = new VectorType(solVecInfo());

  if (u==0) {
    for(int i=0; i<vec->size(); ++i)
      (*vec)[i] = 0;
  } else {
    for(int i=0; i<vec->size(); ++i)
      (*vec)[i] = double(k)*(double(k-1)*(*u[k-1])[i] + 2.0*omega*(*u[k])[i]);
    allOps.M->mult(vec->data(), rhs->data());

    if(allOps.C_deriv) {
      for(int j=0; j<=k-1; ++j) {
        if(allOps.C_deriv[k-j-1]) {
          double ckj = DCombination(k,j);
          for(int i=0; i<vec->size(); ++i) (*vec)[i] = -ckj*(*u[j+1])[i];
          allOps.C_deriv[k-j-1]->multAdd(vec->data(), rhs->data());
        }
      }
    }
  }
  delete vec;
  domain->template buildFreqSweepRHSForce<T>(*rhs, allOps.Muc, allOps.Cuc_deriv, k, omega);
  stopTimerMemory(times->formRhs, times->memoryRhs);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::getRHS(VectorType &rhs, double omega, double deltaomega)
{ 
  VectorType *vec = new VectorType(solVecInfo());
  domain->template buildRHSForce<T>(rhs, *vec, kuc, allOps.Muc, allOps.Cuc_deriv,
                                    omega, deltaomega);
  delete vec;
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::eigmode_projector_prep()
{
  if(Rmem) return; // already done it or requested some other filter

  // Read computed eigenvectors from file EIGENMODES
  // ======================================
  BinFileHandler modefile("EIGENMODES" ,"r");
  if(modefile.get_fileid() <= 0) { fprintf(stderr, " *** Error: Failed to open EIGENMODES file ***\n"); exit(-1); }

  modefile.read(&numR, 1);
  fprintf(stderr," ... Reading %d modes from EIGENMODES file ...\n", numR);

  int eigsize;
  modefile.read(&eigsize, 1);
  if(eigsize != solVecInfo()) { fprintf(stderr, " *** Error: Bad data in EIGENMODES file ***\n"); exit(-1); }

  Rmem = new double[numR*eigsize];
  for(int i = 0; i < numR; ++i)
    modefile.read(Rmem+i*eigsize, eigsize);

/*
  // Check if eigenvectors are M-orthonormal: Phi_i*M*Phi_i = 1 and Phi_i*M*Phi_j = 0
  // ======================================
  GenSparseMatrix<double> *M = (GenSparseMatrix<double> *) allOps.M; // XXXX won't work for complex
  double *tPhiM =  new double[numR*eigsize];
  for(int i = 0; i < numR; ++i)
    M->mult(Rmem+i*eigsize, tPhiM+i*eigsize);  // taking advantage of symmetry of M and computing
                                               // M*Phi_i instead of transpose(Phi_i)*M
  for(int i = 0; i < numR; ++i) {
    for(int j = 0; j < numR; ++j) {
      double PhiMPhi = 0;
      for(int k = 0; k < eigsize; ++k) {
        PhiMPhi += Rmem[k+i*eigsize]*tPhiM[k+j*eigsize];
      }
      fprintf(stderr, "Phi_%d*M*Phi_%d = %19.11e\n", i, j, PhiMPhi);
    }
    fprintf(stderr,"\n");
  }
*/
  // Build U_c(U_c^T*U_c)^{-1} for projector
  // ======================================
  StackFSFullMatrix Rt(numR, eigsize, Rmem);
  FullMatrix R = Rt.transpose();
  FullMatrix RtR(numR, numR);
  Rt.mult(R, RtR);
  FullMatrix RtRinverse = RtR.invert();

  X = new FullMatrix(eigsize, numR); // X = U(U^t*U)^{-1}
  R.mult(RtRinverse,(*X));
}

