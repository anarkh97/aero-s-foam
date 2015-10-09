#include <Utils.d/dbg_alloca.h>
#include <iostream>
#include <cstdio>
#include <cmath>

#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Corotational.d/GeomState.h>
#include <Solvers.d/Rbm.h>

#include <Math.d/FullMatrix.h>
#include <Sfem.d/Sfem.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/MathUtils.h>

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#endif

extern Sfem *sfem;
typedef FSFullMatrix FullMatrix;
template <class Scalar> class GenVector;
typedef GenVector<DComplex> ComplexVector;

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::projector_prep(Rbm *rbms)
{
 numR = rbms->numRBM(); 
 if (!numR) return;

 // KHP: store this pointer to the RBMs to use in the actual
 //      projection step 
 int useProjector = domain->solInfo().filterFlags;
 if(useProjector) {
   filePrint(stderr," ... Building the RBM Projector     ...\n");
   filePrint(stderr," ... Number of RBMs = %-4d          ...\n",numR);
 }

 int useHzemFilter = domain->solInfo().hzemFilterFlag;
 if(useHzemFilter) {
   filePrint(stderr," ... Building the HZEM Projector    ...\n");
   filePrint(stderr," ... Number of HZEMs = %-4d         ...\n",numR);
 }
 int ndof = solVecInfo();
 Rmem = new double[numR*ndof];
 if (useProjector) rbms->getRBMs(Rmem);
 if (useHzemFilter) {
   for(int n=0; n<ndof; ++n) Rmem[n] = 1.;
 }

 int useSlzemFilter = domain->solInfo().slzemFilterFlag;
 if(useSlzemFilter) {
   filePrint(stderr," ... Building the SLZEM Projector    ...\n");
   filePrint(stderr," ... Number of SLZEMs = %-4d         ...\n",numR);
 }

 if (useSlzemFilter) {
   for(int n=0; n<ndof; ++n) Rmem[n] = 1.;
 }

 StackFSFullMatrix Rt(numR, ndof, Rmem);

 FullMatrix R = Rt.transpose();

 FullMatrix RtR(numR,numR);
 Rt.mult(R,RtR);

 FullMatrix RtRinverse = RtR.invert();

 X = new FullMatrix(ndof,numR);
 R.mult(RtRinverse,(*X));
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
 if(domain->solInfo().loadcases.size() > 0)
   filePrint(stderr," ... Building the Force (Case %2d)   ...\n", domain->solInfo().loadcases.front());
 else
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
 filePrint(stderr," ... Building the Force   (inpc)    ...\n");

 // ... BUILD THE RHS FORCE (external + gravity + nonhomogeneous)
 startTimerMemory(times->formRhs, times->memoryRhs);
 rhs.zero();
 for (int i=0; i<(*allOps.rhs_inpc).size(); ++i) rhs[i] = (*allOps.rhs_inpc)[i]; 

 int useProjector=domain->solInfo().filterFlags;
 if(useProjector)
   project(rhs);

 stopTimerMemory(times->formRhs, times->memoryRhs);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::preProcessSA()
{
 domain->buildPreSensitivities<T>(allSens, bcx);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::postProcessSA(GenVector<T> &sol)
{
 domain->buildPostSensitivities<T>(allOps.sysSolver, allOps.K, allOps.spm, allSens, &sol, bcx);
 domain->sensitivityPostProcessing(allSens, &sol, bcx);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::preProcess()
{
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
   FullSquareMatrix *geomKelArray=0, *dummy=0;
   domain->computeGeometricPreStress(allCorot, geomState, kelArray, times,
                                     geomKelArray, dummy);
 }

 stopTimerMemory(times->preProcess, times->memoryPreProcess);

 if(!rigidBodyModes) {
   int useProjector = domain->solInfo().filterFlags;
   int useHzemFilter = domain->solInfo().hzemFilterFlag;
   int useSlzemFilter = domain->solInfo().slzemFilterFlag;

   // ... Construct geometric rigid body modes if necessary
   if(useProjector || domain->solInfo().rbmflg) {
     rigidBodyModes = domain->constructRbm();
     if(useProjector) {
       std::cout << " ... RBM Filter Requested           ..." << std::endl;
       projector_prep(rigidBodyModes);
     }
   }

   // ... Construct "thermal rigid body mode" if necessary
   else if(useHzemFilter || domain->solInfo().hzemFlag) {
     rigidBodyModes = domain->constructHzem();
     if(useHzemFilter) {
       std::cout << " ... HZEM Filter Requested          ..." << std::endl;
       projector_prep(rigidBodyModes);
     }
   }

   // ... Construct "sloshing rigid body mode" if necessary
   else if(useSlzemFilter || domain->solInfo().slzemFlag) {
     rigidBodyModes = domain->constructSlzem();
     if(useSlzemFilter) {
       std::cout << " ... SLZEM Filter Requested         ..." << std::endl;
       projector_prep(rigidBodyModes);
     }
   }
 }

 if(domain->solInfo().rbmflg || domain->solInfo().hzemFlag || domain->solInfo().slzemFlag) {
   domain->template getSolverAndKuc<T>(allOps, kelArray, rigidBodyModes);
 }
 else {
   domain->template getSolverAndKuc<T>(allOps, kelArray, (Rbm*)NULL); 
 }
 solver = allOps.sysSolver;
 kuc = allOps.Kuc;
 kcc = allOps.Kcc;
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
void
SingleDomainStatic<T, VectorType, SolverType>::scaleDisp(VectorType &sol, double alpha)
{
  domain->scaleDisp(sol.data(), alpha);
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
 double solveTime = 0.0;

 memoryUsed = solver->size();
 solveTime  = solver->getSolutionTime();
 if(printTimers) {
   times->printStaticTimers(solveTime, memoryUsed, domain);
 }

 if(ndflag <= 1) filePrint(stderr," --------------------------------------\n");
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
    if(allOps.K_deriv) {
      for(int j=0; j<=k-1; ++j) {
        if(allOps.K_deriv[k-j]) {
          double ckj = DCombination(k,j);
          for(int i=0; i<vec->size(); ++i) (*vec)[i] = -ckj*(*u[j+1])[i];
          allOps.K_deriv[k-j]->multAdd(vec->data(), rhs->data());
        }
      }
    }
  }
  delete vec;
  domain->template buildFreqSweepRHSForce<T>(*rhs, allOps.Muc, allOps.Cuc_deriv,allOps.Kuc_deriv, k, omega);
  stopTimerMemory(times->formRhs, times->memoryRhs);
}

template<class T, class VectorType, class SolverType>
void
SingleDomainStatic<T, VectorType, SolverType>::getRHS(VectorType &rhs, double omega, double deltaomega)
{ 
  VectorType *vec = new VectorType(solVecInfo());
  domain->template buildRHSForce<T>(rhs, *vec, kuc, allOps.Muc,
                                    allOps.Cuc_deriv, allOps.Kuc_deriv,
                                    allOps.Kuc_arubber_l,allOps.Kuc_arubber_m,
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

