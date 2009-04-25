#include <stdio.h>

#include <Interface.d/MpcLocal.h>
#include <Interface.d/FetiValues.h>
#include <Interface.d/FetiParams.h>

#include <Interface.d/SandiaDomain.h>
#include <Utils.d/Connectivity.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/Memory.h>
#include <Feti.d/Feti.h>
#include <Comm.d/Communicator.h>
#include <Driver.d/Communicator.h>

#define D50 1

// global variables
FetiParams *parameters = 0;
SysCom *scom = 0;
// extern FSCommunicator *communicator;
ThreadManager *threadManager;
extern FILE *debugFile;

const int print_timers = (1 << 0);
const int print_rbms   = (1 << 1);
const int prt_dbg      = (1 << 2);
const int prt_dbg2     = (1 << 3);
int problem_type = 0;
double initTime  = 0.0;
double solveTime = 0.0;
long initMemory = 0;
long solveMemory = 0;
int verboseFlag = 0;

bool makeNewInterface, makeNewLHS, makeNewSolver;
bool rebuildKcc, rebuildGtG, rebuildCCt;  
bool updateLHS, isFactored;

bool CU_DEBUG = false;

#ifdef TFLOP
extern "C" int heap_info(int*, int*, int*, int*);
#endif

template<class Scalar>
CU_Feti<Scalar>::CU_Feti(const double *X, const double *Y, const double *Z,
    int num_nodes, int *elem_ptr, int *elem_conn, int num_elems,
    int neq, const int *map, int ndofs_deleted, const int *bc_node_ids,
    int *bc_node_dof, MPI_Comm *salinas_communicator, int numads, int *procad,
    int *com_ptr, int *com_lis, unsigned short *dofmap_on_node,
    int *global_node_nums, FetiParams* params, int num_threads)
{
  // ECHO INPUT
/*
  int i;
  cerr << "num_nodes = " << num_nodes << endl;
  cerr << "XYZ = "; for(i=0;i<num_nodes;++i) cerr << X[i] << "," << Y[i] << "," << Z[i] << " "; cerr << endl;
  cerr << "num_elems = " << num_elems << endl;
  cerr << "elem_ptr = "; for(i=0;i<num_elems+1;++i) cerr << elem_ptr[i] << " "; cerr << endl;
  cerr << "elem_conn = "; for(i=0;i<elem_ptr[num_elems]; ++i) cerr << elem_conn[i] << " "; cerr << endl;
  cerr << "neq = " << neq << endl;
  cerr << "map = "; for(i=0; i<neq; ++i) cerr << map[i] << " "; cerr << endl;
  cerr << "ndofs_deleted = " << ndofs_deleted << endl;
  cerr << "bc_node_ids = "; for(i=0; i<ndofs_deleted; ++i) cerr << bc_node_ids[i] << " "; cerr << endl;
  cerr << "bc_node_dof = "; for(i=0; i<ndofs_deleted; ++i) cerr << bc_node_dof[i] << " "; cerr << endl;
  cerr << "numads = " << numads << endl;
  cerr << "procad = "; for(int i=0; i<numads; ++i) cerr << procad[i] << " "; cerr << endl;
  cerr << "com_ptr = "; for(int i=0; i<numads+1; ++i) cerr << com_ptr[i] << " "; cerr << endl;
  cerr << "com_lis = "; for(int i=0; i<com_ptr[numads]; ++i) cerr << com_lis[i] << " "; cerr << endl;
  cerr << "dofmap_on_node = "; for(int i=0;i<num_nodes;++i) cerr << dofmap_on_node[i] << " "; cerr << endl;
  for(i=0;i<num_nodes;++i) cerr << global_node_nums[i] << " "; cerr << endl;
*/
  //params->print(); 
  initTime -= getTime();
  initMemory -= memoryUsed();

  // CHECK INPUT
  if(num_elems <= 0)
    fprintf(stderr,"ERROR: This subdomain has %d elements\n", num_elems);
  if(num_nodes <= 0)
    fprintf(stdout,"ERROR: This subdomain has %d nodes\n",num_nodes);

  // INITIALIZE VARIABLES
  // 1. scom, structCom, threadManager
  if(salinas_communicator == 0) {
    printf("WARNING: salinas_communicator == 0 \n");
    scom = new SysCom();
  }
  else {
    scom = new SysCom(salinas_communicator);
  }
  structCom = scom;
  threadManager = new ThreadManager(num_threads);
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::CU_Feti() \n"); fflush(stdout); }

  // 2. geoSource
  geoSource = new GeoSource(0);
    
  // 3. domain
  domain = new Domain(0);

  // 4. parameters
  parameters = params;

  // 5. sdomain
  sdomain = new GenSandiaDomain<Scalar>(domain, X, Y, Z, num_nodes, elem_ptr, elem_conn, num_elems,
                                        neq, map, ndofs_deleted, bc_node_ids, bc_node_dof, numads, procad,
                                        com_ptr, com_lis, dofmap_on_node, global_node_nums, params);

  // 6. other
  fSolver = 0; distF = 0; distD = 0;
  verboseFlag = int(params->Verbose_flag() != 0);
  problem_type = params->Solution_type(); // 0 = static, 1 = eigen, 2 = dynamic, 3 = directfrf
                                          // 4 = eigen for modal acceleration (need to compute rbms)

  // MAKE THE SUBDOMAINS FOR THE FETI SOLVER
  sdomain->makeSubD();

  // MAKE THE SUBDOMAIN SHARED NODES CONNECTIVITIES
  sdomain->makeSComm();

  // APPLY THE CONSTRAINTS
  sdomain->constrain();

  makeNewInterface = makeNewLHS = makeNewSolver = true;
  rebuildKcc = rebuildGtG = rebuildCCt = false;
  updateLHS = false;
  isFactored = false;

  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::CU_Feti() \n"); fflush(stdout); }
}

template<class Scalar>
CU_Feti<Scalar>::~CU_Feti() 
{
  //if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::~CU_Feti() \n"); fflush(stdout); }
  if(sdomain) { delete sdomain; sdomain = 0; } 
  if(distF) { delete distF; distF = 0; }
  if(distD) { delete distD; distD = 0; }
  if(fSolver) { delete fSolver; fSolver = 0; }
  if(threadManager) { delete threadManager; threadManager = 0; }
  if(scom) { delete scom; scom = 0; }
  if(domain) { delete domain; domain = 0; }
  //if(CU_DEBUG && (myID == 0)) { printf("exiting CU_Feti::~CU_Feti() \n"); fflush(stdout); }
}

template<class Scalar>
int
CU_Feti<Scalar>::Factor()
{
  // this is used by the salinas eigen solver which require the RBMs from the factorization of Kcc or GtG
  // before the call to CU_Feti::Solve(...)
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::Factor(...) \n"); fflush(stdout); }

/* PJSA 6-15-07
  if(makeNewInterface || makeNewLHS) { printf("ERROR: bad calling sequence \n"); return 1; }
  // PJSA: as currently implemented, salinas must call zeroLHS() and addToLHS() after changing Contact/Tied/MPCs

  // decide here whether to make new solver due to change in problem_type
  // or rebuild Kcc due to change in corners/augmentation/shift/tolerances
  // etc for GtG & CCt
 
  if(problem_type == 3) {
    // always make new LHS & solver for directfrf since Kcc changes with different freq.
    makeNewLHS = true; // PJSA 4-28-05
    makeNewSolver = true; // always make new solver for directfrf since Kcc changes with different freq.
  }
*/
  int new_problem_type = parameters->Solution_type();
  if(new_problem_type != problem_type) { 
    if(scom->myID() == 0) {
      fprintf(stdout,"-----------------------------------------------"
                     "------------------------------\n");
      if(new_problem_type == 0)
        fprintf(stdout," Solution Type                     = Statics\n");
      else if(new_problem_type == 1 || new_problem_type == 4)
        fprintf(stdout," Solution Type                     = Eigen\n");
      else if(new_problem_type == 2)
        fprintf(stdout," Solution Type                     = Implicit Transient\n");
      else if(new_problem_type == 3)
        fprintf(stdout," Solution Type                     = Direct Frequency Response\n");
      fprintf(stdout,"-----------------------------------------------"
                     "------------------------------\n");
    }
    makeNewSolver = true; 
    problem_type = new_problem_type; 
  }

  if(makeNewSolver) {
    if(fSolver) delete fSolver;
    fSolver = sdomain->getSolver();  // this builds new Kcc, GtG and CCt and everything else
    if(problem_type == 3) makeNewLHS = true; // PJSA 6-15-07
    rebuildKcc = rebuildGtG = rebuildCCt = updateLHS = false;
    if(distF) { delete distF; distF = 0; }
    if(distD) { delete distD; distD = 0; }
  }
  else {
    if(updateLHS) { sdomain->updateSysMatrixInSolver(); updateLHS = false; }
    if(rebuildKcc) { fSolver->deleteKcc(); fSolver->makeKcc(); rebuildKcc = false; }
    if(rebuildGtG) { fSolver->deleteGtG(); fSolver->makeGtG(); rebuildGtG = false; }
    if(rebuildCCt) { fSolver->deleteCCt(); fSolver->buildCCt(); rebuildCCt = false; }
  }
  /*if(problem_type != 3)*/ makeNewSolver = false; // PJSA 6-15-07 need mrhs for freq sweep
 
  isFactored = true;
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::Factor(...) \n"); fflush(stdout); }
  return 0;
} 

template<class Scalar>
int
CU_Feti<Scalar>::Solve(Scalar *Force, Scalar *U_real, FetiValues &returnvals)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::Solve(...) \n"); fflush(stdout); }

  int err = 0;
  if(!isFactored) err = Factor();
  if(err) return 1;
  isFactored = false;  // always go into Factor() for every subsequent solve to check whether solver needs to be rebuilt or updated

  if(distF) distF->zero(); 
  else distF = new GenDistrVector<Scalar>(fSolver->localInfo());

  if(distD) distD->zero();
  else distD = new GenDistrVector<Scalar>(fSolver->localInfo());

  initMemory += memoryUsed();
  initTime   += getTime();
  solveTime   -= getTime();
  solveMemory -= memoryUsed();

  if(verboseFlag && (scom->myID() == 0)) fprintf(stdout," Making force,");
  sdomain->makeF(Force, distF);
  if(verboseFlag && (scom->myID() == 0)) fprintf(stdout," Begin Solving\n");
  fSolver->solve(*distF, *distD);
  if(verboseFlag && (scom->myID() == 0)) fprintf(stdout," Unifying Solution, ");
  sdomain->getD(U_real, distD);
  if(verboseFlag &&(scom->myID() == 0)) fprintf(stdout,"Returning to Salinas\n");

  solveTime   += getTime();
  solveMemory += memoryUsed();

  Timings &times    = fSolver->getTimers();
  int systemNumber  = times.numSystems-1;
  if(systemNumber < 0) systemNumber = 0;

  // Set return values to Salinas
  returnvals.number_iterations    = times.numIter; // PJSA DEBUG (for mrhs) 2-17-6
  returnvals.primal_residual      = times.iterations[systemNumber].finalPrimal;
  returnvals.dual_residual        = times.iterations[systemNumber].finalDual;

  if(systemNumber == 0) 
    returnvals.number_reorthog_used = (parameters->Max_N_orthog_vecs() > times.numIter) 
                                     ? times.numIter : parameters->Max_N_orthog_vecs();
  else {
    if(parameters->Max_N_orthog_vecs() > times.numIter) {
      returnvals.number_reorthog_used = times.iterations[systemNumber].numFetiIter;
    } 
    else {
      int i;
      returnvals.number_reorthog_used = parameters->Max_N_orthog_vecs();
      for(i=0; i<systemNumber; ++i)
        returnvals.number_reorthog_used -= times.iterations[i].numFetiIter;
      if(returnvals.number_reorthog_used < 0) returnvals.number_reorthog_used = 0;
    }
  }

  returnvals.memory_use_this_sub  = solveMemory + initMemory;
  returnvals.time_in_feti         = solveTime;

  if(times.converged)
    returnvals.feti_status = 0; // Success
  else
    returnvals.feti_status = 1; // Stagnation

  if(problem_type == 1 || problem_type == 4) { // EigenValue problem
    if(!parameters->Multiple_rhs()) fSolver->resetOrthoSet(); // PJSA 2-16-06
  }

  // By selecting a negative number of orthog vectors, FETI
  // will NOT use the MRHS acceleration techniques but will
  // apply full reorthogonalization for each rhs vector.
  if(parameters->Max_N_orthog_vecs() < 0)
    fSolver->resetOrthoSet();

  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::Solve(...) \n"); fflush(stdout); }
  return 0;
}

template<class Scalar>
int 
CU_Feti<Scalar>::addToLHS(const Scalar multiplier, const double *A_aa, const int *A_aj,
                          const int *A_ai, const bool isK)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::addToLHS(...) \n"); fflush(stdout); }

  // PJSA: should check here for change in LHS pattern & if detected set makeNewLHS = true
  if(makeNewInterface) { 
    sdomain->preProcess(); makeNewInterface = false; 
    if(distF) { delete distF; distF = 0; }
    if(distD) { delete distD; distF = 0; }
  }
  if(makeNewLHS) { 
    if(sdomain->numWaveDir() > 0) sdomain->setWaveNumbers(parameters); // PJSA 2-17-05
    sdomain->initSysMatrix(); 
    if(problem_type == 3) makeNewSolver = true;
    makeNewLHS = false; 
  }
  sdomain->addSysMatrix(A_aa, A_aj, A_ai, multiplier, isK);  
  rebuildKcc = true;
  updateLHS = true;
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::addToLHS(...) \n"); fflush(stdout); }
  return 0;
}

template<class Scalar>
void
CU_Feti<Scalar>::zeroLHS()
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::zeroLHS(...) \n"); fflush(stdout); }
  sdomain->zeroSysMatrix();
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::zeroLHS(...) \n"); fflush(stdout); }
}

template<class Scalar>
int 
CU_Feti<Scalar>::getNumRBMs() 
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::getNumRBMs(...) \n"); fflush(stdout); }
  if(!isFactored) Factor();
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::getNumRBMs(...) \n"); fflush(stdout); }
  return fSolver->numRBM();
}

template<class Scalar>
void 
CU_Feti<Scalar>::getRBMs(Scalar *rbms)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::getRBMs(...) \n"); fflush(stdout); }
  sdomain->getRBMs(rbms);
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::getRBMs(...) \n"); fflush(stdout); }
}

template<class Scalar>
int 
CU_Feti<Scalar>::NumberMpcConstraints()
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::NumberMpcConstraints() \n"); fflush(stdout); }
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::NumberMpcConstraints() \n"); fflush(stdout); }
  return sdomain->getNumberMpc();
}

template<class Scalar>
MpcLocal*
CU_Feti<Scalar>::GetLocalMpcs()
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::GetLocalMpcs() \n"); fflush(stdout); }
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::GetLocalMpcs() \n"); fflush(stdout); }
  return sdomain->getLocalMpcs();

}

template<class Scalar>
void 
CU_Feti<Scalar>::MpcForces(double *lambdas)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::MpcForces(...) \n"); fflush(stdout); }
  sdomain->mpcForces(lambdas);
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::MpcForces(...) \n"); fflush(stdout); }
}

template<class Scalar>
void 
CU_Feti<Scalar>::ConstraintProduct(int num_vect, const double* R[], Scalar** V, int trans)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::ConstraintProduct(...) \n"); fflush(stdout); }
  sdomain->constraintProduct(num_vect, R, V, trans);
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::ConstraintProduct(...) \n"); fflush(stdout); }
}

template<class Scalar>
void 
CU_Feti<Scalar>::setContactSurfaces(int numContactNeigh, int *contactNeighb, int *neighbPtr,
                                    enum face_type *ftype, int *face_ptr, int *face_conn,
                                    double (*face_normal)[3], double (*gap_vector)[3], int *glbl_id,
                                    double *normal_tol, double *tangent_tol)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::setContactSurfaces(...) \n"); fflush(stdout); }
  sdomain->setSurfaceInteractions(1, numContactNeigh, contactNeighb, neighbPtr,
                                  ftype, face_ptr, face_conn, face_normal, gap_vector, glbl_id,
                                  normal_tol, tangent_tol);

  makeNewInterface = true; rebuildGtG = true; rebuildCCt = true;
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::setContactSurfaces(...) \n"); fflush(stdout); }
}

template<class Scalar>
void 
CU_Feti<Scalar>::setTiedSurfaces(int numTiedNeigh, int *tiedNeighb, int *neighbPtr,
                                 enum face_type *ftype, int *face_ptr, int *face_conn,
                                 double (*face_normal)[3], double (*gap_vector)[3], int *glbl_id,
                                 double *normal_tol, double *tangent_tol)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::setTiedSurfaces(...) \n"); fflush(stdout); }
  sdomain->setSurfaceInteractions(0, numTiedNeigh, tiedNeighb, neighbPtr, ftype,
                                  face_ptr, face_conn, face_normal, gap_vector, glbl_id,
                                  normal_tol, tangent_tol);
  if(parameters->Mpc_method().method == MPC_Method::Dual) { 
    makeNewInterface = true; rebuildGtG = true; rebuildCCt = true; 
  }
  else rebuildKcc = true; 
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::setTiedSurfaces(...) \n"); fflush(stdout); }
}

template<class Scalar>
void
CU_Feti<Scalar>::setMPCs(int NumberMpc, MpcLocal* MpcVector)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::setMPCs(...) \n"); fflush(stdout); }
  sdomain->setMPCs(NumberMpc, MpcVector);
  if(parameters->Mpc_method().method == MPC_Method::Dual) { 
    makeNewInterface = true; rebuildGtG = true; rebuildCCt = true; 
  } 
  else rebuildKcc = true; 
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::setMPCs(...) \n"); fflush(stdout); }
}

template<class Scalar>
int
CU_Feti<Scalar>::updateLocalMpcs(MpcLocal** MpcVector)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti::updateLocalMpcs(...) \n"); fflush(stdout); }
  
  int ret = (sdomain->updateLocalMpcs(MpcVector));
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::updateLocalMpcs(...) \n"); fflush(stdout); }
  return ret;
}

template<class Scalar>
int
CU_Feti<Scalar>::setContactNormal(const double (*face_normal)[3])
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti:setContactNormal(...)\n"); fflush(stdout); }
  sdomain->setContactNormal(face_normal);
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::setContactNormal(...)\n"); fflush(stdout); }
  return(0);
}

template<class Scalar>
int
CU_Feti<Scalar>::setContactGap(/*const double (*gap_vectors)[3]*/ const double *disp)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti:setContactGap(...)\n"); fflush(stdout); }
  //sdomain->setContactGap(gap_vectors); // PJSA 10-27-05
  sdomain->setContactGap(disp); // PJSA 9-18-2007
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::setContactGap(...)\n"); fflush(stdout); }
  return(0);
}

template<class Scalar>
int
CU_Feti<Scalar>::updateContactDisp(const double *disp_updated)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti:updateContactDisp(...)\n"); fflush(stdout); }
  sdomain->updateContactGap(disp_updated);
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::updateContactDisp(...)\n"); fflush(stdout); }
  return(0);
}


template<class Scalar>
int
CU_Feti<Scalar>::getContactForce(double *contact_forces)
{
  if(CU_DEBUG && (scom->myID() == 0)) { printf("entering CU_Feti:getContactForce(...)\n"); fflush(stdout); }
  sdomain->getContactForces(contact_forces);
  if(CU_DEBUG && (scom->myID() == 0)) { printf("exiting CU_Feti::getContactForce(...)\n"); fflush(stdout); }
  return(0);
}

