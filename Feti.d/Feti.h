#ifndef _FETI_H_
#define _FETI_H_

#include <Feti.d/DistrVector.h>
#include <Threads.d/Paral.h>
#include <Timers.d/Timing.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Communicator.h>
#include <Utils.d/MyComplex.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/SComm.h>
#include <Utils.d/dbg_alloca.h>
#include <Feti.d/DistrVectorSet.h>
#include <Solvers.d/ParallelSolver.h>

template <class Scalar> class GenFetiOp;
typedef GenFetiOp<double> FetiOp;
template <class Scalar> class GenFetiOpControler;
typedef GenFetiOpControler<double> FetiOpControler;
template <class Scalar> class GenFetiSolver;
typedef GenFetiSolver<double> FetiSolver;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
class IntFullM;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
class OrthoSet;
template <class Scalar> class GenGMRESOrthoSet;
template <class Scalar> class GenGCROrthoSet;
template <class Scalar> class GenCGOrthoSet;
template <class Scalar> class GenFetiWorkSpace;
typedef GenFetiWorkSpace<double> FetiWorkSpace;
template <class Scalar> class GenSymFullMatrix;
typedef GenSymFullMatrix<double> SymFullMatrix;
template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
class FetiInfo;
class DistrGeomState;
template <class Scalar> class GenBigMatrix;
typedef GenBigMatrix<double> BigMatrix;
template <class Scalar> class GenSkyMatrix;
typedef GenSkyMatrix<double> SkyMatrix;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
class Rbm;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
template <class Scalar> class GenBlockSky;
typedef GenBlockSky<double> BlockSky;
template <class Scalar> class GenBLKSparseMatrix;
typedef GenBLKSparseMatrix<double> BLKSparseMatrix;
template <class Type> class FSCommPattern;
template <class Scalar> class CCtSolver;

// The FetiSolver class is too big. Future versions need it to be leaner
template <class Scalar>
class GenFetiSolver  : public GenParallelSolver<Scalar>
{
  protected:
    GenSubDomain<Scalar> **sd;
    int nsub;
    FetiInfo *fetiInfo;
    FSCommunicator *fetiCom; // PJSA
    int myCPU, numCPUs; // PJSA
    Connectivity *subToSub, *mpcToSub, *mpcToSub_primal;
    Connectivity *edgeToSub, *subToEdge;
    Connectivity *coarseConnect;  // first level coarse prob. connectivity
    compStruct *renum;
    compStruct renumber;
    SimpleNumberer *eqNums;
    SimpleNumberer *PFcNums;
    int gOffset;
    int mOffset;
    GenSparseMatrix<Scalar> *singleCoarse;
    GenSolver<Scalar> *singleCoarseSolver;
    double epsilon2;
    int maxiter;
    GenCGOrthoSet<Scalar> *oSetCG;
    GenGMRESOrthoSet<Scalar> *oSetGMRES;
    GenGCROrthoSet<Scalar> *oSetGCR;
    int numP;
    int numrbms, halfSize; // number of  rbms, half interface size
    int glNumRBM;
    DistrInfo internalDI, interface;
    FSCommPattern<Scalar> *vPat;
    FSCommPattern<Scalar> *rbmPat;
    FSCommPattern<int> *sPat;
    int *glSubToLoc;
    int crns; // sum of corner lambdas of all subdomains
    TaskDescr **fetiTasks;
    GenFetiOp<Scalar> **fetiOps;
    GenFetiOpControler<Scalar> *opControl;
    GenSkyMatrix<Scalar> *GtGSkyMatrix; // used for nonlinear rebuilding of GtG
    GenSolver<Scalar> *GtGsolver;
    GenBigMatrix<Scalar> *PCtFPC;
    GenFetiWorkSpace<Scalar> *wksp;
    int QGisLocal;  // Whether or not QG is local
    int isDynamic;  // Whether or not we are in dynamics
    int isFeti2;    // Whether or not we are using FETI2
    int numSystems; // Nonlinear additions
    Timings times;
    GenSymFullMatrix<Scalar> *GtQGs;
    GenFullM<Scalar> *GtFCs;
    Scalar *gtqglocal;
    int numNodes;
    Connectivity *cpuToSub;
    int glNumMpc, glNumMpc_primal;

    void makeGtG();
    void makeDistGtG(int *glSubToLocal);
    void assembleDistGtQGs(int i, int *);
#ifdef BOOL_NOT_DEFINED
    int isLowestLocalNeighbor(int subI, int subJ);
#else
    bool isLowestLocalNeighbor(int subI, int subJ);
#endif
    void addNonLocalGtQG(int subI, int subJ);
    void addNonLocalGContrib(int subI, int subJ);
    void addNonLocalCContrib(int subI, int subJ);
    void getNonLocalGtQMult(int myNum, int neighbN, Scalar *va,
                            GenDistrVector<Scalar> *dv);
    void getNonLocalFCtMult(int myNum, int neighbN, Scalar *va,
                            GenDistrVector<Scalar> *dv);
    void getNonLocalSubAlphaGtQ(int subI, int subJ, Scalar *va, 
                                GenDistrVector<Scalar> *dv);
    void getNonLocalGtQMult(int subI, int subJ);
    void getGtQMult(int iSub, Scalar *, GenDistrVector<Scalar> *);
    void getFCMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv);
    void reBuildGtG();
    void makePCtFPC();
    void reBuildPCtFPC();
    void makelocalFcoarse();
    int  collectIntGlobalSum();
    void subRgcMult(int i, int nThreads, GenVector<Scalar>* alpha, GenVector<Scalar>* result);
    void subRgcTransMult(int i, int nThreads, GenVector<Scalar>* alpha, GenVector<Scalar>* result);
    void fSend(int i,  GenDistrVector<Scalar> &);
    void fScale(int i, GenDistrVector<Scalar> &);
    void fSplit(int iSub, GenDistrVector<Scalar> &force);
    void fSendCoupled(int iSub, GenDistrVector<Scalar> &force, GenDistrVector<Scalar> &fw);
    void fScaleCoupled(int iSub, GenDistrVector<Scalar> &force, GenDistrVector<Scalar> &fw);

    // Routines to create the single coarse problem
    void makeGandFG();
    void makeSingleCoarse();
    void singleCoarseAssembly();
    void singleCoarseAssembleG(int isub);
    void singleCoarseAssembleMPCs(int iSub);
    Connectivity * getCoarseToSubConnect();

    Connectivity * makeSingleConnect(Connectivity *coarseConnect, 
                                     Connectivity *coarseToSub,
                                     Connectivity *subToCoarse, int gOffset);
    void singleCoarseSolve(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &);
    void setAllOffsets(int iSub, int gOffset);
    void getSRMult(int iSub, GenDistrVector<Scalar> *r, GenDistrVector<Scalar> *lambda, 
                   Scalar *alpha);
    void getSGtMult(int iSub, GenDistrVector<Scalar> *r, Scalar *alpha);
    void getSCtMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv);
    void getSQtMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv);
    void getFGMult(int iSub, GenDistrVector<Scalar> *r, Scalar *alpha);
    void addG(int iSub, GenDistrVector<Scalar> *r, Scalar *sv);
    void addSG(int iSub, GenDistrVector<Scalar> *r, Scalar *sv);
    void addGs(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &w, GenVector<Scalar> &);
    void singlePr(GenDistrVector<Scalar> &y, GenDistrVector<Scalar> &p, GenVector<Scalar> &);
    void addMpcRhs(int iMPC, Scalar *singleC);
    void preprocessMPCs(int iSub);
    // corner preprocessing
    void preProcessCorners();
    
    void addC(int iSub, GenDistrVector<Scalar> *lambda, Scalar *sv);
    void getQtKpBMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv);
    void getNeighbFGs(int iSub);

  private:
    void initialize();

  public:
    GenFetiSolver() { initialize(); };
    GenFetiSolver(int nsub, GenSubDomain<Scalar> **, Connectivity *,
                  FetiInfo *finfo, FSCommunicator *fetiCom, int *glToLoc, 
                  Connectivity *mpcToSub, Connectivity *cpuToSub,
                  GenSolver<Scalar> **sysMatrices=0, GenSparseMatrix<Scalar> **sysMat = 0, 
                  Rbm **_rbms=0);
    GenFetiSolver(int _nsub, int _numThreads) : internalDI(_nsub), 
                  interface(_nsub), times(_numThreads,_nsub) { initialize(); }
    virtual ~GenFetiSolver();

    void sendDeltaF(int iSub, GenDistrVector<Scalar>& deltaF);
    void normDeltaF(int iSub, double * subDots, GenDistrVector<Scalar>* deltaF);
    void assembleMatrices(int isub);
    void factorMatrices(int isub);
    void constructMatrices(int isub);
    void sendScale(int isub);
    void collectScale(int isub);
    void makeSubdomainStaticLoad(int isub,GenDistrVector<Scalar>& f,DistrGeomState *gs=0);//HB: add DistrGeomState
    void makeSubdomainStaticLoadGalPr(int isub,GenDistrVector<Scalar>& f, GenDistrVector<Scalar>& tmp, double *, DistrGeomState *gs=0);//HB: add DistrGeomState
                                                                                          //follower load (pressure)
    void getErrorEstimator(int iSub, GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &es);
    void interfaceDiff(int iSub, GenDistrVector<Scalar> &v);
    void rebuildInterface(int iSub, GenDistrVector<Scalar> &v);
    void multKbb(int iSub, GenDistrVector<Scalar>& v, GenDistrVector<Scalar> &interfvec,
		 GenDistrVector<Scalar>& deltaU, GenDistrVector<Scalar> &deltaF, bool &errorFlag);
    void localSolve(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
                    GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4);
    void localSolve2(int iSub, GenDistrVector<Scalar> *v1, GenDistrVector<Scalar> *v2,
                     GenVector<Scalar> *beta, GenDistrVector<Scalar> *v3, 
                     GenDistrVector<Scalar> *v4);
    void interfSend(int iSub, GenDistrVector<Scalar> &dv1);
    //void interfDiff(int iSub, GenDistrVector<Scalar> &dv1);
    void interfDiffAndDot(int iSub, GenDistrVector<Scalar> &dv1, GenDistrVector<Scalar> &dv2);
    void getRMult(int iSub, GenDistrVector<Scalar> *localvec, Scalar *alpha);
    void getGtMult(int iSub, GenDistrVector<Scalar> *localvec, Scalar *alpha);
    void addRP(int iSub, GenDistrVector<Scalar> * localvec, Scalar *alpha);
    void addRS(int iSub, GenDistrVector<Scalar> * localvec, Scalar *alpha);
    virtual void solve(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &);
    void distributeForce(GenDistrVector<Scalar> &force);
    void distributeForce(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fw); 
    void reSolve(GenDistrVector<Scalar> &);
    Timings& getTimers() { return times; }
    void makeRbmPat();
    void makeSingleIntPat();

    // NonLinear Functions
    virtual void reBuild(FullSquareMatrix **kel, DistrGeomState& gs, int iter=0,
                         int step = 1);
    void reBuildMatrices(int isub, FullSquareMatrix **kel);
    void reBuildErrorEstimator(GenFullSquareMatrix<Scalar> **kel);
    void subdomainReBuild(int isub, FullSquareMatrix **kel, DistrGeomState *gs);
    int nlPreCondition(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &);
    double getSolutionTime() { return times.solve; }
    void reSendInterfaceRBM(int iSub);
    void reGetNeighbQGs(int iSub);
    void reMakeLocalFcoarse();
    void reGetNeighbFC(int iSub);
    void reComputeFiBC(int iSub);
    void makeCompatible(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &deltaU);
    Scalar localSolveAndJump(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &,
                             GenDistrVector<Scalar> &, GenDistrVector<Scalar> &);
    Scalar localSolveAndJump(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &,
                             GenVector<Scalar> &, GenDistrVector<Scalar> &, 
                             GenDistrVector<Scalar> &);
    void project(GenDistrVector<Scalar> &, GenVector<Scalar> &, GenDistrVector<Scalar> &, int isDirect = 1);
    void tProject(GenDistrVector<Scalar> &, GenVector<Scalar> &, GenDistrVector<Scalar> &, int isDirect = 1);
    void computeDynamL0(GenDistrVector<Scalar> &, GenVector<Scalar> &, 
                        GenDistrVector<Scalar> &, GenDistrVector<Scalar> &);
    void computeL0(GenDistrVector<Scalar> &, GenVector<Scalar> &, GenDistrVector<Scalar> &);
    void computeL0(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &, GenDistrVector<Scalar> &, 
                   GenVector<Scalar> &, GenDistrVector<Scalar> &);
    void addR(GenDistrVector<Scalar> &, GenVector<Scalar> &);
    void addRSingle(GenDistrVector<Scalar> &, GenVector<Scalar> &);
    void outputPrimalResidual(int iter, GenDistrVector<Scalar> &deltaF);

    // preconditioner returns the error estimator
    double preCondition(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Pv, bool errorFlag = true);

    void orthoAdd(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &, Scalar pFp = 0);

    void resetOrthoSet();
    void orthoAddCG(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar pFp);
    void orthogonalize(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &); 

    // GMRES functions
    void initGMRES(GenDistrVector<Scalar> &p);
    double orthoAddGMRES(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp);
    void GMRESSolution(GenDistrVector<Scalar> &p);

    // GCR functions
    int predictGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &lambda0);
    void orthogonalizeGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fr,
                          GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp);
    void orthoAddGCR(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar FpFp);

    int predict(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &);

    DistrInfo &interfInfo() { return interface; }
    DistrInfo &localInfo()  { return internalDI; }
    GenDistrVector<Scalar> & getLambda() { return wksp->ret_lambda(); }
    void makeStaticLoad(GenDistrVector<Scalar> &, DistrGeomState *gs=0); //HB: add DistrGeomState for 
                                                                         //follower forces (i.e. pressure)
    void makeStaticLoad(GenDistrVector<Scalar> &, double, double, DistrGeomState *gs=0); //HB: add DistrGeomState for 
    int neq() { return internalDI.len; }

    void updateFeti2lambda(GenDistrVector<Scalar> &, GenDistrVector<Scalar>&, 
                           GenDistrVector<Scalar> &,
                           GenVector<Scalar>&, GenVector<Scalar>&);
    void updateFeti2y(GenDistrVector<Scalar> &, GenDistrVector<Scalar>&,
                      GenVector<Scalar>&, GenVector<Scalar>&);
    void computeRgcTransMult(GenVector<Scalar>&, GenVector<Scalar>&);
    void computeRgcMult(GenVector<Scalar>&, GenVector<Scalar>&);
    void updateFeti2Vector(GenDistrVector<Scalar> &, GenVector<Scalar>&, GenVector<Scalar>&);

    void addAllFcoarse(GenFullM<Scalar> &);
    void finishRgc(int, int);

    void getCtMult(GenDistrVector<Scalar> &w, GenVector<Scalar> &gamma);

    void scatterHalfInterface(int iSub, Scalar *v1, GenDistrVector<Scalar> *v2);
    void gatherHalfInterface(int iSub, GenDistrVector<Scalar> *v1, 
                             GenDistrVector<Scalar> *v2,
                             Scalar *v3, Scalar *v4);

    void setAndStoreInfo(int iter, double finalPrimal2, double finalDual2 );
    void findProfileSizes(int iSub, int *subSizes);

    // For eigen problem 
    virtual int numRBM();
    virtual void getRBMs(GenDistrVectorSet<Scalar> &); //CBM
    virtual void getRBMs(Scalar *);
    void Ksolve(int iSub, GenStackDistVector<Scalar> &R);

    int halfOffset(int iSub)      { return fetiOps[iSub]->halfOffset; }
    int numNeighbor(int iSub); // { return sd[iSub]->getSComm()->numNeighb; }
    Scalar *interfaceBuffer(int iSub) { return fetiOps[iSub]->interfBuff; }
    virtual void clean_up();
    virtual double getFNormSq(GenDistrVector<Scalar> &f);

    virtual void getLocalMpcForces(int iSub, double *mpcLambda) { };  // only implemented for DP

    // for helmholtz freqency sweep
    void getFreqSweepRHS(GenDistrVector<Scalar> *rhs, GenDistrVector<Scalar> **u, int k);
    void multM(int iSub, GenDistrVector<Scalar> *rhs,  GenDistrVector<Scalar> **u, int k);
    void multMCoupled1(int iSub, GenDistrVector<Scalar> *rhs, GenDistrVector<Scalar> **sol_prev, int k);
    void multMCoupled2(int iSub, GenDistrVector<Scalar> *rhs);
    void pade(GenDistrVector<Scalar> *sol, GenDistrVector<Scalar> **u, double *h, double x);
    void subPade(int iSub, GenDistrVector<Scalar> *sol, GenDistrVector<Scalar> **u, double *h, double x);
    virtual void rebuildLHSfreq();  // for multiple LHS frequency sweep

    GenSolver<Scalar> * newSolver(int type, Connectivity *con, EqNumberer *nums, double tol, GenSparseMatrix<Scalar> *&sparse); // PJSA 2-23-2007

    virtual void initNewton() { cerr << "initNewton undefined for Feti-1 \n"; }

  protected:
    FSCommPattern<Scalar> *wiPat;
};

template<class Scalar>
class GenFetiDPSolver : public GenFetiSolver<Scalar> 
{
    DistrInfo internalR, internalC, internalWI;
    GenSolver<Scalar>       *KccSolver;
    GenSparseMatrix<Scalar> *KccSparse;
    int glNumCorners;
    Connectivity *cornerToSub;
    DofSetArray *cornerEqs;
    ConstrainedDSA *c_cornerEqs;
    int mpcOffset; // mpc equation offset for coarse grid
    void initialize();
    bool computeRbms;
    enum StepType { CG, PROPORTIONING, EXPANSION };
    double delta_lag, delta_alphabar, alphabar_prev, alphabar_lower, alphabar_upper, delta_lag_max, alphabar_max;

 public:
    GenFetiDPSolver(int nsub, GenSubDomain<Scalar> **sd, Connectivity *subToSub,
                    FetiInfo *finfo, FSCommunicator *fetiCom, int *glToLoc, Connectivity *mpcToSub, Connectivity *mpcToSub_primal,
                    Connectivity *mpcToMpc, Connectivity *mpcToCpu, Connectivity *cpuToSub, 
                    Connectivity *bodyToSub = 0, GenSolver<Scalar> **sysMatrices = 0,
                    GenSparseMatrix<Scalar> **sysMat = 0, Rbm **rbms = 0, int sandiaFlag = 0,
                    bool computeRbms = true);
    virtual ~GenFetiDPSolver();

    int ngrbm;
    void makeKbb(int iSub);
    void makeFc(int iSub, GenDistrVector<Scalar> &fr, /*GenVector<Scalar> &fc,*/ GenDistrVector<Scalar> &lambda);
    void makeFcB(int iSub, GenDistrVector<Scalar> &bf);
    void reBuild(FullSquareMatrix **kel, DistrGeomState& gs, int iter=0,
                 int step = 1);
    void reBuildMatrices(int isub, GenFullSquareMatrix<Scalar> **kel);
    void subdomainReBuild(int isub, FullSquareMatrix **kel, DistrGeomState *gs);    
    void KrrReSolve(int iSub, GenDistrVector<Scalar> &ur);
    void makeKcc();
    void deleteKcc() { /* not implemented */ };
    void reBuildKcc(); // this is used for nonlinear
    void setSysMatrices(GenSolver<Scalar> **sysMatrices, GenSparseMatrix<Scalar> **sysMat) { /* not implemented */ };
    void assembleFcStar(GenVector<Scalar> &FcStar);
    void mergeSolution(GenDistrVector<Scalar> &ur, GenVector<Scalar> &uc, GenDistrVector<Scalar> &u,
                       GenDistrVector<Scalar> &lambda);
    void mergeUr(int iSub, GenDistrVector<Scalar> &ur, GenVector<Scalar> &uc, GenDistrVector<Scalar> &u,
                 GenDistrVector<Scalar> &lambda);
    void multKrc(int iSub, GenDistrVector<Scalar> &fr, GenVector<Scalar> &uc);
    void extractFr(int iSub, GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fr);
    void extractFc(int iSub, GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fc);
    void extractFw(int iSub, GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fw);
    void getFc(GenDistrVector<Scalar> &f, GenVector<Scalar> &fc);
    void makeEdgeConnectivity();
    void countEdges(int iSub, int *edges);
    void numberEdges(int iSub, int *eP, int *ep2, int *edges, FSCommPattern<int> *sPat);
    void receiveNeighbEdgeNums(int iSub, int *eP, int *edges, FSCommPattern<int> *sPat);
    void factorLocalMatrices(int isub);
    double extractForceVectors(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fr, GenVector<Scalar> &fc, GenDistrVector<Scalar> &fw);
    void printSummary(int iter);
    double getFNormSq(GenDistrVector<Scalar> &f);
    void getRBMs(Scalar *);
    void getRBMs(GenDistrVectorSet<Scalar> &); //CBM
    void getGlobalRBM(int iSub, int &iRBM, GenDistrVector<Scalar> &R); //CBM
    int numRBM();
    void clean_up();
    void subdomainSolve(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
                        GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4, 
                        GenVector<Scalar> &v5);
    void subdomainSolveCoupled1(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
                        GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4,
                        GenVector<Scalar> &v5);
    void subdomainSolveCoupled2(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
                        GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4,
                        GenVector<Scalar> &v5, GenDistrVector<Scalar> &fw);
    void subdomainSolveCoupled2(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
                        GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4,
                        GenVector<Scalar> &v5);
    void localSolveAndJump(GenDistrVector<Scalar> &fr, GenDistrVector<Scalar> &lambda,
                           GenDistrVector<Scalar> &ur, GenVector<Scalar> &fc,
                           GenVector<Scalar> &uc, GenDistrVector<Scalar> &r,
                           GenDistrVector<Scalar> &fw);
    Scalar localSolveAndJump(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &dur, 
                             GenVector<Scalar> &duc, GenDistrVector<Scalar> &Fp);
    void solve(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
    void solveCG(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
    void solveGMRES(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
    void solveGCR(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
    void init_solve();
    void init_iter();
    void init_linesearch();
    void orthoAdd(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar pFp);
    Scalar lagrangian(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &d);
    double preCondition(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Pv, bool errorFlag = true);
    // augmented Lagrangian
    void solveCG_augLag(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
    void localSolveAndJump(GenDistrVector<Scalar> &fr, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &ur, GenVector<Scalar> &fc,
                           GenVector<Scalar> &uc, GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &fw, double rho, GenVector<Scalar> &mu);
    Scalar localSolveAndJump(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &dur, GenVector<Scalar> &duc,
                             GenDistrVector<Scalar> &Fp, double rho);
    Scalar lagrangian(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &r, double rho, GenVector<Scalar> &mu);
    double preCondition(GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Pv, double rho, bool errorFlag = false);
    double primalError();
    double dualError(GenDistrVector<Scalar> &r, bool &proportional);

    // note: these functions could belong to a FetiInterfaceVector class (GenDistrVector sub-class with pointer to **sub)
    Scalar Max(GenDistrVector<Scalar> &v, int flag);
    void subMax(int iSub, GenDistrVector<Scalar> &v, Scalar *partial, int flag);
    int Equal(GenDistrVector<Scalar> &v, Scalar val, int flag);
    void subEqual(int iSub, GenDistrVector<Scalar> &v, Scalar val, int *ipartial, int flag);

    void subtractKGap( int iSub, GenDistrVector<Scalar> &referenceRHS);
    void subtractMpcRhs(int iSub, GenDistrVector<Scalar> &dv1);
    bool updateActiveSet(GenDistrVector<Scalar> &v, int flag);
    bool updateActiveSet_one(GenDistrVector<Scalar> &v, Scalar tol, int flag);
    int selectOne(GenDistrVector<Scalar> &v, Scalar tol, int flag);
    void subUpdateActiveSet(int iSub, GenDistrVector<Scalar> &v, double tol, int flag, bool &statusChange);

    void projectActiveIneq(GenDistrVector<Scalar> &v);
    void subProjectActiveIneq(int iSub, GenDistrVector<Scalar> &v);
    bool checkInequalities(GenDistrVector<Scalar> &v, int flag = 1);
    void subCheckInequalities(int iSub, GenDistrVector<Scalar> &v, bool &ret, int flag);
    void split(int iSub, GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &v_f, GenDistrVector<Scalar> &v_c, GenDistrVector<Scalar> &v_p);
    void chop(int iSub, GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &v_c, double tol, int chop_flag);
    void getn_u(int iSub, GenDistrVector<Scalar> &n_u, int mpcid);
    void update(Scalar nu, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fp,
                GenDistrVector<Scalar> &ur, GenDistrVector<Scalar> &dur, GenVector<Scalar> &uc, GenVector<Scalar> &duc, int l, double alphabar, double rho = 0.0);
    Scalar alpha_f();
    void expansionStep(GenDistrVector<Scalar> &lambda, Scalar nu, GenDistrVector<Scalar> &p, double alphabar);
    bool linesearch(Scalar &nu, GenDistrVector<Scalar> &p, int &l, double &alphabar);
    bool checkWolfe(Scalar nu, GenDistrVector<Scalar> &p, int which, bool &wolfe1, bool &wolfe2);
    void saveStep();
    void restoreStep();
    void quotient(GenDistrVector<Scalar> &q, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &p);
    void subQuotient(int iSub, GenDistrVector<Scalar> &q, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &p);

  private:
    bool globalFlagCtc;
    double t7;
    int *ngrbmGr;
    int nGroups, nGroups1;
    int *groups;
    Connectivity *groupToSub, *bodyToSub, *subToGroup;
    GenSolver<Scalar> *GtG, *GtGtilda;
    GenSparseMatrix<Scalar> *GtGsparse;
    int *glCrnGroup;
    Connectivity *subToBody;
    int nSubIterDual, nSubIterPrimal, nExtraFp, nRebuildGtG, nRebuildCCt;
    int nLinesearchIter, nStatChDual, nStatChPrimal, nStatChDual_gtg, nStatChPrimal_gtg;
    bool dualStatusChange, primalStatusChange, stepLengthChange;
    int ngrbms;
    int firstAlpha, nRowAlpha;
    Connectivity *coarseConnectGtG;
    SimpleNumberer *eqNumsGtG;
    Connectivity *mpcToMpc;
    CCtSolver<Scalar> *CCtsolver;
    Scalar lambdaTolerance;
    bool mpcPrecon;  // mpc preconditioner flag, true = use generalized preconditioner with CCt
                     // false = use scaling method (diagonal CCt) or no preconditioning
    SimpleNumberer *mpcEqNums;
    Connectivity *mpcToCpu;
    ResizeArray<GenSubDomain<Scalar> *> *subsWithMpcs;
    int numSubsWithMpcs;
    int *mpcSubMap;
    void singularValueDecomposition(FullM &A, FullM &U, int ncol, int nrow, int &rank, double tol, FullM *V = 0);
    double spectralNorm(FullM &A);
    FSCommPattern<int> *mpcPat;
    FSCommPattern<Scalar> *mpcSPat;
    Scalar alpha_l, beta_l; // bisection linesearch parameters
    GenVector<Scalar> e_copy;
    double ee;
    GenDistrVector<Scalar> *lambda_total;
    int newton_iter;

    // Contact functions
  public:
    void makeGtG();
    void deleteGtG() { cerr << "deleteGtG is not implemented\n"; };
    void trMultC(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f);
    void multC(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &cu);
  private:
    void makeE(GenDistrVector<Scalar> &f);
    void assembleE(int iGroup, GenVector<Scalar> &e, GenDistrVector<Scalar> &f);
    void assembleGtG(int iGroup, int flag = 1);
    void rebuildGtGtilda();  
    void computeL0(GenDistrVector<Scalar> &lambda0, GenDistrVector<Scalar> &f);
    void computeL00(int iSub, GenDistrVector<Scalar> &lambda00, GenDistrVector<Scalar> &fr);
    void initRand(GenDistrVector<Scalar> &p);
    void normalizeC();
    void subTrMultC(int iSub, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f);
    void subMultC(int iSub, GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &cu);
    double computeFNorm();
    void project(GenDistrVector<Scalar> &z, GenDistrVector<Scalar> &y, int project_level = 1, bool eflag = false);
    void tProject(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &w, double &dual_error, int project_level = 1, bool pflag = true);
    void multG(GenVector<Scalar> &x, GenDistrVector<Scalar> &y, double alpha, double beta, int flag = 0);
    void subMultG(int iSub, GenVector<Scalar> &x, GenDistrVector<Scalar> &y, double alpha, int flag);
    void trMultG(GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha, double beta, int flag = 0);
    void subTrMultG(int iGroup, GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha, int flag);
    void reSolveGtG(GenVector<Scalar> &x, int flag = 0);
    bool redundant(bool print_flag = false); 
    bool inconsistent(GenVector<Scalar> &e, bool print_warning = true);
    bool feasible(GenDistrVector<Scalar> &lambda, bool check_ieq = true, bool print_warning_eq = true);
    void addRalpha(int iSub, GenDistrVector<Scalar> &u, GenVector<Scalar> &alpha);
    void computeProjectedDisplacement(GenDistrVector<Scalar> &u);
    void addRstar_gT(int iGroup, GenDistrVector<Scalar> &u, GenVector<Scalar> &beta);
    void subtractRstar_g(int iSub, GenDistrVector<Scalar> &u, GenVector<Scalar> &beta);

    // MPC & WI functions
  public:
    void buildCCt();
    void rebuildCCt();
    void deleteCCt() { if(CCtsolver) delete CCtsolver; CCtsolver = 0; }
    Connectivity * getBlockToMpc();
    void cctSolveMpc(GenDistrVector<Scalar> &v);
    void getLocalMpcForces(int iSub, double *mpcLambda);
  private:
    void setBodyRBMoffset(int iSub, int *zColOffset);
    void addMpcRHS(int iMPC, Scalar *singleC);
    void wetInterfaceComms();  // coupled_dph
    void computeLocalWaveNumbers();
  public:
    void rebuildLHSfreq();     // for multiple LHS frequency sweep, assuming no changes except freq/k
                               // ie no change in mpcs, contact, etc.
    void reconstruct();
    void refactor();
    void initNewton() { newton_iter = 0; }
    GenDistrVector<Scalar> * getLambdaTotal() { return lambda_total; }
};


// Class FetiWorkSpace is used to allocate and 
// organize Distributed Vectors for FETI.
template<class Scalar>
class GenFetiWorkSpace  
{
   GenDistrVector<Scalar> *r;
   GenDistrVector<Scalar> *lambda;
   GenDistrVector<Scalar> *w;
   GenDistrVector<Scalar> *y;
   GenDistrVector<Scalar> *z;
   GenDistrVector<Scalar> *p;
   GenDistrVector<Scalar> *pr;
   GenDistrVector<Scalar> *Fp;
   GenDistrVector<Scalar> *Fr;
   GenDistrVector<Scalar> *du;
   GenDistrVector<Scalar> *uzero;
   GenDistrVector<Scalar> *fr;
   GenDistrVector<Scalar> *fr2;
   GenDistrVector<Scalar> *ur;
   GenDistrVector<Scalar> *fw;
   GenDistrVector<Scalar> *deltaL;
   GenDistrVector<Scalar> *q;

   GenDistrVector<Scalar> *wrk1;
   GenDistrVector<Scalar> *wrk2;
   GenDistrVector<Scalar> *deltaU;
   GenDistrVector<Scalar> *deltaF;

   // Extra Distributed Vectors for Nonlinear 
   GenDistrVector<Scalar> *zz; 
   GenDistrVector<Scalar> *rCompare;
   GenDistrVector<Scalar> *uu; 

   // Added Vectors
   GenVector<Scalar> *alpha;
   GenVector<Scalar> *beta;
   GenVector<Scalar> *gamma;
   GenVector<Scalar> *working;
   GenVector<Scalar> *fc;
   GenVector<Scalar> *uc;
   GenVector<Scalar> *duc;
   GenVector<Scalar> *e;

   // copies for contact
   GenDistrVector<Scalar> *lambda_copy;
   GenDistrVector<Scalar> *p_copy;
   GenDistrVector<Scalar> *r_copy;
   GenDistrVector<Scalar> *Fp_copy;
   GenDistrVector<Scalar> *ur_copy;
   GenDistrVector<Scalar> *du_copy;
   GenVector<Scalar> *uc_copy;
   GenVector<Scalar> *duc_copy;

   // initial state for augmented lagrangian
   GenDistrVector<Scalar> *r0;
   GenDistrVector<Scalar> *lambda0;
   GenVector<Scalar> *gamma0;

public:
   GenFetiWorkSpace(DistrInfo& interface, DistrInfo& local, int isNonlinear,
                    int numrbms, int numcrns);
   GenFetiWorkSpace(DistrInfo& interface, DistrInfo& local, DistrInfo& wet, int ngrbms, int numC, bool contact);
   ~GenFetiWorkSpace();

   GenDistrVector<Scalar>& ret_r()       { return *r; }
   GenDistrVector<Scalar>& ret_lambda()  { return *lambda; }
   GenDistrVector<Scalar>& ret_w()       { return *w; }
   GenDistrVector<Scalar>& ret_y()       { return *y; }
   GenDistrVector<Scalar>& ret_z()       { return *z; }
   GenDistrVector<Scalar>& ret_p()       { return *p; }
   GenDistrVector<Scalar>& ret_pr()      { return *pr; }
   GenDistrVector<Scalar>& ret_Fp()      { return *Fp; }
   GenDistrVector<Scalar>& ret_Fr()      { return *Fr; }
   GenDistrVector<Scalar>& ret_du()      { return *du; }
   GenDistrVector<Scalar>& ret_uzero()   { return *uzero; }
   GenDistrVector<Scalar>& ret_wrk1()    { return *wrk1; }
   GenDistrVector<Scalar>& ret_wrk2()    { return *wrk2; }
   GenDistrVector<Scalar>& ret_deltaU()  { return *deltaU; }
   GenDistrVector<Scalar>& ret_deltaF()  { return *deltaF; }
   GenDistrVector<Scalar>& ret_fr()      { return *fr;  }
   GenDistrVector<Scalar>& ret_fr2()     { return *fr2; }
   GenDistrVector<Scalar>& ret_ur()      { return *ur;  }
   GenDistrVector<Scalar>& ret_fw()      { return *fw;  }
   GenDistrVector<Scalar>& ret_r0()      { return *r0;  }
   GenDistrVector<Scalar>& ret_lambda0() { return *lambda0; }
   GenDistrVector<Scalar>& ret_deltaL()  { return *deltaL; }
   GenDistrVector<Scalar>& ret_q()       { return *q; }
   GenDistrVector<Scalar>& ret_zz()      { return *zz; }
   GenDistrVector<Scalar>& ret_rCompare(){ return *rCompare; }
   GenDistrVector<Scalar>& ret_uu()      { return *uu; }
   GenDistrVector<Scalar>& ret_lambda_copy() { return *lambda_copy; }
   GenDistrVector<Scalar>& ret_p_copy()  { return *p_copy; }
   GenDistrVector<Scalar>& ret_r_copy()  { return *r_copy; }
   GenVector<Scalar>& ret_alpha()        { return *alpha; }
   GenVector<Scalar>& ret_beta()         { return *beta;  }
   GenVector<Scalar>& ret_gamma()        { return *gamma; }
   GenVector<Scalar>& ret_working()      { return *working; }
   GenVector<Scalar>& ret_fc()           { return *fc;  }
   GenVector<Scalar>& ret_uc()           { return *uc; }
   GenVector<Scalar>& ret_duc()          { return *duc; }
   GenVector<Scalar>& ret_e()            { return *e; }
   GenVector<Scalar>& ret_gamma0()       { return *gamma0; }

   // extra functions for contact 
   void save();
   void restore(bool flag = true);
   void save_initial();
   void save_lambda();
   void save_p();

   void clean_up();
   void zeroPointers();
};

typedef GenFetiSolver<double> FetiSolver;
typedef GenFetiDPSolver<double> FetiDPSolver;
typedef GenFetiWorkSpace<double> FetiWorkSpace;

// HB: MPC stuff -> for creating "superblocks"
struct BlockPair {

  int Id;
  double cost;
  double bandwidth;

  BlockPair() { Id = 0; cost = 0.0; bandwidth = 0.0; }
  BlockPair(int _Id, double _cost, double _bandwidth)
    { Id = _Id; cost = _cost; bandwidth = _bandwidth; }
  BlockPair(int _Id, double _cost)
    { Id = _Id; cost = _cost; bandwidth = 0.0; }
  BlockPair(int _Id)
    { Id = _Id; cost = 0.0; bandwidth = 0.0; }
  BlockPair(double _cost)
    { Id = 0;  cost = _cost; bandwidth = 0.0; }
  bool operator <(const BlockPair &bp) const {
    return (cost < bp.cost) || ( (cost == bp.cost) && (Id < bp.Id) ) ;
   }
  void print() {
    filePrint(stderr," --- block ---\n");
    filePrint(stderr,"  # Id       : %d\n",Id);
    filePrint(stderr,"  # cost     : %e\n",cost);
    filePrint(stderr,"  # bandwidth: %e\n",bandwidth);
  }
};

#ifdef _TEMPLATE_FIX_
  #include <Feti.d/Feti.C>
  #include <Feti.d/FetiDP.C>
  #include <Feti.d/NLFeti.C>
#ifdef DISTRIBUTED
  #include <Dist.d/DistFeti.C>
#endif
#endif

#endif


