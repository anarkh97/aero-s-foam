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
#include <iostream>

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

template <typename Scalar>
class FetiSub;

// The FetiSolver class is too big. Future versions need it to be leaner
template <class Scalar>
class GenFetiSolver  : public GenParallelSolver<Scalar>
{
protected:
	std::vector<FetiSub<Scalar> *> subdomains;
	GenSubDomain<Scalar> **sd;
	int nsub, glNumSub;
	FetiInfo *fetiInfo;
	int verboseFlag;
	FSCommunicator *fetiCom;
	int myCPU, numCPUs;
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
	mutable GenCGOrthoSet<Scalar> *oSetCG; // Workspace
	mutable GenGMRESOrthoSet<Scalar> *oSetGMRES;
	mutable GenGCROrthoSet<Scalar> *oSetGCR;
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
	mutable GenFetiOpControler<Scalar> *opControl;
	GenSolver<Scalar> *GtGsolver;
	GenBigMatrix<Scalar> *PCtFPC;
	GenFetiWorkSpace<Scalar> *wksp;
	int QGisLocal;  // Whether or not QG is local
	int isDynamic;  // Whether or not we are in dynamics
	int isFeti2;    // Whether or not we are using FETI2
	mutable int numSystems; // Nonlinear additions
	mutable Timings times;
	GenSymFullMatrix<Scalar> *GtQGs;
	GenFullM<Scalar> *GtFCs;
	Scalar *gtqglocal;
	int numNodes;
	Connectivity *cpuToSub;
	int glNumMpc, glNumMpc_primal;

	void makeGtG();
	void makeDistGtG(int *glSubToLocal);
	void assembleDistGtQGs(int i, int *);
	bool isLowestLocalNeighbor(int subI, int subJ);
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
	void fSend(int i,  GenDistrVector<Scalar> &) const;
	void fScale(int i, GenDistrVector<Scalar> &) const;
	void fSplit(int iSub, GenDistrVector<Scalar> &force);
	void fSendCoupled(int iSub, GenDistrVector<Scalar> &force, GenDistrVector<Scalar> &fw) const;
	void fScaleCoupled(int iSub, GenDistrVector<Scalar> &force, GenDistrVector<Scalar> &fw) const;

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
	void singleCoarseSolve(const GenDistrVector<Scalar> &rhs, GenDistrVector<Scalar> &u) const;
	void setAllOffsets(int iSub, int gOffset);
	void getSRMult(int iSub, GenDistrVector<Scalar> *r, GenDistrVector<Scalar> *lambda,
	               Scalar *alpha);
	void getSGtMult(int iSub, GenDistrVector<Scalar> *r, Scalar *alpha);
	void getSCtMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv);
	void getSQtMult(int iSub, GenDistrVector<Scalar> *r, Scalar *sv);
	void getFGMult(int iSub, GenDistrVector<Scalar> *r, Scalar *alpha);
	void addG(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const;
	void addSG(int iSub, GenDistrVector<Scalar> *r, Scalar *sv) const;
	void addGs(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &w, GenVector<Scalar> &) const;
	void singlePr(GenDistrVector<Scalar> &y, GenDistrVector<Scalar> &p, GenVector<Scalar> &) const;
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
	//GenFetiSolver() { initialize(); };
	GenFetiSolver(int nsub, GenSubDomain<Scalar> **, Connectivity *,
	              FetiInfo *finfo, FSCommunicator *fetiCom, int *glToLoc,
	              Connectivity *mpcToSub, Connectivity *cpuToSub,
	              GenSolver<Scalar> **sysMatrices = 0, GenSparseMatrix<Scalar> **sysMat = 0,
	              Rbm **_rbms = 0, int verboseFlag = 0);
	GenFetiSolver(int _nsub, GenSubDomain<Scalar> **subs, int _numThreads, int _verboseFlag);
	virtual ~GenFetiSolver();

	void sendDeltaF(int iSub, GenDistrVector<Scalar>& deltaF);
	void normDeltaF(int iSub, double * subDots, GenDistrVector<Scalar>* deltaF);
	void factorMatrices(int isub);
	void sendScale(int isub);
	void collectScale(int isub);
	void getErrorEstimator(int iSub, GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &es);
	void interfaceDiff(int iSub, GenDistrVector<Scalar> &v) const;
	void rebuildInterface(int iSub, GenDistrVector<Scalar> &v);
	void multKbb(int iSub, const GenDistrVector<Scalar>& v, GenDistrVector<Scalar> &interfvec,
	             GenDistrVector<Scalar>& deltaU, GenDistrVector<Scalar> &deltaF, bool &errorFlag) const;
	void localSolve(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
	                GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4) const;
	void localSolve2(int iSub, GenDistrVector<Scalar> *v1, GenDistrVector<Scalar> *v2,
	                 GenVector<Scalar> *beta, GenDistrVector<Scalar> *v3,
	                 GenDistrVector<Scalar> *v4) const;
	void interfSend(int iSub, GenDistrVector<Scalar> &dv1);
	//void interfDiff(int iSub, GenDistrVector<Scalar> &dv1);
	void interfDiffAndDot(int iSub, GenDistrVector<Scalar> &dv1, GenDistrVector<Scalar> &dv2) const;
	void getRMult(int iSub, GenDistrVector<Scalar> *localvec, Scalar *alpha);
	void getGtMult(int iSub, GenDistrVector<Scalar> *localvec, Scalar *alpha);
	void addRP(int iSub, GenDistrVector<Scalar> * localvec, Scalar *alpha);
	void addRS(int iSub, GenDistrVector<Scalar> * localvec, Scalar *alpha);
	void solve(const GenDistrVector<Scalar> &rhs, GenDistrVector<Scalar> &x) override;
	void distributeForce(GenDistrVector<Scalar> &force) const;
	void distributeForce(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fw) const;
	void reSolve(GenDistrVector<Scalar> &) override;
	Timings& getTimers() { return times; }
	void makeRbmPat();
	void makeSingleIntPat();

	// NonLinear Functions
	virtual void reBuild(FullSquareMatrix **kel, DistrGeomState& gs, int iter=0,
	                     int step = 1);
	void reBuildMatrices(int isub, FullSquareMatrix **kel);
	void reBuildErrorEstimator(GenFullSquareMatrix<Scalar> **kel);
	void subdomainReBuild(int isub, FullSquareMatrix **kel, DistrGeomState *gs);
	int nlPreCondition(const GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Pr) const;
	double getSolutionTime() override { return times.solve; }
	void reSendInterfaceRBM(int iSub);
	void reGetNeighbQGs(int iSub);
	void reMakeLocalFcoarse();
	void reGetNeighbFC(int iSub);
	void reComputeFiBC(int iSub);
	void makeCompatible(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &deltaU) const;
	Scalar localSolveAndJump(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &,
	                         GenDistrVector<Scalar> &, GenDistrVector<Scalar> &) const;
	Scalar localSolveAndJump(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &,
	                         GenVector<Scalar> &, GenDistrVector<Scalar> &,
	                         GenDistrVector<Scalar> &) const;
	void project(GenDistrVector<Scalar> &, GenVector<Scalar> &, GenDistrVector<Scalar> &, int isDirect = 1) const;
	void tProject(GenDistrVector<Scalar> &, GenVector<Scalar> &, GenDistrVector<Scalar> &, int isDirect = 1) const;
	void computeDynamL0(GenDistrVector<Scalar> &, GenVector<Scalar> &,
	                    GenDistrVector<Scalar> &, GenDistrVector<Scalar> &) const;
	void computeL0(GenDistrVector<Scalar> &, GenVector<Scalar> &, GenDistrVector<Scalar> &) const;
	void computeL0(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &, GenDistrVector<Scalar> &,
	               GenVector<Scalar> &, GenDistrVector<Scalar> &) const;
	void addR(GenDistrVector<Scalar> &, GenVector<Scalar> &) const;
	void addRSingle(GenDistrVector<Scalar> &, GenVector<Scalar> &) const;
	void outputPrimalResidual(int iter, GenDistrVector<Scalar> &deltaF) const;

	// preconditioner returns the error estimator
	double preCondition(const GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Pv, bool errorFlag = true) const;

	void resetOrthoSet();
	void orthoAddCG(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar pFp) const;
	void orthogonalize(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &) const;

	// GMRES functions
	void initGMRES(GenDistrVector<Scalar> &p);
	double orthoAddGMRES(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp);
	void GMRESSolution(GenDistrVector<Scalar> &p);

	// GCR functions
	int predictGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &lambda0);
	void orthogonalizeGCR(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fr,
	                      GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp);
	void orthoAddGCR(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &Fp, Scalar FpFp);

	/** \brief Create a guess for lambda0.
	 *
	 * @param r RHS for which a guess is desired.
	 * @param lambda0 The best guess lambda.
	 * @return Whether a prediction was made.
	 */
	bool predict(const GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &lambda0) const;

	DistrInfo &interfInfo() { return interface; }
	DistrInfo &localInfo()  { return internalDI; }
	GenDistrVector<Scalar> & getLambda() { return wksp->ret_lambda(); }
	int neq() { return internalDI.len; }

	void updateFeti2lambda(GenDistrVector<Scalar> &, GenDistrVector<Scalar>&,
	                       GenDistrVector<Scalar> &,
	                       GenVector<Scalar>&, GenVector<Scalar>&) const;
	void updateFeti2y(GenDistrVector<Scalar> &, GenDistrVector<Scalar>&,
	                  GenVector<Scalar>&, GenVector<Scalar>&) const;
	void computeRgcTransMult(GenVector<Scalar>&, GenVector<Scalar>&) const;
	void computeRgcMult(GenVector<Scalar>& beta, GenVector<Scalar>& result) const;
	void updateFeti2Vector(GenDistrVector<Scalar> &, GenVector<Scalar>&, GenVector<Scalar>&) const;

	void addAllFcoarse(GenFullM<Scalar> &);
	void finishRgc(int, int);

	void getCtMult(GenDistrVector<Scalar> &w, GenVector<Scalar> &gamma) const;

	void scatterHalfInterface(int iSub, Scalar *v1, GenDistrVector<Scalar> *v2) const;
	void gatherHalfInterface(int iSub, const GenDistrVector<Scalar> *v1,
	                         const GenDistrVector<Scalar> *v2,
	                         Scalar *v3, Scalar *v4) const;

	void setAndStoreInfo(int iter, double finalPrimal2, double finalDual2 ) const;
	void findProfileSizes(int iSub, int *subSizes);

	// For eigen problem
	virtual int numRBM() override;
	virtual void getRBMs(GenDistrVectorSet<Scalar> &);
	virtual void getRBMs(Scalar *);
	void Ksolve(int iSub, GenStackDistVector<Scalar> &R);

	int halfOffset(int iSub) const { return fetiOps[iSub]->halfOffset; }
	int numNeighbor(int iSub) const;
	Scalar *interfaceBuffer(int iSub) { return fetiOps[iSub]->interfBuff; }
	virtual void clean_up();
	double getFNormSq(GenDistrVector<Scalar> &f) override;

	virtual void getLocalMpcForces(int iSub, double *mpcLambda) { };  // only implemented for DP

protected:
	FSCommPattern<Scalar> *wiPat;
};

template<class Scalar>
class GenFetiDPSolver : public GenFetiSolver<Scalar>
{
	using GenFetiSolver<Scalar>::fetiInfo;
	using GenFetiSolver<Scalar>::verboseFlag;

	DistrInfo internalR, internalC, internalWI;
	DistrInfo *coarseInfo;
	GenSolver<Scalar>         *KccSolver;
	GenParallelSolver<Scalar> *KccParallelSolver;
	Scalar *kccrbms;
	GenSparseMatrix<Scalar> *KccSparse;
	int glNumCorners;
	Connectivity *cornerToSub;
	DofSetArray *cornerEqs;
	int mpcOffset, augOffset; // mpc equation offset for coarse grid
	void initialize();
	bool rbmFlag;
	bool geometricRbms;
	enum StepType { CG, PROPORTIONING, EXPANSION };
	bool proportional;

public:
	GenFetiDPSolver(int nsub, int glNumSub, GenSubDomain<Scalar> **sd, Connectivity *subToSub,
	                FetiInfo *finfo, FSCommunicator *fetiCom, int *glToLoc, Connectivity *mpcToSub, Connectivity *mpcToSub_primal,
	                Connectivity *mpcToMpc, Connectivity *mpcToCpu, Connectivity *cpuToSub,
	                Connectivity *bodyToSub = 0, GenSolver<Scalar> **sysMatrices = 0,
	                GenSparseMatrix<Scalar> **sysMat = 0, Rbm **rbms = 0, bool rbmFlag = 0,
	                bool geometricRbms = true, int verboseFlag = 0);
	virtual ~GenFetiDPSolver();

	//int ngrbm;
	void makeKbb(int iSub);
	void makeFc(int iSub, GenDistrVector<Scalar> &fr, GenDistrVector<Scalar> &lambda);
	void makeFcB(int iSub, GenDistrVector<Scalar> &bf);
	void KrrReSolve(int iSub, GenDistrVector<Scalar> &ur);
	void makeKcc();
	void deleteKcc() { /* not implemented */ };
	void setSysMatrices(GenSolver<Scalar> **sysMatrices, GenSparseMatrix<Scalar> **sysMat) { /* not implemented */ };
	void assembleFcStar(GenVector<Scalar> &FcStar) const;
	void mergeSolution(GenDistrVector<Scalar> &ur, GenVector<Scalar> &uc, GenDistrVector<Scalar> &u,
	                   GenDistrVector<Scalar> &lambda) const;
	void mergeUr(int iSub, GenDistrVector<Scalar> &ur, GenVector<Scalar> &uc, GenDistrVector<Scalar> &u,
	             GenDistrVector<Scalar> &lambda) const;
	/** \brief Compute \f$ f_r = K_{rc} u_c \f$ for subdomain iSub */
	void multKrc(int iSub, GenDistrVector<Scalar> &fr, const GenVector<Scalar> &uc) const;
	void extractFr(int iSub, const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fr) const;
	void extractFc(int iSub, const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fc) const;
	void extractFw(int iSub, const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fw) const;
	void getFc(const GenDistrVector<Scalar> &f, GenVector<Scalar> &fc) const;
	void makeEdgeConnectivity();
	void countEdges(int iSub, int *edges) const;
	void numberEdges(int iSub, int *eP, int *ep2, int *edges, FSCommPattern<int> *sPat);
	void receiveNeighbEdgeNums(int iSub, int *eP, int *edges, FSCommPattern<int> *sPat);
	void factorLocalMatrices(int isub);
	/** \brief Extract vector, modifies f.
	 *
	 * @param f
	 * @param fr
	 * @param fc
	 * @param fw
	 * @return
	 */
	double extractForceVectors(GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &fr,
	                           GenVector<Scalar> &fc, GenDistrVector<Scalar> &fw) const;
	void printSummary(int iter) const;
	double getFNormSq(GenDistrVector<Scalar> &f);
	void getRBMs(Scalar *);
	void getRBMs(GenDistrVectorSet<Scalar> &);
	void getGlobalRBM(int iSub, int &iRBM, GenDistrVector<Scalar> &R);
	int numRBM();
	void clean_up();
	void subdomainSolve(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
	                    GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4,
	                    GenVector<Scalar> &v5) const;
	void subdomainSolveCoupled1(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
	                            GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4,
	                            GenVector<Scalar> &v5) const;
	void subdomainSolveCoupled2(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
	                            GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4,
	                            GenVector<Scalar> &v5, GenDistrVector<Scalar> &fw) const;
	void subdomainSolveCoupled2b(int iSub, GenDistrVector<Scalar> &v1, GenDistrVector<Scalar> &v2,
	                            GenDistrVector<Scalar> &v3, GenDistrVector<Scalar> &v4,
	                            GenVector<Scalar> &v5) const;
	void localSolveAndJump(GenDistrVector<Scalar> &fr, GenDistrVector<Scalar> &lambda,
	                       GenDistrVector<Scalar> &ur, GenVector<Scalar> &fc,
	                       GenVector<Scalar> &uc, GenDistrVector<Scalar> &r,
	                       GenDistrVector<Scalar> &fw) const;
	Scalar localSolveAndJump(GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &dur,
	                         GenVector<Scalar> &duc, GenDistrVector<Scalar> &Fp) const;
	void solve(const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
	void solveCG(const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
	void solveGMRES(const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
	void solveGCR(const GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &u);
	double preCondition(const GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &Pv, bool errorFlag = true) const;

	void subtractMpcRhs(int iSub, GenDistrVector<Scalar> &dv1);
	bool updateActiveSet(GenDistrVector<Scalar> &v, int flag, double tol = 0.0);
	void subUpdateActiveSet(int iSub, GenDistrVector<Scalar> &v, double tol, int flag, bool *statusChange);
	void subRecvMpcStatus(int iSub, FSCommPattern<int> *mpcPat, int flag, bool *statusChange);

	void projectActiveIneq(const GenDistrVector<Scalar> &x, GenDistrVector<Scalar> &y) const;
	void subProjectActiveIneq(int iSub, GenDistrVector<Scalar> &v) const;
	void split(int iSub, GenDistrVector<Scalar> &v, GenDistrVector<Scalar> &v_f, GenDistrVector<Scalar> &v_c);
	void update(Scalar nu, GenDistrVector<Scalar> &lambda,
	            GenDistrVector<Scalar> &p, GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &Fp,
	            GenDistrVector<Scalar> &ur, GenDistrVector<Scalar> &dur,
	            GenVector<Scalar> &uc, GenVector<Scalar> &duc) const;
	void saveStep() const;
	void restoreStep() const;

private:
	bool globalFlagCtc;
	mutable double t7;
	int *ngrbmGr;
	int nGroups, nGroups1;
	int *groups;
	Connectivity *groupToSub, *bodyToSub, *subToGroup;
	GenSolver<Scalar> *GtGtilda;
	GenSparseMatrix<Scalar> *GtGsparse;
	Connectivity *subToBody;
	/// Statistic variables.
	mutable int nSubIterDual, nSubIterPrimal, nMatVecProd, nRebuildGtG, nRebuildCCt;
	mutable int nLinesearch, nLinesearchIter, nStatChDual, nStatChPrimal;
	mutable bool dualStatusChange, primalStatusChange, stepLengthChange;
	int ngrbms;
	Connectivity *coarseConnectGtG;
	SimpleNumberer *eqNumsGtG;
	Connectivity *mpcToMpc;
	CCtSolver<Scalar> *CCtsolver;
	bool mpcPrecon;  // mpc preconditioner flag, true = use generalized preconditioner with CCt
	// false = use scaling method (diagonal CCt) or no preconditioning
	SimpleNumberer *mpcEqNums;
	Connectivity *mpcToCpu;
	ResizeArray<GenSubDomain<Scalar> *> *subsWithMpcs;
	int numSubsWithMpcs;
	int *mpcSubMap;
	void singularValueDecomposition(FullM &A, FullM &U, int ncol, int nrow, int &rank, double tol, FullM *V = 0);
	FSCommPattern<int> *mpcPat;
	mutable double lastError;

	// Contact functions
public:
	void makeGtG();
	void deleteGtG() { std::cerr << "deleteGtG is not implemented\n"; };
	void trMultC(GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f);
	void multC(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &cu);
private:
	void makeE(GenDistrVector<Scalar> &f) const; // Modifies the workspace.
	void assembleE(int iGroup, GenVector<Scalar> &e, GenDistrVector<Scalar> &f);
	void assembleGtG(int iGroup);
	void rebuildGtGtilda();
	void computeL0(GenDistrVector<Scalar> &lambda0, GenDistrVector<Scalar> &f) const;
	void normalizeC();
	void subTrMultC(int iSub, GenDistrVector<Scalar> &lambda, GenDistrVector<Scalar> &f);
	void subMultC(int iSub, GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &cu);
	double computeFNorm();
	void project(GenDistrVector<Scalar> &z, GenDistrVector<Scalar> &y, int eflag = 0) const;
	double tProject(GenDistrVector<Scalar> &r, GenDistrVector<Scalar> &w) const;
	void multG(const GenVector<Scalar> &x, GenDistrVector<Scalar> &y, double alpha, double beta) const;
	void subMultG(int iSub, const GenVector<Scalar> &x, GenDistrVector<Scalar> &y, double alpha) const;
	void trMultG(const GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha, double beta) const;
	void subTrMultG(int iGroup, const GenDistrVector<Scalar> &x, GenVector<Scalar> &y, double alpha) const;
	void addRalpha(int iSub, GenDistrVector<Scalar> &u, GenVector<Scalar> &alpha) const;
	void computeProjectedDisplacement(GenDistrVector<Scalar> &u) const;
	void addRstar_gT(int iGroup, GenDistrVector<Scalar> &u, GenVector<Scalar> &beta) const;
	void subtractRstar_g(int iSub, GenDistrVector<Scalar> &u, GenVector<Scalar> &beta) const;
	bool checkStoppingCriteria(int iter, double error, double fnorm) const;

	// MPC & WI functions
public:
	void buildCCt();
	void rebuildCCt();
	void deleteCCt() { if(CCtsolver) delete CCtsolver; CCtsolver = 0; }
	Connectivity * getBlockToMpc();
	void cctSolveMpc(GenDistrVector<Scalar> &v) const;
	void getLocalMpcForces(int iSub, double *mpcLambda);
private:
	void setBodyRBMoffset(int iSub, int *zColOffset);
	void addMpcRHS(int iMPC, Scalar *singleC);
	void wetInterfaceComms();  // coupled_dph
	void computeLocalWaveNumbers();
public:
	void reconstruct();
	void refactor();
	void reconstructMPCs(Connectivity *_mpcToSub, Connectivity *_mpcToMpc, Connectivity *_mpcToCpu);
	void zeroG();
	void deleteG();
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
	GenDistrVector<Scalar> *du;
	GenDistrVector<Scalar> *uzero;
	GenDistrVector<Scalar> *fr;
	GenDistrVector<Scalar> *fr2;
	GenDistrVector<Scalar> *ur;
	GenDistrVector<Scalar> *fw;
	GenDistrVector<Scalar> *wrk1;
	GenDistrVector<Scalar> *wrk2;
	GenDistrVector<Scalar> *deltaU;
	GenDistrVector<Scalar> *deltaF;

	// Extra Distributed Vectors for Nonlinear (not used for feti-dp)
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

	// Contact
	GenDistrVector<Scalar> *lambda_copy;
	GenDistrVector<Scalar> *p_copy;
	GenDistrVector<Scalar> *r_copy;
	GenDistrVector<Scalar> *Fp_copy;
	GenDistrVector<Scalar> *du_copy;
	GenVector<Scalar> *uc_copy;
	GenVector<Scalar> *duc_copy;
	GenDistrVector<Scalar> *gc;
	GenDistrVector<Scalar> *gf;

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
	//GenDistrVector<Scalar>& ret_Fr()      { return *Fr; }
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
	GenDistrVector<Scalar>& ret_zz()      { return *zz; }
	GenDistrVector<Scalar>& ret_rCompare(){ return *rCompare; }
	GenDistrVector<Scalar>& ret_uu()      { return *uu; }
	GenDistrVector<Scalar>& ret_lambda_copy() { return *lambda_copy; }
	GenDistrVector<Scalar>& ret_p_copy()  { return *p_copy; }
	GenDistrVector<Scalar>& ret_r_copy()  { return *r_copy; }
	GenDistrVector<Scalar>& ret_gc()      { return *gc; }
	GenDistrVector<Scalar>& ret_gf()      { return *gf; }
	GenVector<Scalar>& ret_alpha()        { return *alpha; }
	GenVector<Scalar>& ret_beta()         { return *beta;  }
	GenVector<Scalar>& ret_gamma()        { return *gamma; }
	GenVector<Scalar>& ret_working()      { return *working; }
	GenVector<Scalar>& ret_fc()           { return *fc;  }
	GenVector<Scalar>& ret_uc()           { return *uc; }
	GenVector<Scalar>& ret_duc()          { return *duc; }
	GenVector<Scalar>& ret_e()            { return *e; }
	// extra functions for contact
	void save();
	void restore();
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
#include <Feti.d/FetiDPCore.C>
#endif

#endif
