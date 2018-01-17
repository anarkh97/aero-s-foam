//
// Created by Michel Lesoinne on 11/3/17.
//

#ifndef FEM_FETUSUB_H
#define FEM_FETUSUB_H
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Math.d/matrix.h>
#include <Utils.d/GlobalToLocalMap.h>
#include <Driver.d/SComm.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/MpcSparse.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/FsiSparse.h>

class FSCommStructure;
template <typename Scalar>
class FSCommPattern;
template <typename Scalar>
class GenVector;
template <typename Scalar>
class GenFullM;
template <typename Scalar>
class GenSolver;
template <typename Scalar>
class GenSparseMatrix;
template <typename Scalar>
class GenAssembledFullM;
template <typename Scalar>
class GenCuCSparse;
template <typename Scalar>
class SubLMPCons;
template <typename Scalar>
class GenSparseSet;

class Connectivity;
class CoordSet;
class DofSet;
class FetiInfo;


/** \brief Pure Interface of what a the notion of Subdomain provides for FETI solver. */
class FetiBaseSub {
public:
	/// \brief Obtain the solver settings. TODO Get rid of this. Why should the subdomain data know all the solver details?
	virtual const FetiInfo &getFetiInfo() const = 0;
	/** \brief Obtain the size of the interface of this subdomain. */
	virtual int interfLen() const = 0;
	/** \brief Obtain the size of the half interface for which this subdomain is the master. */
	virtual int halfInterfLen() const = 0;
	/** \brief Obtain the index of this subdomain in the global system. */
	virtual int subNum() const = 0;
	/** \brief Set the communication size in the pattern.
	 * TODO abstract FSCommPattern to isolate the scalar type.
	 *
	 * @param numRBM The number of RBMs for which the setup should be done.
	 * @param pattern The pattern to adjust.
	 */
	virtual void setRbmCommSize(int numRBM, FSCommStructure *pattern) const = 0;
	/** \brief Set the communication pattern from this subdomain to its neighbors for a RHS or solution vector. */
	virtual void setDofCommSize(FSCommStructure *) const = 0;

	virtual void setCommSize(FSCommStructure *pat, int size) const = 0;

	virtual void setMpcCommSize(FSCommStructure *mpcPat) const = 0;

	virtual void setMpcNeighbCommSize(FSCommPattern<int> *pt, int size) const = 0;

	virtual void setGCommSize(FSCommStructure *pat) const = 0;

	virtual int localSubNum() const = 0;
	virtual int localLen() const = 0;
	virtual int localRLen() const = 0;

	virtual const int *getGlNodes() const = 0;

	virtual int getNumUncon() const = 0;

	/** \brief Make the subdomain determine its master flag, i.e. whether it is the primary holder of a variable. */
	virtual void computeMasterFlag(const Connectivity &mpcToSub)  = 0;
	virtual const bool* getMasterFlag() const = 0;

	/// \brief Obtain the number of MPC constraints.
	virtual int numMPCs() const = 0;

	virtual int numMPCs_primal() const = 0;
	/// \brief Obtain the number of coarse dofs. This method computes a cached value.
	virtual int numCoarseDofs() = 0;
	/// \brief Obtain the number of corner nodes.
	virtual int numCorners() const = 0;

	const std::vector<int> &getLocalCornerNodes() const { return cornerNodes; };

	virtual int numWetInterfaceDofs() const = 0;

	virtual int getLocalMPCIndex(int globalMpcIndex) const = 0;

	virtual int getGlobalMPCIndex(int localMpcIndex) const = 0;

	virtual const CoordSet& getNodeSet() const = 0;

	virtual bool isEdgeNeighbor(int neighb) const = 0;

	virtual double getShiftVal() const = 0;

	const std::vector<int> &getCornerNodes() const { return glCornerNodes; }
	std::vector<int> &getCornerNodes() { return glCornerNodes; }
	void markCornerDofs(int *glCornerDofs) const;
	void makeKccDofs(DofSetArray *cornerEqs, int augOffset, Connectivity *subToEdge, int mpcOffset = 0);

	int numEdgeDofs(int i) const { return edgeDofSize[i]; }

	const std::vector<int> &getWeights() const { return weight; }
	std::vector<int> &getWeights() { return weight; }
	int dofWeight(int i) const { return weight[i]; }
	/* missing:
	 * splitInterf
	 * setMpcNeighbCommSize
	 * setGCommSize
	 * assembleTrbmE
	 * assembleE
	 */
	std::vector<DofSet> cornerDofs;
	DofSet *edgeDofs;      // JAT 112113

	SComm *scomm = nullptr;

	const auto &getCornerEqNums() const { return cornerEqNums; }
	int getGroup() const { return group; }

	SComm *getSComm() { return scomm; }
	void setSComm(SComm *sc);
	void setSendData(int neighb, void *data)
	{ scomm->setExchangeData(neighb, data); }
	void *getExchangeData(int neighbID) // Communication mechanism
	{ return scomm->getExchangeData(neighbID); }
	void *getExchangePointer(int neighbID) // Communication mechanism
	{ return scomm->exchangeData[neighbID]; }
	/** \brief Obtain the number of neighbor subdomains. */
	int numNeighbors() const { return scomm->numNeighb;}
	int numEdgeNeighbors() const { return scomm->numEdgeNeighb; }

	bool* getMpcMaster() const { return mpcMaster; }
	// Multiple Point Constraint (MPC) functions
	int getNumMpc() const       { return numMPC; }

	virtual ConstrainedDSA *get_c_dsa() const = 0;
	ConstrainedDSA *getCCDSA() const;

	void addMPCsToGlobalZstar(FullM *globalZstar, int startRow, int startCol, int numCol);
	void addSPCsToGlobalZstar(FullM *globalZstar, int &zRow, int zColOffset);

	void setWIoneCommSize(FSCommStructure *pat) const;
	void setWICommSize(FSCommStructure *wiPat);
	void setWImapCommSize(FSCommPattern<int> *pat);
	bool isWetInterfaceDof(int d) const { return (wetInterfaceMap.size() > 0) ? (wetInterfaceMap[d] > -1) : false; }
	void GramSchmidt(double *Q, bool *isUsed, DofSet desired, int nQPerNeighb, bool isPrimalAugmentation);

protected:
	bool isCoupled = false; // TODO Ensure this is set or derived from some other info.
	int boundLen = 0;
	int internalLen = 0;
	int totalInterfSize;
	std::vector<int> allBoundDofs;
	std::vector<int> boundMap;
	std::vector<int> internalMap;
	/// \brief Corner nodes in local numbering.
	std::vector<int> cornerNodes;
	std::vector<bool> isCornerNode;   //<! \brief True for node which is a corner node; false otherwise.
	std::vector<int> glCornerNodes; //!< \brief Corner nodes in global numbering.
	int numCRN = 0;
	int numCRNdof = 0;
	std::vector<std::vector<DofSet>> boundaryDOFs;
	std::vector<int> edgeDofSizeTmp;   // XXXX
	std::vector<int> edgeDofSize;      //<! \brief Number of edge DOF per neighbor.
	std::vector<int> cornerEqNums; //<! \brief unique equation numbers for subdomain corner dofs
	std::unique_ptr<ConstrainedDSA> cc_dsa;
	std::vector<int> ccToC; //!< Mapping from cc_dsa to c_dsa. All indices are >= 0
	std::vector<int> cToCC; //!< Mapping from c_dsa to cc_dsa. Indices for corner DOFs are < 0.
	bool isMixedSub = false;

	std::vector<int> weight; ///!< \brief DOF weights (i.e. number of subd sharing that dof).
	std::vector<int> weightPlus; ///!< \brief DOF weights (i.e. number of subd sharing that dof) including corner DOFs.

	// MPC related data
	Connectivity *localMpcToMpc = nullptr;
	Connectivity *localMpcToGlobalMpc = nullptr;
	bool *faceIsSafe = nullptr;
	int *localToGroupMPC = nullptr;
	int *boundDofFlag = nullptr;  // boundDofFlag[i] = 0 -> perfect interface dof  (not contact or mpc)
	// boundDofFlag[i] = 1 -> node-to-node contact interface dof
	// boundDofFlag[i] = 2 -> mpc dof, only used for rixen method, domain->mpcflag = 1
	// note: boundDofFlag[i] > 2 can be used for future extensions, eg mortar contact
	bool *masterFlag = nullptr; // masterFlag[i] = true if this sub is the "master" of allBoundDofs[i]
	bool *internalMasterFlag = nullptr;
	int masterFlagCount = 0;
	bool *mpcMaster = nullptr;  // mpcMaster[i] = true if this subd contains masterdof for mpc i
	Connectivity *mpcToDof = nullptr;
	Connectivity *localMpcToBlock = nullptr;
	Connectivity *blockToLocalMpc = nullptr;
	Connectivity *blockToBlockMpc = nullptr;
	Connectivity *localMpcToBlockMpc = nullptr;
	Connectivity *mpcToBoundDof = nullptr;
	double *localLambda = nullptr;  // used for contact pressure output
	std::vector<int> invBoundMap;
	int *mpclast = nullptr;

	mutable int *mpcStatus;
	mutable bool *mpcStatus1, *mpcStatus2;
public:
	void setGroup(Connectivity *subToGroup) { this->group = (*subToGroup)[subNum()][0]; }
	void setNumGroupRBM(int *ngrbmGr);
	void getNumGroupRBM(int *ngrbmGr);
	void makeLocalToGroupMPC(Connectivity *groupToMPC);

	GlobalToLocalMap &getGlobalToLocalNode() { return glToLocalNode; }

	int group = 0;
	// Multiple Point Constraint (MPC) Data
	int numMPC = 0;             // number of local Multi-Point Constraints
	int *localToGlobalMPC = nullptr;  // local to global MPC numbering
	GlobalToLocalMap globalToLocalMPC; // alternative data structure for global to local MPC numbering
	// not a pointer so don't have to de-reference before using [] operator

	int numMPC_primal = 0;
	int *localToGlobalMPC_primal = nullptr;
	GlobalToLocalMap globalToLocalMPC_primal;

	int *cornerMap = nullptr;

	void sendNumNeighbGrbm(FSCommPattern<int> *pat);
	void recvNumNeighbGrbm(FSCommPattern<int> *pat);

	void sendNumWIdof(FSCommPattern<int> *sPat) const;
	void recvNumWIdof(FSCommPattern<int> *sPat);
	void sendWImap(FSCommPattern<int> *pat);
	void recvWImap(FSCommPattern<int> *pat);

	void sendNeighbGrbmInfo(FSCommPattern<int> *pat);
	void receiveNeighbGrbmInfo(FSCommPattern<int> *pat);

	bool isWetInterfaceNode(int n) const { return (wetInterfaceNodeMap.size() > 0) ? (wetInterfaceNodeMap[n] > -1) : false; }

	// variables and routines for parallel GRBM algorithm and floating bodies projection
	// and MPCs (rixen method)
protected:
	int nGrbm = 0;
	int *neighbNumGRBMs = nullptr;

	int numGroupRBM = 0, groupRBMoffset = 0;
	int *neighbNumGroupGrbm = nullptr;
	int *neighbGroupGrbmOffset = nullptr;
	int numGlobalRBMs = 0;
	int *dualToBoundary = nullptr;

	std::unique_ptr<Rbm> rigidBodyModesG;

	int numWIdof = 0;  // number of dofs on the wet interface (both fluid and structure)
	int numWInodes = 0;  // number of nodes on the wet interface (both fluid and structure)
	std::vector<int> wetInterfaceMap;  // dof map
	std::vector<int> wetInterfaceNodeMap;
	std::vector<int> wetInterfaceNodes;
	std::vector<int> numNeighbWIdof;
	std::vector<int> wiInternalMap;
	std::vector<DofSet> wetInterfaceDofs;

	GlobalToLocalMap glToLocalWImap;

	GlobalToLocalMap *neighbGlToLocalWImap = nullptr;

	GlobalToLocalMap glToLocalNode; // This seems to be for coarse problem only.

	/// \brief store indices for possible rebuild (multiple LHS freq sweep)
	int edgeQindex[2] = {-1, -1};

	double k_f = 0.0, k_p = 0.0, k_s = 0.0, k_s2 = 0.0;  // wave numbers for FETI-DPH for this subdomain
	double *neighbK_p = nullptr, *neighbK_s = nullptr, *neighbK_s2 = nullptr, *neighbK_f = nullptr;  // neighbors' wave numbers
	double Ymod, Prat = 0.0, Dens = 0.0, Thih = 0.0, Sspe = 0.0;  // Young's modulus, Poisson ration, density, thickness, speed of sound
	double *neighbYmod = nullptr, *neighbPrat = nullptr, *neighbDens = nullptr, *neighbThih = nullptr, *neighbSspe = nullptr;  // neighbor's values

};

template <typename Scalar>
class _AVMatrix : public Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> {
public:
	_AVMatrix() = default;

	const Scalar *operator[](int i) const { return &(*this)(0, i); }
	Scalar *operator[](int i) { return &(*this)(0, i); }
};

/** \brief Pure Interface of what a the notion of Subdomain provides for FETI solver. */
template <typename Scalar>
class FetiSub : virtual public FetiBaseSub {
public:
	virtual void multMFi(GenSolver<Scalar> *s, Scalar *, Scalar *, int numRHS) const = 0;
	virtual void getQtKQ(GenSolver<Scalar> *s) = 0;
	virtual void getQtKQ(int iMPC, Scalar *QtKQ) = 0;
	Scalar getMpcRhs(int iMPC) const;
	Scalar getMpcRhs_primal(int iMPC) const;
	// TODO Figure out how to make this const.
	void sendDiag(GenSparseMatrix<Scalar> *s, FSCommPattern<Scalar> *vPat);
	void factorKii();
	virtual void sendInterf(const Scalar *interfvec, FSCommPattern<Scalar> *vPat) const = 0;
	virtual void scatterHalfInterf(const Scalar *s, Scalar *loc) const = 0;
	virtual void getHalfInterf(const Scalar *s, Scalar *t) const = 0;
	virtual void getHalfInterf(const Scalar *s, Scalar *t, const Scalar *ss, Scalar *tt) const = 0;
	/** \brief Basic FETI operation.
	 * \details Computes localvec = K-1 (localvec -B interfvec)
	 * then    interfvec = B^T localvec and sends local data to neighbors
	 *
	 * @param s
	 * @param localvec
	 * @param interfvec
	 */
	virtual void fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec) const = 0;
	void fetiBaseOp(Scalar *uc,GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec) const;
	virtual void fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec, Scalar *beta) const = 0;
	void fetiBaseOpCoupled2(const Scalar *uc, const Scalar *localvec, Scalar *interfvec,
	                        FSCommPattern<Scalar> *wiPat, const Scalar *fw = nullptr) const;
	void fetiBaseOpCoupled1(GenSolver<Scalar> *s, Scalar *localvec, const Scalar *interfvec,
	                                FSCommPattern<Scalar> *wiPat) const;
	virtual void multQt(int glMPCnum, const Scalar *V, int numV, Scalar *QtV) const = 0;
	virtual void multQtKBt(int glNumMPC, const Scalar *G, Scalar *QtKBtG, Scalar alpha=1.0, Scalar beta=1.0) const = 0;
	virtual int numRBM() const = 0;
	virtual const Scalar *getQtKpBt() const = 0;
	virtual void split(const Scalar *v, Scalar *v_f, Scalar *v_c) const = 0;

	/// \brief Generate the B matrices for the problem
	void makeBs();

	void makeLocalMpcToDof(); //HB: create the LocalMpcToDof connectivity for a given DofSetArray
	void makeLocalMpcToMpc();
	void updateActiveSet(Scalar *v, double tol, int flag, bool &statusChange);

	// (G^T*G) matrix assembly
	void assembleGtGsolver(GenSparseMatrix<Scalar> *GtGsolver);
	void getLocalMpcForces(double *mpcLambda, DofSetArray *cornerEqs,
	                       int mpcOffset, GenVector<Scalar> &uc);
	void useKrrNullspace();
	// R matrix construction and access
	void makeLocalRstar(FullM **Qtranspose); // this is used by decomposed domain GRBM algorithm
	// R matrix-vector multiplication
	void addRalpha(Scalar *u, GenVector<Scalar> &alpha) const;  // u += R_g*alpha

	int zColDim() { return rigidBodyModesG->Zstar->numCol(); }
	int zRowDim() { return rigidBodyModesG->Zstar->numRow(); }

	void deleteLocalRBMs() { rigidBodyModesG.reset(nullptr); }
	void setBodyRBMoffset(int _boff) { bodyRBMoffset = _boff; }
	void assembleE(GenVector<Scalar> &e, Scalar *f) const; // e = R^T*f
	// G matrix construction and destruction
	void makeG();
	void makeTrbmG(Scalar *rbms, int nrbms, int glNumCDofs);

	void setGCommSize(FSCommStructure *pat) const override;
	void sendG(FSCommPattern<Scalar> *rbmPat);
	void receiveG(FSCommPattern<Scalar> *rbmPat);
	void zeroG();
	void deleteG();

	void getFr(const Scalar *f, Scalar *fr) const;
	void getFc(const Scalar *f, Scalar *fc) const;
	virtual void getFw(const Scalar *f, Scalar *fw) const;
	// G matrix-vector multiplication
	void multG(const GenVector<Scalar> &x, Scalar *y, Scalar alpha) const;  // y = alpha*G*x
	void trMultG(const Scalar *x, GenVector<Scalar> &y, Scalar alpha) const; // y = alpha*G^T*x


	// R_g matrix construction and access
	void buildGlobalRBMs(GenFullM<Scalar> &Xmatrix, const Connectivity *cornerToSub); // use null space of (G^T*P_H*G) ... trbm method !!!
	void getGlobalRBM(int iRBM, Scalar *Rvec) const;
	// R_g matrix-vector multiplication
	void subtractRstar_g(Scalar *u, GenVector<Scalar> &beta) const; // u -= R_g*beta
	void addRstar_gT(Scalar *u, GenVector<Scalar> &beta) const; // u += R_g*beta
	// (R_g^T*R_g) matrix assembly
	void assembleRtR(GenFullM<Scalar> &RtRu);

	virtual void makeKbbMpc() = 0;
	virtual void makeKbb(DofSetArray *dofsetarray=0) = 0;
	void rebuildKbb();

	// new B operators
	void multBr(const Scalar *localvec, Scalar *interfvec, const Scalar *uc = 0, const Scalar *uw = 0) const;
	void multAddBrT(const Scalar *interfvec, Scalar *localvec, Scalar *uw = nullptr) const;

	void multAddCT(const Scalar *interfvec, Scalar *localvec) const;
	void multC(const Scalar *localvec, Scalar *interfvec) const;

	/** \brief Compute \f$ f_r = K_{rc} u_c \f$. */
	void multKrc(Scalar *fr, const Scalar *uc) const;

	void multKcc();
	void multKbb(const Scalar *u, Scalar *Pu, Scalar *delta_u = 0, Scalar * delta_f= 0, bool errorFlag = true);
	void multKbbCoupled(const Scalar *u, Scalar *Pu, Scalar *deltaF, bool errorFlag = true);
	void multDiagKbb(const Scalar *u, Scalar *Pu) const;

	void collectScaling(FSCommPattern<Scalar> *vPat);
	void fScale(Scalar *locF, FSCommPattern<Scalar> *vPat, Scalar *locFw = 0);
	void initMpcScaling();
	void initScaling();
	void weightEdgeGs();
	void makeQ();

	void applyBtransposeAndScaling(const Scalar *u, Scalar *v, Scalar *deltaU = 0, Scalar *localw = 0) const;
	void applyScalingAndB(const Scalar *res, Scalar *Pu, Scalar *localw = 0) const;
	void setMpcDiagCommSize(FSCommStructure *mpcDiagPat) const;
	void sendMpcDiag(FSCommPattern<Scalar> *mpcDiagPat);
	void collectMpcDiag(FSCommPattern<Scalar> *mpcDiagPat);
	void sendMpcScaling(FSCommPattern<Scalar> *mpcPat);
	void collectMpcScaling(FSCommPattern<Scalar> *mpcPat);
	void assembleMpcIntoKcc();

	void makeEdgeVectorsPlus(bool isFluidSub = false, bool isThermalSub = false,
	                         bool isUndefinedSub = false);
	void makeAverageEdgeVectors();
	void extractInterfRBMs(int numRBM, Scalar *locRBMs, Scalar *locInterfRBMs);
	void sendInterfRBMs(int numRBM, Scalar *locInterfRBMs, FSCommPattern<Scalar> *rbmPat);
	void recvInterfRBMs(int iNeighb, int numNeighbRBM, Scalar *neighbInterfRBMs, FSCommPattern<Scalar> *rbmPat);
	void sendInterfaceGrbm(FSCommPattern<Scalar> *rbmPat);
	void receiveInterfaceGrbm(FSCommPattern<Scalar> *rbmPat);
	// templated R and G functions
	// note #1: we use feti to solve global domain problem: min 1/2 u_g^T*K_g*u_g - u_g^T*f_g subj. to C_g*u_g <= g
	//          by solving an equivalent decomposed domain problem: min 1/2 u^T*K*u - u^T*f subj to B*u = 0, C*u <= g
	// the columns of R_g span the left null space of [ K_g & C_gtilda^T // C_gtilda & 0 ] ... C_gtilda = gtilda are the active constraints
	// the columns of R span the left null space of K
	// G = [B^T C^T]^T*R
	// e = R^T*f
	// note #2: the null space of a matrix must be templated (i.e. it is real if the matrix is real or complex if the matrix is complex)
	// note #3: the geometric rigid body modes (GRBMs) or heat zero energy modes (HZEMs) can be used to construct the null space SOMETIMES, NOT ALWAYS !!!!
	// note #4: the GRBMs are not always computed correctly when there are mechanisms
	// note #5: the GRBMs/HZEMs are always real

	void addTrbmRalpha(Scalar *rbms, int nrbms, int glNumCDofs, Scalar *alpha, Scalar *ur) const; // u += R_g*alpha
	void assembleTrbmE(Scalar *rbms, int nrbms, int glNumCDofs, Scalar *e, Scalar *fr) const; // e = R^T*f

	void projectActiveIneq(Scalar *v) const;
	void normalizeCstep1(Scalar *cnorm);
	void normalizeCstep2(Scalar *cnorm);
	void recvMpcStatus(FSCommPattern<int> *mpcPat, int flag, bool &statusChange);

	void mergeUr(Scalar *ur, Scalar *uc, Scalar *u, Scalar *lambda);

	void multfc(const VectorView<Scalar> &fr, /*Scalar *fc,*/ const VectorView<Scalar> &lambda) const;
	void multFcB(Scalar *bf);

	void subtractMpcRhs(Scalar *interfvec);

	void setLocalLambda(Scalar *_localLambda);
	double getMpcError() const;
	void applyMpcSplitting();

	void initMpcStatus();
	void saveMpcStatus();
	void restoreMpcStatus();
	void saveMpcStatus1() const; // const is a lie but we have mutable.
	void saveMpcStatus2();
	void cleanMpcData();

	void constructKcc();

	void makeKccDofsExp2(int nsub, FetiBaseSub **sd, int augOffset,
	                     Connectivity *subToEdge);

	const std::vector<Scalar> &getfc() const { return fcstar; }

	/// \brief Solver for the remainder DOFs.
	std::unique_ptr<GenSolver<Scalar>> Krr;
	/// \brief Sparse view of the solver. Typically used to fill the matrix before calling factor.
	GenSparseMatrix<Scalar>   *KrrSparse = nullptr; //!< Alias to Krr.
	std::unique_ptr<GenSparseMatrix<Scalar>> KiiSparse;
	GenSolver<Scalar>         *KiiSolver = nullptr;
	std::unique_ptr<GenCuCSparse<Scalar> >     Kib;
	std::unique_ptr<GenAssembledFullM<Scalar>> Kcc;
	std::unique_ptr<GenCuCSparse<Scalar>>      Krc;
	std::unique_ptr<GenCuCSparse<Scalar>>      Grc;
	_AVMatrix<Scalar> Ave;
	_AVMatrix<Scalar> Eve;

public:
	std::vector<Scalar> rbms;
	std::vector<Scalar> interfaceRBMs;
	// MPC related data
	mutable std::vector<std::unique_ptr<SubLMPCons<Scalar>>> mpc; // multiple point constraints
	std::vector<std::unique_ptr<SubLMPCons<Scalar>>> mpc_primal;
	std::unique_ptr<GenSparseSet<Scalar>> Src;
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> BKrrKrc;
	std::unique_ptr<GenDBSparseMatrix<Scalar>> Kbb;    //!< Boundary to boundary stiffness matrix.

	GenVector<Scalar> diagCCt;
protected:
	GenVector<Scalar> scaling;
	mutable GenVector<Scalar> deltaFmpc;
	GenVector<Scalar> deltaFwi;
	GenVector<Scalar> wweight;
	GenVector<Scalar> kweight; //!< stiffness weights (i.e. sum of Kii for all subd sharing that dof)

	// templated RBMs
	GenFullM<Scalar> Rstar;
	GenFullM<Scalar> Rstar_g;
	std::unique_ptr<GenFullM<Scalar>> sharedRstar_g;
	std::unique_ptr<GenFullM<Scalar>> tmpRstar_g;
	std::vector<std::unique_ptr<GenFullM<Scalar>>> G;
	std::vector<std::unique_ptr<GenFullM<Scalar>>> neighbG;

	int bodyRBMoffset = 0;
	std::unique_ptr<Rbm> rigidBodyModes;

	mutable std::vector<Scalar> fcstar; // TODO Move this out!

	std::unique_ptr<GenFsiSparse<Scalar>> neighbKww;
	mutable std::vector<Scalar> localw_copy;
	// coupled_dph
	std::unique_ptr<GenDBSparseMatrix<Scalar>> Kww;
	std::unique_ptr<GenCuCSparse<Scalar>>      Kcw;
	std::unique_ptr<GenMpcSparse<Scalar>> Kcw_mpc;
	std::unique_ptr<GenCuCSparse<Scalar>> Krw;
	mutable std::vector<Scalar> localw;

	Eigen::SparseMatrix<double> B, Bw;
	Eigen::SparseMatrix<Scalar> Bm, Bc;

};

#endif //FEM_FETUSUB_H
