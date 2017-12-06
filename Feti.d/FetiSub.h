//
// Created by Michel Lesoinne on 11/3/17.
//

#ifndef FEM_FETUSUB_H
#define FEM_FETUSUB_H
#include <vector>

class FSCommStructure;
template <typename Scalar>
class FSCommPattern;
template <typename Scalar>
class GenSolver;
template <typename Scalar>
class GenSparseMatrix;
template <typename Scalar>
class GenAssembledFullM;
template <typename Scalar>
class GenCuCSparse;

class Connectivity;
class CoordSet;
class DofSet;

/** \brief Pure Interface of what a the notion of Subdomain provides for FETI solver. */
class FetiBaseSub {
public:
	/** \brief Obtain the number of neighbor subdomains. */
	virtual int numNeighbors() const = 0;
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

	virtual void setMpcDiagCommSize(FSCommStructure *mpcDiagPat) const = 0;

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

	const std::vector<int> &getCornerNodes() const { return glCornerNodes; }
	std::vector<int> &getCornerNodes() { return glCornerNodes; }
	void markCornerDofs(int *glCornerDofs) const;
	/* missing:
	 * splitInterf
	 * setMpcNeighbCommSize
	 * setGCommSize
	 * assembleTrbmE
	 * assembleE
	 */
	std::vector<DofSet> cornerDofs;
	DofSet *edgeDofs;      // JAT 112113

protected:
	/// \brief Corner nodes in local numbering.
	std::vector<int> cornerNodes;
	std::vector<bool> isCornerNode;   // true for node which is a corner node; false otherwise
	std::vector<int> glCornerNodes; // corner nodes in global numbering
	int numCRN = 0;
	int numCRNdof = 0;
};

/** \brief Pure Interface of what a the notion of Subdomain provides for FETI solver. */
template <typename Scalar>
class FetiSub : virtual public FetiBaseSub {
public:
	virtual void multMFi(GenSolver<Scalar> *s, Scalar *, Scalar *, int numRHS) const = 0;
	virtual void getQtKQ(GenSolver<Scalar> *s) = 0;
	virtual void getQtKQ(int iMPC, Scalar *QtKQ) = 0;
	virtual Scalar getMpcRhs(int iMPC) const = 0;
	virtual Scalar getMpcRhs_primal(int iMPC) const = 0;
	// TODO Figure out how to make this const.
	virtual void sendDiag(GenSparseMatrix<Scalar> *s, FSCommPattern<Scalar> *vPat) = 0;
	virtual void factorKii() = 0;
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
	virtual void fetiBaseOp(Scalar *uc,GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec) const = 0;
	virtual void fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec, Scalar *beta) const = 0;
	virtual void fetiBaseOpCoupled2(const Scalar *uc, const Scalar *localvec, Scalar *interfvec,
	                                FSCommPattern<Scalar> *wiPat, const Scalar *fw = nullptr) const = 0;
	virtual void fetiBaseOpCoupled1(GenSolver<Scalar> *s, Scalar *localvec, const Scalar *interfvec,
	                                FSCommPattern<Scalar> *wiPat) const = 0;
	virtual void multQt(int glMPCnum, const Scalar *V, int numV, Scalar *QtV) const = 0;
	virtual void multQtKBt(int glNumMPC, const Scalar *G, Scalar *QtKBtG, Scalar alpha=1.0, Scalar beta=1.0) const = 0;
	virtual int numRBM() const = 0;
	virtual const Scalar *getQtKpBt() const = 0;
	virtual void split(const Scalar *v, Scalar *v_f, Scalar *v_c) const = 0;
	// Missing:
	/*
	 * split
	 * addRstar_gT
	 * subtractRstar_g
	 * addRalpha
	 getfc
	 multKrc
	 multfc
	 multFcB
	 multG
	 multAddCT
	 trMultG
	 getMpcError
	 subtractMpcRhs
	 getLocalMpcForces
	 getMpcRhs_primal
	 setBodyRBMoffset
	 getGlobalRBM
	 setWI...
	 projectActivIneq
	 normalizeCstep1
	 assembleRtR
	 * */
	/// \brief Solver for the remainder DOFs.
	std::unique_ptr<GenSolver<Scalar>> Krr;
	/// \brief Sparse view of the solver. Typically used to fill the matrix before calling factor.
	GenSparseMatrix<Scalar>   *KrrSparse = nullptr; //!< Alias to Krr.
	std::unique_ptr<GenAssembledFullM<Scalar>> Kcc;
	std::unique_ptr<GenCuCSparse<Scalar>>      Krc;
	std::unique_ptr<GenCuCSparse<Scalar>>      Grc;
};

#endif //FEM_FETUSUB_H
