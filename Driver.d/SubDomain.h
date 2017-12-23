#ifndef _SUB_DOMAIN_H_
#define _SUB_DOMAIN_H_

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/ControlLawInfo.h>
#include <Feti.d/DistrVector.h>
#include <Feti.d/FetiSub.h>
#include <Corotational.d/Corotator.h>
#include <Math.d/DistVector.h>
#include <Utils.d/MyComplex.h>
#include <Math.d/FsiSparse.h>
#include <Driver.d/SComm.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/GlobalToLocalMap.h>
#include <Utils.d/MathUtils.h>
#include <vector>
#include <list>

extern GeoSource *geoSource;

template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
template <class Scalar> class GenDistrDomain;
typedef GenDistrDomain<double> DistrDomain;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
template <class Scalar> class GenFetiSolver;
typedef GenFetiSolver<double> FetiSolver;
template <class Scalar> class GenFetiDPSolver;
typedef GenFetiDPSolver<double> FetiDPSolver;
template <class Scalar> class GenFetiOp;
typedef GenFetiOp<double> FetiOp;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
class IntFullM;
class GeomState;
class DistrGeomState;
template <class Scalar> class GenBLKSparseMatrix;
typedef GenBLKSparseMatrix<double> BLKSparseMatrix;
template <class Scalar> class GenSkyMatrix;
typedef GenSkyMatrix<double> SkyMatrix;
template <class Scalar> class GenDomainGroupTask;
typedef GenDomainGroupTask<double> DomainGroupTask;
template <class Scalar> class GenAssembledFullM;
typedef GenAssembledFullM<double> AssembledFullM;
template <class Scalar> class GenSparseSet;
typedef GenSparseSet<double> SparseSet;
class SubCornerHandler;
template <class Type> class FSCommPattern;
class DistrComplexVector;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenMpcSparse;

class BaseSub : virtual public Domain , virtual public FetiBaseSub
{
protected:
	int subNumber;
	int localSubNumber; // relevant when running in distributed

	GlobalToLocalMap glToLocalNode;
	GlobalToLocalMap glToLocalElem;
	int *glNums = nullptr;
	int *glElems = nullptr;
	int glNumNodes;
	std::vector<int> weight; ///!< \brief DOF weights (i.e. number of subd sharing that dof).
	double *bcx = nullptr;
	DComplex *bcxC = nullptr; // FETI-H
	double *vcx = nullptr, *acx = nullptr;
	int *locToGlSensorMap = nullptr;
	int *locToGlActuatorMap = nullptr;
	int *locToGlUserDispMap = nullptr;
	int *locToGlUserForceMap = nullptr;
	int boundLen = 0;
	int *boundMap = nullptr;
	int *dualToBoundary = nullptr;
	int internalLen = 0;
	int *internalMap = nullptr;
	int crnDofSize = 0;
	int *crnPerNeighb = nullptr;
	long memK = 0;       // memory necessary to store K(s)
	long memPrec = 0;    // memory necessary to store Preconditioner
#ifdef DISTRIBUTED
	// for distributed output of single nodes
    int numNodalOutput = 0;
    int *outputNodes = nullptr;
    int *outIndex = nullptr;
#endif
	int globalNMax;  // highest global node number of all the nodes in this subdomain
	int globalEMax;  // highest global element number of all the elements in this subdomain


public:
	BaseSub();
	BaseSub(Domain &dom, int sn, Connectivity &con, Connectivity &nds, int gn);
	BaseSub(Domain &dom, int sn, int nNodes, int *nds,
	        int nElems, int *elems, int gn);
	BaseSub(Domain &dom, int sn, CoordSet* nodes, Elemset* elems, int *glNodeNums,
	        int *glElemNums, int gn); // PJSA: for new sower
	virtual ~BaseSub();

	const FetiInfo &getFetiInfo() const override { return solInfo().getFetiInfo(); }

//	int *glCrnGroup;    // group of each corner node (global numbering)
	int nGrbm = 0;
	DofSet **boundaryDOFs = nullptr;
	int nCDofs = -1;
	int *neighbNumGRBMs = nullptr;
	int *edgeDofSizeTmp = nullptr;   // XXXX
	double k_f = 0.0, k_p = 0.0, k_s = 0.0, k_s2 = 0.0;  // wave numbers for FETI-DPH for this subdomain
	double *neighbK_p = nullptr, *neighbK_s = nullptr, *neighbK_s2 = nullptr, *neighbK_f = nullptr;  // neighbors' wave numbers
	double Ymod, Prat = 0.0, Dens = 0.0, Thih = 0.0, Sspe = 0.0;  // Young's modulus, Poisson ration, density, thickness, speed of sound
	double *neighbYmod = nullptr, *neighbPrat = nullptr, *neighbDens = nullptr, *neighbThih = nullptr, *neighbSspe = nullptr;  // neighbor's values

	int dofWeight(int i) { return weight[i]; }
	int crnDofLen() const  { return crnDofSize; }
	IntFullM* getC(int &crnDofSize, FSCommPattern<int> *sPat);
	void showExchangeData();
	void countCornerDofs(int *cWeight);
	void applyAuxData();
	void distributeBCs(int *);
	void setControlData(ControlLawInfo *_claw, int *, int *, int *, int *);
#ifdef DISTRIBUTED
	void setOutputNodes(int, int *, int *);
	int *getOutputNodes()       { return outputNodes; }
	int *getOutIndex()          { return outIndex; }
	int getNumNodalOutput() const { return numNodalOutput; }
#endif
	int getBC(BCond *, int, int *, BCond *&);
	void setGlNodes(int *globalNodeNums) { glNums = globalNodeNums; }
	const int *getGlNodes() const { return glNums; }
	int *getGlElems() const    { return glElems; }
	int *getGlMPCs()  const     { return localToGlobalMPC; }
	int glToPackElem(int e) const { return (geoSource->glToPackElem(e) > globalEMax) ? -1 : glToLocalElem[geoSource->glToPackElem(e)]; }
	int *getSensorDataMap() const { return locToGlSensorMap; }
	int *getUserDispDataMap() const { return locToGlUserDispMap; }
	int countElemNodes();
	int numMPCs() const override { return numMPC; }
	int numMPCs_primal() const  { return numMPC_primal; }
	int globalNumNodes();
	int numNodes() const        { return numnodes; }
	int findProfileSize();
	int renumberBC(int *);
	void makeGlobalToLocalNodeMap();
	void makeGlobalToLocalElemMap();
	int globalToLocal(int i)    { return (i < 0 || i > globalNMax) ? -1 : glToLocalNode[i]; }
	GlobalToLocalMap &getGlobalToLocalNode() { return glToLocalNode; }
	int localToGlobal(int i)    { return glNums[i]; }
	int globalToLocalElem(int i) { return (i < 0 || i > globalEMax) ? -1 : glToLocalElem[i]; }
	int localToGlobalElem(int i) { return glElems[i]; }
	int getGlobalNMax()         { return globalNMax; }
	int* makeBMaps(const DofSetArray *dofsetarray=0);
	int* makeIMaps(const DofSetArray *dofsetarray=0);
	int subNum() const override  { return subNumber; }
	int localSubNum() const override { return localSubNumber; }
	int getNumUncon() const override { return numUncon(); }
	int localLen() const override { return (cc_dsa) ? cc_dsa->size() : c_dsa->size(); }
	ConstrainedDSA * getCCDSA()  { return (cc_dsa) ? cc_dsa.get() : c_dsa; }
	int localRLen() const override { return cc_dsa->size(); }
	void sendNumNeighbGrbm(FSCommPattern<int> *pat);
	void recvNumNeighbGrbm(FSCommPattern<int> *pat);
	void putNumMPC(int *ptr) { ptr[subNumber] = numMPC; }
	void putLocalToGlobalMPC(int *ptr, int *tg) { for(int i=0; i<numMPC; ++i) tg[ptr[subNumber]+i] = localToGlobalMPC[i]; }
	void putNumMPC_primal(int *ptr) { ptr[subNumber] = numMPC_primal; }
	void putLocalToGlobalMPC_primal(int *ptr, int *tg) { for(int i=0; i<numMPC_primal; ++i) tg[ptr[subNumber]+i] = localToGlobalMPC_primal[i]; }

	void makeLocalToGroupMPC(Connectivity *groupToMPC);
	void findEdgeNeighbors();
	void makeMpcInterface(Connectivity *subToMpc, const Connectivity &lmpcToSub,
	                      Connectivity *subToSub_mpc);
	void makeFsiInterface(Connectivity *subToFsi, Connectivity *fsiToSub,
	                      Connectivity *subToSub_fsi);

	bool checkForColinearCrossPoints(int numCornerPoints, int *localCornerPoints);
	void addCornerPoints(int *glCornerList);
	int numCorners() const override { return numCRN; }

	int numCornerDofs()	const { return numCRNdof; }
	int numCoarseDofs();
	int nCoarseDofs()  const { return nCDofs; }

	ConstrainedDSA *get_c_dsa() const { return c_dsa; }
	const std::vector<int> &getWeights() const { return weight; }

public:
	/// \copydoc
	int getLocalMPCIndex(int globalMpcIndex) const override;
	/// \copydoc
	int getGlobalMPCIndex(int localMpcIndex) const override;
	void makeLocalMpcToGlobalMpc(Connectivity *mpcToMpc);
	void setLocalMpcToBlock(Connectivity *mpcToBlock, Connectivity *blockToMpc);
	void setGroup(Connectivity *subToGroup) { this->group = (*subToGroup)[subNumber][0]; }
	void setNumGroupRBM(int *ngrbmGr);
	void getNumGroupRBM(int *ngrbmGr);
	void addNodeXYZ(double *centroid, double* nNodes);
	void sendNeighbGrbmInfo(FSCommPattern<int> *pat);
	void receiveNeighbGrbmInfo(FSCommPattern<int> *pat);
	void setCommSize(FSCommStructure *pat, int size) const override;
	void setMpcNeighbCommSize(FSCommPattern<int> *pt, int size) const override;

public:
	void setCorners(int nCorners, int *crnList);

	bool isEdgeNeighbor(int neighb) const { return scomm->isEdgeNeighb[neighb]; }
	int interfLen() const override; //<! \brief Total length for the local interface
	int halfInterfLen() const override; //<! \brief Length of the "half interface"
	void computeMasterFlag(const Connectivity &mpcToSub) override;
	bool* getMasterFlag() { return masterFlag; }
	const bool* getMasterFlag() const override { return masterFlag; }
	const bool* getInternalMasterFlag();

protected:
	void computeInternalMasterFlag();

public:
	void setNodeCommSize(FSCommStructure *, int d = 1) const ;
	/// \copydoc
	void setDofCommSize(FSCommStructure *) const override;
	void setDofPlusCommSize(FSCommStructure *) const;
	void setRbmCommSize(int numRBM, FSCommStructure *) const override;
	void setMpcCommSize(FSCommStructure *mpcPat) const override;

	// for timing file
	double getSharedDofCount();
	int getTotalDofCount();
	SubCornerHandler *getCornerHandler();

	void initHelm(Domain &dom);

	// DPH functions
	int isFluid(int i=0);
	void setWaveNumbers(double *waveNumbers);
	void computeWaveNumbers();
	void sendWaveNumbers(FSCommPattern<double> *kPat);
	void collectWaveNumbers(FSCommPattern<double> *kPat);
	void getDirections(int numDirec, int numWaves, double *&wDir_x, double *&wDir_y, double *&wDir_z);
	void getOneDirection(double d, int i, int j, int k, int &nnum, int numWaves, double *wDir_x, double *wDir_y, double *wDir_z);
	void getOneDirection(double x, double y, double z, int &nnum, int numWaves, double *wDir_x, double *wDir_y, double *wDir_z);
	void getDirections13(int numDirec, double *wDir_x, double *wDir_y, double *wDir_z);
	void GramSchmidt(double *Q, bool *isUsed, DofSet desired, int nQPerNeighb, bool isPrimalAugmentation);
	void averageMatProps();
	void sendMatProps(FSCommPattern<double> *matPat);
	void collectMatProps(FSCommPattern<double> *matPat);

	void setDirichletBC(std::list<BCond *> *_list);
	void setNeumanBC(std::list<BCond *> *_list);
	void setInitialDisplacement(std::list<BCond *> *_list);
	void setInitialDisplacement6(std::list<BCond *> *_list);
	void setInitialVelocity(std::list<BCond *> *_list);
	void setSensor(std::list<BCond *> *_list);
	void setActuator(std::list<BCond *> *_list);
	void setUsdd(std::list<BCond *> *_list);
	void setUsdf(std::list<BCond *> *_list);
	void setClaw(char* _fileName, char* _routineName) {
		claw = new ControlLawInfo;
		claw->fileName = _fileName;
		claw->routineName = _routineName;
	}
	void setComplexDirichletBC(std::list<ComplexBCond *> *_list);
	void setComplexNeumanBC(std::list<ComplexBCond *> *_list);
	void setDnb(std::list<SommerElement *> *_list);
	void setScat(std::list<SommerElement *> *_list);
	void setArb(std::list<SommerElement *> *_list);
	// coupled_dph
protected:
	bool *wetInterfaceMark = nullptr;
	bool *wetInterfaceFluidMark = nullptr;
	bool *wetInterfaceStructureMark = nullptr;
	bool *wetInterfaceCornerMark = nullptr;

	int *wDofToNode = nullptr; //HB
	DofSet *wetInterfaceDofs = nullptr;
	Connectivity *drySharedNodes = nullptr;
	bool *wiMaster = nullptr;
	GlobalToLocalMap *neighbGlToLocalWImap = nullptr;
	int numFsiNeighb;
	int *fsiNeighb = nullptr;
	int *wiInternalMap = nullptr;
	bool isMixedSub = false;
	int edgeQindex[2] = {-1, -1};
	double prev_cscale_factor;

public:
	Connectivity *nodeToSub;
	void setnodeToSubConnectivity(Connectivity *nTsubConn) { nodeToSub = nTsubConn; }
	void markWetInterface(int nWI, int *wiNum);
	bool onWetInterface(int iNode) { return wetInterfaceMark[iNode]; }
	bool onWetInterfaceFluid(int iNode) { return wetInterfaceFluidMark[iNode]; }
	bool onWetInterfaceStructure(int iNode) { return wetInterfaceStructureMark[iNode]; }
	bool isWetInterfaceCorner(int iNode) { return wetInterfaceCornerMark[iNode]; }
	void setWetInterface(int nWI, int *wiNum);
	void sendNumWIdof(FSCommPattern<int> *sPat) const;
	void recvNumWIdof(FSCommPattern<int> *sPat);
	void sendWImap(FSCommPattern<int> *pat);
	void recvWImap(FSCommPattern<int> *pat);
	void makeDSA();
	void makeCDSA();
	void makeCCDSA();
	int numWetInterfaceDofs() const { return numWIdof; }
	GlobalToLocalMap& getGlToLocalWImap() { return glToLocalWImap; }
	GlobalToLocalMap& getNeighbGlToLocalWImap(int i) { return neighbGlToLocalWImap[i]; }
	void zeroEdgeDofSize();
	void mergeInterfaces();

#ifdef HB_COUPLED_PRECOND
	Connectivity* precNodeToNode;
#endif
};


template<class Scalar>
class GenSubDomain : public BaseSub , public FetiSub<Scalar>
{
private:
	Scalar *kweight; // stiffness weights (i.e. sum of Kii for all subd sharing that dof)
	mutable Scalar *deltaFmpc;
	int *cornerWeight;
	void applyBtransposeAndScaling(const Scalar *u, Scalar *v, Scalar *deltaU = 0, Scalar *localw = 0) const;
	void applyScalingAndB(const Scalar *res, Scalar *Pu, Scalar *localw = 0) const;
	void initialize();

protected:
	Scalar *scaling;
	void sendDOFList(FSCommPattern<int> *pat); // Send to neighbors the list of DOFs on the shared nodes

public:

	Scalar                    *rbms;
	Scalar                    *interfaceRBMs;
	GenFullM<Scalar>          *qtkq;
	std::unique_ptr<GenSparseMatrix<Scalar>> KiiSparse;
	GenSolver<Scalar>         *KiiSolver;
	std::unique_ptr<GenCuCSparse<Scalar> >     Kib;
	GenSparseMatrix<Scalar>   *MPCsparse;
	std::unique_ptr<GenDBSparseMatrix<Scalar>> Kbb;    // for preconditioning
	Corotator           	    **corotators;
	Scalar                    *QtKpBt;

	int *glBoundMap;
	int *glInternalMap;

private:
	GenSolver<Scalar> *localCCtsolver;
	GenSparseMatrix<Scalar> *localCCtsparse;
	Scalar *diagCCt;
	int lengthCCtData;
	int *CCtrow, *CCtcol;
	Scalar* CCtval;
	Scalar *bcx_scalar;




public:
	GenSubDomain(int, int);
	GenSubDomain(Domain &, int sn, Connectivity &con, Connectivity &nds, int gn);
	GenSubDomain(Domain &, int sn, int nNodes, int *nds,
	             int nElems, int *elems, int gn);
	GenSubDomain(Domain &, int sn, CoordSet* nodes, Elemset* elems, int *glNodeNums,
	             int *glElemNums, int gn); // PJSA: for new sower
	~GenSubDomain();

	long getMemoryK() { Scalar s; return memK*long(sizeof(s))/long(1024); /* PJSA return memory usage in KB */ }
	long getMemoryPrec() { Scalar s; return memPrec*long(sizeof(s))/long(1024); /* PJSA return memory usage in KB */ }
	void extractControlData(Scalar *, Scalar *, Scalar *,
	                        Scalar *, Scalar *, Scalar *);
	void addUserForce(Scalar *, Scalar *);
	void addCtrl(Scalar *, Scalar *);
	Scalar *getBcx()  { if(!bcx_scalar) makeBcx_scalar(); return bcx_scalar; }
	double *getVcx()  { return vcx; }
	double *getAcx()  { return acx; }
	void setUserDefBC(double *, double *, double *, bool nlflag);
	void reBuildKbb(FullSquareMatrix *kel);
	void addDMass(int glNum, int dof, double m);
	// computes localvec = K-1 (localvec -B interfvec)
	// then    interfvec = B^T localvec and sends local data to neighbors
	void fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec) const override;
	void fetiBaseOp(Scalar *uc,GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec) const;
	void fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec, Scalar *beta) const override;
	void interfaceJump(Scalar *iterfData, FSCommPattern<Scalar> *vPat) const;
	void sendInterf(const Scalar *interfvec, FSCommPattern<Scalar> *vPat) const override;
	void extractAndSendInterf(const Scalar *subvec, FSCommPattern<Scalar> *pat) const;
	void assembleInterf(Scalar *subvec, FSCommPattern<Scalar> *pat) const;
	void splitInterf(Scalar *subvec) const;
	void assembleInterfInvert(Scalar *subvec, FSCommPattern<Scalar> *pat) const;
	void getHalfInterf(const Scalar *s, Scalar *t) const override;
	void getHalfInterf(const Scalar *s, Scalar *t, const Scalar *ss, Scalar *tt) const override;
	void scatterHalfInterf(const Scalar *s, Scalar *loc) const override;
	void rebuildInterf(Scalar *v, FSCommPattern<Scalar> *vPat) const;
	void renumberElements();
	void renumberElementsGlobal();
	void renumberSharedNodes();
	void renumberDirichlet();
	void renumberBCsEtc();
	void renumberControlLaw();
	void renumberMPCs();
	void extractInterfRBMs(int numRBM, Scalar *locRBMs, Scalar *locInterfRBMs);
	void sendInterfRBMs(int numRBM, Scalar *locInterfRBMs, FSCommPattern<Scalar> *rbmPat);
	void recvInterfRBMs(int iNeighb, int numNeighbRBM, Scalar *neighbInterfRBMs, FSCommPattern<Scalar> *rbmPat);
	void sendNode(Scalar (*subvec)[11], FSCommPattern<Scalar> *pat);
	void collectNode(Scalar (*subvec)[11], FSCommPattern<Scalar> *pat);

	void getSRMult(const Scalar *lvec, const Scalar *lbvec, int nRBM, const double *locRBMs, Scalar *alpha) const;
	void sendInterfaceGrbm(FSCommPattern<Scalar> *rbmPat);
	void receiveInterfaceGrbm(FSCommPattern<Scalar> *rbmPat);
	void sendDeltaF(const Scalar *deltaF, FSCommPattern<Scalar> *vPat);
	double collectAndDotDeltaF(Scalar *deltaF, FSCommPattern<Scalar> *vPat);
	void makeKbbMpc();
	void rebuildKbb();
	void makeKbb(DofSetArray *dofsetarray=0);
	void factorKii() override;
	void multKbb(const Scalar *u, Scalar *Pu, Scalar *delta_u = 0, Scalar * delta_f= 0, bool errorFlag = true);
	void multDiagKbb(const Scalar *u, Scalar *Pu) const;
	void multFi(GenSolver<Scalar> *s, Scalar *, Scalar *);
	void multMFi(GenSolver<Scalar> *s, Scalar *, Scalar *, int numRHS) const override;
	void assembleLocalComplexEls(GenSparseMatrix<Scalar> *Kas, GenSolver<Scalar> *smat = 0);
	void mergePrimalError(Scalar* error, Scalar* primal);
	void mergeStress(Scalar *stress, Scalar *weight,
	                 Scalar *globStress, Scalar *globWeight, int glNumNodes);
	void mergeElemStress(Scalar *loc, Scalar *glob, Connectivity *);
	void mergeDisp(Scalar (*xyz)[11], GeomState* locGS, Scalar (*xyz_loc)[11] = NULL);
	void mergeAllDisp(Scalar (*xyz)[11], Scalar *d, Scalar (*xyz_loc)[11] = NULL);
	void mergeAllVeloc(Scalar (*xyz)[11], Scalar *v, Scalar (*xyz_loc)[11] = NULL);
	void mergeAllAccel(Scalar (*xyz)[11], Scalar *a, Scalar (*xyz_loc)[11] = NULL);
	void forceContinuity(Scalar *locdisp, Scalar (*xyz)[11]);
	void mergeDistributedNLDisp(Scalar (*xyz)[11], GeomState* u, Scalar (*xyz_loc)[11] = NULL);
	void mergeForces(Scalar (*mergedF)[6], Scalar *subF);
	void mergeReactions(Scalar (*mergedF)[11], Scalar *subF);
	void mergeDistributedForces(Scalar (*mergedF)[6], Scalar *subF);
	void mergeDistributedReactions(Scalar (*mergedF)[11], Scalar *subF);
	void mergeElemProps(double* props, double* weights, int propType);
	void sendExpDOFList(FSCommPattern<int> *pat);
	template<class Scalar1> void dispatchNodalData(FSCommPattern<Scalar> *pat, NewVec::DistVec<Scalar1> *);
	template<class Scalar1> void addNodalData(FSCommPattern<Scalar> *pat, NewVec::DistVec<Scalar1> *);
	void dispatchInterfaceGeomState(FSCommPattern<double> *geomStatePat, GeomState *geomState);
	void collectInterfaceGeomState(FSCommPattern<double> *geomStatePat, GeomState *geomState);
	void dispatchInterfaceGeomStateDynam(FSCommPattern<double> *geomStatePat, GeomState *geomState);
	void collectInterfaceGeomStateDynam(FSCommPattern<double> *geomStatePat, GeomState *geomState);
	void dispatchInterfaceNodalInertiaTensors(FSCommPattern<double> *pat);
	void collectInterfaceNodalInertiaTensors(FSCommPattern<double> *pat);
	void dispatchGeomStateData(FSCommPattern<double> *, GeomState *);
	void collectGeomStateData(FSCommPattern<double> *, GeomState *);
	void computeElementForce(int, Scalar *u, int Findex, Scalar *force);
	void computeStressStrain(int, Scalar *u, int Findex,
	                         Scalar *stress, Scalar *weight = 0);
	void computeStressStrain(GeomState *gs, Corotator **allCorot,
	                         int, int Findex, Scalar *glStress, Scalar *glWeight = 0, GeomState *refState = NULL);
	void initScaling();
	void sendDiag(GenSparseMatrix<Scalar> *s, FSCommPattern<Scalar> *vPat) override;
	void collectScaling(FSCommPattern<Scalar> *vPat);
	void fSend(const Scalar *locF, FSCommPattern<Scalar> *vPat, Scalar *locFw = 0) const;
	void fScale(Scalar *locF, FSCommPattern<Scalar> *vPat, Scalar *locFw = 0);
	void fSplit(Scalar *locF);
	void updatePrescribedDisp(GeomState *geomState, Scalar deltaLambda);
	void updatePrescribedDisp(GeomState *geomState);
	Scalar displacementNorm(Scalar *displacement);
	void firstAssemble(GenSparseMatrix<Scalar> *K);
	void initMpcScaling();
	void makeZstarAndR(double *centroid);  // makes Zstar and R
	void makeKccDofsExp2(int nsub, GenSubDomain<Scalar> **sd, int augOffset,
	                     Connectivity *subToEdge);
//  void makeKccDofsExp2(int nsub, GenSubDomain<Scalar> **sd);
	void makeKccDofs(DofSetArray *cornerEqs, int augOffset, Connectivity *subToEdge, int mpcOffset = 0);
	void deleteKcc();
	void multKbbMpc(const Scalar *u, Scalar *Pu, Scalar *deltaU, Scalar *deltaF, bool errorFlag = true);
	void getQtKQ(GenSolver<Scalar> *s) override;
	void getQtKQ(int iMPC, Scalar *QtKQ) override;
	void multQt(int glMPCnum, const Scalar *V, int numV, Scalar *QtV) const override;
	void multQt(int glMPCnum, const Scalar *x, Scalar *result) const;
	void multQtKBt(int glNumMPC, const Scalar *G, Scalar *QtKBtG, Scalar alpha=1.0, Scalar beta=1.0) const override;
	void gatherDOFList(FSCommPattern<int> *pat);
	void gatherDOFListPlus(FSCommPattern<int> *pat);

	const Scalar *getQtKpBt() const override { return QtKpBt; }

	friend class GenDistrDomain<Scalar>;
	friend class GenDecDomain<Scalar>;
	friend class GenFetiOp<Scalar>;
	friend class GenFetiSolver<Scalar>;

	void makeQ();
	void precondGrbm();
	void setMpcSparseMatrix();
	void getFw(const Scalar *f, Scalar *fw) const override;
	int numRBM() const { return nGrbm; }
	void makeEdgeVectorsPlus(bool isFluidSub = false, bool isThermalSub = false,
	                         bool isUndefinedSub = false);
	void makeAverageEdgeVectors();
	void weightEdgeGs();
	void constructKcc();
	void constructKrc();
	void initSrc();
	void clean_up();

	// MPC and contact functions
	void extractMPCs(int glNumMPC, ResizeArray<LMPCons *> &lmpc);
	void extractMPCs_primal(int glNumMPC, ResizeArray<LMPCons *> &lmpc);
	void printLMPC();
	void applySplitting();
	void applyDmassSplitting();
	void applyForceSplitting();
	void applyMpcSplitting();
	void constraintProduct(int num_vect, const double* R[], Scalar** V, int trans);
	void addConstraintForces(std::map<std::pair<int,int>, double> &mu, std::vector<double> &lambda, GenVector<Scalar> &f);
	void addCConstraintForces(std::map<std::pair<int,int>, double> &mu, std::vector<double> &lambda, GenVector<Scalar> &fc, double s);
	void locateMpcDofs();
	void deleteMPCs();

	void split(const Scalar *v, Scalar *v_f, Scalar *v_c) const override;
	void bmpcQualify(std::vector<LMPCons *> *bmpcs, int *pstatus, int *nstatus);
	void assembleGlobalCCtsolver(GenSolver<Scalar> *CCtsolver, SimpleNumberer *mpcEqNums);
	void computeSubContributionToGlobalCCt(SimpleNumberer *mpcEqNums); //HB: only compute the subdomain contribution to global CCt
	void assembleGlobalCCtsolver(GenSolver<Scalar> *CCtsolver); //HB: add the subdomain contributions to global CCt
	void assembleBlockCCtsolver(int iBlock, GenSolver<Scalar> *CCtsolver, SimpleNumberer *blockMpcEqNums);
	void constructLocalCCtsolver();
	void assembleLocalCCtsolver();
	void setCCtCommSize(FSCommPattern<Scalar> *cctPat);
	void sendNeighbCCtsolver(FSCommPattern<Scalar> *cctPat, Connectivity *mpcToSub);
	void recNeighbCCtsolver(FSCommPattern<Scalar> *cctPat, Connectivity *mpcToSub);
	void factorLocalCCtsolver();
	void zeroLocalCCtsolver();
	void deleteLocalCCtsolver() { if(localCCtsolver) delete localCCtsolver; }
	void setMpcDiagCommSize(FSCommStructure *mpcDiagPat) const override;
	void sendMpcDiag(FSCommPattern<Scalar> *mpcDiagPat);
	void collectMpcDiag(FSCommPattern<Scalar> *mpcDiagPat);

	void extractMpcResidual(Scalar *subv, GenVector<Scalar> &mpcv, SimpleNumberer *mpcEqNums);
	void insertMpcResidual(Scalar *subv, GenVector<Scalar> &mpcv, SimpleNumberer *mpcEqNums);
	void solveLocalCCt(Scalar *subv);
	void extractBlockMpcResidual(int block, Scalar *subv, GenVector<Scalar> *mpcv,
	                             SimpleNumberer *blockMpcEqNums);
	void insertBlockMpcResidual(Scalar *subv, GenVector<Scalar> **mpcv, Connectivity *mpcToBlock,
	                            SimpleNumberer **blockMpcEqNums);
	void sendMpcInterfaceVec(FSCommPattern<Scalar> *mpcPat, Scalar *interfvec);
	void combineMpcInterfaceVec(FSCommPattern<Scalar> *mpcPat, Scalar *interfvec);
	void sendMpcScaling(FSCommPattern<Scalar> *mpcPat);
	void collectMpcScaling(FSCommPattern<Scalar> *mpcPat);
	void sendMpcStatus(FSCommPattern<int> *mpcPat, int flag);
	void printMpcStatus();
	void initMpcStatus();
	void saveMpcStatus();
	void restoreMpcStatus();
	void saveMpcStatus1();
	void saveMpcStatus2();
	void cleanMpcData();
	void computeContactPressure(Scalar *globStress, Scalar *globWeight);
	void computeLocalContactPressure(Scalar *stress, Scalar *weight);

	void getConstraintMultipliers(std::map<std::pair<int,int>,double> &mu, std::vector<double> &lambda);
	void getLocalMultipliers(std::vector<double> &lambda);
	void setMpcRhs(Scalar *interfvec, double t, int flag);
	void updateMpcRhs(Scalar *interfvec);

	const CoordSet &getNodeSet() const override { return getNodes(); }

protected:
	double *mpcForces;


	GenFsiSparse<Scalar> *neighbKww;
	mutable std::vector<Scalar> localw_copy;
	Scalar Bcx(int i);
	void makeBcx_scalar();
	Scalar *deltaFwi;
	bool isWetInterfaceNode(int n) { return (wetInterfaceNodeMap.size() > 0) ? (wetInterfaceNodeMap[n] > -1) : false; }
	Scalar *wweight;
#ifdef HB_COUPLED_PRECOND
	Scalar* kSumWI;
#endif

public:
	void addSingleFsi(LMPCons *localFsi);
	void constructKww();
	void constructKcw();
	void fetiBaseOpCoupled1(GenSolver<Scalar> *s, Scalar *localvec, const Scalar *interfvec,
	                        FSCommPattern<Scalar> *wiPat) const override;
	void fetiBaseOpCoupled2(const Scalar *uc, const Scalar *localvec, Scalar *interfvec,
	                        FSCommPattern<Scalar> *wiPat, const Scalar *fw = nullptr) const override;
	void multKbbCoupled(const Scalar *u, Scalar *Pu, Scalar *deltaF, bool errorFlag = true);
	void scaleAndSplitKww();
	void reScaleAndReSplitKww();
	void addSommer(SommerElement *ele); // XDEBUG

	// frequency sweep
	GenSparseMatrix<Scalar> *M;
	GenCuCSparse<Scalar> *Muc;
	GenSparseMatrix<Scalar> *C;
	GenCuCSparse<Scalar> *Cuc;
	GenSparseMatrix<Scalar> **C_deriv;
	GenSparseMatrix<Scalar> **Cuc_deriv;
	int numC_deriv;
	GenSparseMatrix<Scalar> **K_deriv;
	GenSparseMatrix<Scalar> **Kuc_deriv;
	int numK_deriv;
	int num_K_arubber;
	GenSparseMatrix<Scalar> **K_arubber_l;
	GenSparseMatrix<Scalar> **K_arubber_m;
	GenSparseMatrix<Scalar> **Kuc_arubber_l;
	GenSparseMatrix<Scalar> **Kuc_arubber_m;

private:
	// frequency sweep
	void makeFreqSweepLoad(Scalar *load, int iRHS, double omega);
	GenVector<Scalar> **a, **b;  // pade P, Q coefs
	int ia, ib;
	GenVector<Scalar> *P, *Q; // pade P(x), Q(x)
	bool rebuildPade;

public:
	void multM(Scalar *localrhs, GenStackVector<Scalar> **u, int k);
	void multMCoupled1(Scalar *localrhs, GenStackVector<Scalar> **u, int k,
	                   FSCommPattern<Scalar> *wiPat);
	void multMCoupled2(Scalar *localrhs, FSCommPattern<Scalar> *wiPat);
	void multWCAWE(Scalar *localrhs, GenStackVector<Scalar> **u,
	               Scalar *pU, Scalar *pb, int maxRHS, int iRHS);

	void pade(GenStackVector<Scalar> *sol,  GenStackVector<Scalar> **u, double *h, double x);
	void setRebuildPade(bool _rebuildPade) { rebuildPade = _rebuildPade; }

	// new B operators
	void multAddBrT(const Scalar *interfvec, Scalar *localvec, Scalar *uw = 0) const;
	void multBr(const Scalar *localvec, Scalar *interfvec, const Scalar *uc = 0, const Scalar *uw = 0) const;

	int *l2g;
	void makeLocalToGlobalDofMap();
	void multAddLT(const Scalar *localvec, Scalar *globalvec);
	void multAddLinv(const Scalar *localvec, Scalar *globalvec);
	void multLTinv(const Scalar *globalvec, Scalar *localvec);
	void multL(const Scalar *globalvec, Scalar *localvec);

};

typedef GenSubDomain<double> SubDomain;


#endif
