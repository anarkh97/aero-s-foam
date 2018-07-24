#ifndef _DEC_DOMAIN_H_
#define _DEC_DOMAIN_H_

#include <memory>
#include <Feti.d/DistrVector.h>
#include <Feti.d/DistrVectorSet.h>
#include <Paral.d/SubDOp.h>
#include <vector>
#include <Rom.d/EiGalerkinProjectionSolver.h>

extern Communicator * structCom;
extern GeoSource * geoSource;

class Domain;
template <class Scalar> class GenSubDomain;
template <class Scalar> class GenParallelSolver;
template <class Scalar> class GenFetiSolver;
template <class Scalar> class FetiBaseClass;
template <class Scalar> class GenDomainGroupTask;
template <class V> class SysState;
template <class Scalar> class DiagParallelSolver;
template <class Scalar> class GenSolver;
template <class Scalar> class GenMDDynamMat;
class DistrGeomState;
class MatrixTimers;
class Connectivity;
class Corotator;
class FetiSubCornerHandler;
class Rbm;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class MultiDomainRbm;

template<class Scalar>
class GenDecDomain 
{
 public:
	std::unique_ptr<Connectivity> mpcToSub_dual;
  Connectivity *mpcToMpc;
  Connectivity *mpcToCpu;
 protected:
  Domain *domain;
  MatrixTimers &mt;
  Connectivity *subToNode;
  Connectivity *subToSub;
  Connectivity *mpcToSub_primal;
  Connectivity *subToElem;
  const Connectivity *elemToNode;
  Connectivity *nodeToSub = nullptr;
  Connectivity *nodeToSub_copy = nullptr;
  Connectivity *cpuToSub;
  Connectivity *elemToSub;
  std::vector<GenSubDomain<Scalar> *> subDomain;
  DistrInfo *internalInfo, *internalInfo2;
  DistrInfo *masterSolVecInfo_;
  DistrInfo *nodeInfo;
  DistrInfo *nodeVecInfo, *eleVecInfo, *bcVecInfo;
  std::vector<DistrInfo*> vecInfoStore;
  Connectivity *grToSub;
  FILE *primalFile; // file to store primal residual

  // only used if user requests stress outputs
  GenDistrVector<Scalar> *stress;
  GenDistrVector<Scalar> *weight;
  Scalar *globalStress;
  Scalar *globalWeight;

    // TODO Switch to a std::map
	lc_sub_idx *glSubToLocal;            //<! \brief Local numbering of global subdomains.
	std::vector<gl_sub_idx> localSubToGl;            //<! \brief Local to global subdomain index map.
	int globalNumSub;             //<! \brief Global number of subdomains on all cpus.
	int numSub;                   //<! \brief Local number of subdomains on this process.

	FSCommunicator *communicator; //<! Communicator for all the structural processes.
  FSCommPattern<Scalar> *wiPat;

  int numCPU, myCPU;            // number of CPUs, my CPU number
  Connectivity *cpuToCPU;       // global problem connectivity
  //int *subToCPU;                // subdomain to cpu mapping  
  int numDualMpc, numPrimalMpc;

  int numWetInterfaceNodes;
  int *wetInterfaceNodes;
  int outFreqCount;
  int outEigCount;

  bool soweredInput;
  bool coarseLevel;

  int sizeSfemStress;
  bool firstOutput; 

  GenBasicAssembler<Scalar> *ba;
  GenBasicAssembler<Scalar> *ba2;

 public:
  GenDecDomain(Domain *d, Communicator * = structCom, bool = geoSource->binaryInput, bool = false);
  virtual ~GenDecDomain();

  virtual void clean();
  Domain *getDomain() { return domain; }

  auto getAllSubDomains() { return subDomain.data(); }
  GenSubDomain<Scalar>* getSubDomain(int isub) { return subDomain[isub]; }
  void setElemToNode(const Connectivity *_elemToNode) { elemToNode = _elemToNode; }
  void setSubToElem(Connectivity *_subToElem) { subToElem = _subToElem; }
  void setCPUMap(Connectivity *);
  Connectivity * getSubToSub() { return subToSub; }
  Connectivity * getElemToSub() { return elemToSub; }
  const Connectivity * getElemToNode() { return elemToNode; }
  Connectivity * getNodeToSub() { return nodeToSub; }
  Connectivity * getGroupToSub() { return grToSub; }
  int *getGlSubToLocal() { return glSubToLocal; }
  GenParallelSolver<Scalar> * getFetiSolver(GenDomainGroupTask<Scalar> &);
  void buildOps(GenMDDynamMat<Scalar>&, double, double, double, Rbm **rbm = 0, FullSquareMatrix **kelArray = 0,
				bool make_feti = true, FullSquareMatrix **melArray = 0, FullSquareMatrix **celArray = 0, bool factor = true);
  DiagParallelSolver<Scalar> *getDiagSolver(int nSub, GenSubDomain<Scalar> **, GenSolver<Scalar> **);
  void rebuildOps(GenMDDynamMat<Scalar>&, double, double, double, FullSquareMatrix** = 0, FullSquareMatrix** = 0, FullSquareMatrix** = 0);
  void subRebuildOps(int iSub, GenMDDynamMat<Scalar>&, double, double, double, FullSquareMatrix**, FullSquareMatrix**, FullSquareMatrix**);
  int getNumSub() { return numSub; }
  int getGlobalNumSub() { return globalNumSub; }
  int getNumCPU() { return numCPU; }
  Connectivity *getCpuToSub() const { return cpuToSub; }
  const Connectivity &getMpcToSub() const { return *mpcToSub_dual; }
  Connectivity *getMpcToSub_primal() { return mpcToSub_primal; }
  virtual void preProcess();
  virtual void postProcessing(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &f, double eigV = 0.0,
							  GenDistrVector<Scalar> *aeroF = 0, int x = 0, GenMDDynamMat<Scalar> *dynOps = 0,
							  SysState<GenDistrVector<Scalar> > *distState = 0, int ndflag = 0);
  virtual void postProcessing(DistrGeomState *u, GenDistrVector<Scalar> &, Corotator ***, double x = 0,
							  SysState<GenDistrVector<Scalar> > *distState = 0, GenDistrVector<Scalar> *aeroF = 0,
							  DistrGeomState *refState = 0, GenDistrVector<Scalar> *reactions = 0,
							  GenMDDynamMat<Scalar> *dynOps = 0, GenDistrVector<Scalar> *resF = 0);
  DistrInfo &solVecInfo() { return *internalInfo; } // unconstrained dofs
  const DistrInfo &masterSolVecInfo() const;
  DistrInfo &sysVecInfo() { return *internalInfo2; } // all dofs
  DistrInfo &ndVecInfo(); // all nodes
  // user defined control functions
  void extractControlData(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &,
						  GenDistrVector<Scalar> &,
						  Scalar *, Scalar *, Scalar *);
  void addUserForce(GenDistrVector<Scalar> &, Scalar*);
  void addCtrl(GenDistrVector<Scalar> &, Scalar*);
  // Non linear functions
  DistrInfo* elementVectorInfo();
  DistrInfo* pbcVectorInfo();
  void scaleDisp(GenDistrVector<Scalar> &u);
  void scaleInvDisp(GenDistrVector<Scalar> &u);
  void scaleDisp(GenDistrVector<Scalar> &u, double alpha);
  virtual void forceContinuity(GenDistrVector<Scalar> &u);
  virtual void forceAssemble(GenDistrVector<Scalar> &u);
  void setNewProperties(int);
  void assignRandMat();
  void retrieveElemset();
  virtual void setsizeSfemStress(int fileNumber);
  virtual int getsizeSfemStress() { return sizeSfemStress; }
  virtual Scalar* getSfemStress(int fileNumber) { return globalStress; }
  virtual void updateSfemStress(Scalar* str, int fileNumber);
  double computeStabilityTimeStep(GenMDDynamMat<Scalar>&);
  void extractSubDomainMPCs(int iSub);
  void reProcessMPCs();
	void setConstraintGap(DistrGeomState *geomState, DistrGeomState *refState, GenFetiSolver<Scalar> *fetisolver, double _lambda);
	void setConstraintGap(DistrGeomState *geomState, DistrGeomState *refState, FetiBaseClass<Scalar> *fetisolver, double _lambda);
  FSCommPattern<Scalar> * getWiCommPattern();
  GenAssembler<Scalar> * getSolVecAssembler();
  void exchangeInterfaceGeomState(DistrGeomState *geomState);
  void assembleNodalInertiaTensors(FullSquareMatrix **melArray);
  FSCommunicator * getCommunicator() { return communicator; }

 protected:
	/** \brief Build the basic subdomain data.
	 * \details Read the elements renumber their nodes with local numbers and form the unconstrained DSA.
	 * If there is FSI, build additional connectivities.
	 */
  void makeSubDomains();
  void renumberElements(int iSub);
  void createElemToNode();
  /** \brief Build the SComm for each subdomain, based on the input connectivities.
   *
   * @tparam ConnectivityType1
   * @tparam ConnectivityType2
   * @param nodeToSub
   * @param subToNode
   */
  template<class ConnectivityType1, class ConnectivityType2>
  void getSharedNodes(const ConnectivityType1 *nodeToSub, const ConnectivityType2 *subToNode) const;
  void addBMPCs();
  void makeSubToSubEtc();
  void preProcessBCsEtc();
  void preProcessMPCs();
  void distributeBCs();
  void distributeControlLawData();
  void distributeDiscreteMass();
  void renumberBC();
  void getSharedDOFs();
  void getSharedMPCs();
  void makeCorners();
  void makeInternalInfo();
  void makeSolVecInfo();
  void makeSysVecInfo();
  void makeNodeInfo();
  void setNonTrivialMasterFlag(DistrInfo &);
  void getCPUMap();
  void makeSubDMaps();
  void constructSubDomains(int iSub);
  void makeMpcToMpc();
  void makeGlobalMpcToMpc(Connectivity *procMpcToMpc);
  void makeMpcToSub();
  void buildFFP(GenDistrVector<Scalar> &u, FILE *fffp, bool direction);
  void makeCornerHandler(int iSub, FetiSubCornerHandler **cornerHandler);
  void setLocalCorners(int iSub, FetiSubCornerHandler **cornerHandler);
  void deleteMPCs();
  void extractPosition(int iSub, DistrGeomState &geomState, GenDistrVector<Scalar> &x);
  virtual void setMpcRhs(int iSub, GenDistrVector<Scalar> &cu, double t, int flag);
  void dispatchInterfaceGeomState(int isub, FSCommPattern<double> *geomStatePat, DistrGeomState *geomState);
  void collectInterfaceGeomState(int isub, FSCommPattern<double> *geomStatePat, DistrGeomState *geomState);
  void dispatchInterfaceNodalInertiaTensors(int isub, FSCommPattern<double> *pat, FullSquareMatrix **melArray);
  void collectInterfaceNodalInertiaTensors(int isub, FSCommPattern<double> *pat);

 public:
  void printLMPC();
  void makeBlockCyclicDistrInfo(DistrInfo &, int globalLen, int blockSize);
  void makeNonOverlappingDistrInfo(DistrInfo &info);
  void makeBasicDistrInfo(DistrInfo &, int(Domain::*)() const);

 private:
  void initialize();
  // output functions
  void getElementForce(GenDistrVector<Scalar>&, int, int, double);
 public:
  virtual void getStressStrain(GenDistrVector<Scalar>&, int, int, double, int printFlag=0);
  void getEnergies(GenDistrVector<Scalar> &disp, GenDistrVector<Scalar> &extF, int fileNumber, double time,
				   SysState<GenDistrVector<Scalar> > *distState, GenMDDynamMat<Scalar> *dynOps,
				   GenDistrVector<Scalar> *aeroF);
  void getEnergies_b(DistrGeomState *geomState, GenDistrVector<Scalar> &extF, Corotator ***allCorot, int fileNumber,
				   double time, SysState<GenDistrVector<Scalar> > *distState, GenMDDynamMat<Scalar> *dynOps,
				   GenDistrVector<Scalar> *aeroF);
  void getDissipatedEnergy(DistrGeomState *geomState, Corotator ***allCorot, int fileNumber, double time);
  MultiDomainRbm<Scalar>* constructRbm(bool printFlag = true);
 private:
  void getElementStressStrain(GenDistrVector<Scalar>&, int, int, double, int printFlag=0);
  void getElementAttr(int, int, double);
  void getPrincipalStress(GenDistrVector<Scalar>&, int, int, double);
  void getElementPrincipalStress(GenDistrVector<Scalar> &u, int, int, double);
  void getStressStrain(DistrGeomState *u, Corotator ***, int, int, double, DistrGeomState *refState);
  void getPrincipalStress(DistrGeomState *u, Corotator ***, int, int, double, DistrGeomState *refState);
  void getElementPrincipalStress(DistrGeomState *u, Corotator ***, int, int, double, DistrGeomState *refState);
  void computeSubdElemForce(int iSub, Scalar *globForce,
							GenDistrVector<Scalar> *u, int fileNumber, int Findex);
  void computeSubdStress(int, GenDistrVector<Scalar>*, GenDistrVector<Scalar>*,
						 GenDistrVector<Scalar>*, int, int) const;
  void computeSubdElemStress(int, Scalar *, GenDistrVector<Scalar> *, int, int) const;
  void computeSubdStress_NL(int iSub, GenDistrVector<Scalar> *globStress,
						 GenDistrVector<Scalar> *globWeight, DistrGeomState *u,
						 Corotator ***allCorot, int *, int *Findex, DistrGeomState *refState) const;
  void getElementStressStrain(DistrGeomState *gs, Corotator ***allCorot,
							  int fileNumber, int Findex, double time, DistrGeomState *refState) const;
  void computeSubdElemStress_NL(int iSub, Scalar *glElemStress,
							 DistrGeomState *u, Corotator ***allCorot,
							 int fileNumber, int Findex, DistrGeomState *refState) const;
  void outputPrimal(GenDistrVector<Scalar>& primal, int iter);
  void getPrimalVector(int fileNumber, Scalar (*xyz)[11], int numNodes,
					   int ndof, double time);//DofSet::max_known_nonL_dof
  void getPrimalScalar(int fileNumber, Scalar (*xyz)[11], int numNodes,
					   int dof, double time);//DofSet::max_known_nonL_dof
  void getAeroForceScalar(int fileNumber, Scalar (*mergedAeroF)[6],
						  int numNodes, int dof, double time);
  void scaleSubDisp(int iSub, GenDistrVector<Scalar> &u);
  void scaleInvSubDisp(int iSub, GenDistrVector<Scalar> &u);
  void scaleSubDisp_2(int iSub, GenDistrVector<Scalar> &u, double alpha) const;
  void subGetEnergies(int iSub, GenDistrVector<Scalar> &disp, GenDistrVector<Scalar> &extF, double time,
					  SysState<GenDistrVector<Scalar> > *distState, GenMDDynamMat<Scalar> *dynOps,
					  GenDistrVector<Scalar> *aeroF, double *subW);
  void subGetEnergies_b(int iSub, DistrGeomState *geomState, GenDistrVector<Scalar> &extF,
					  Corotator ***allCorot, double time, SysState<GenDistrVector<Scalar> > *distState,
					  GenMDDynamMat<Scalar> *dynOps, GenDistrVector<Scalar> *aeroF, double *subW);
  void subGetDissipatedEnergy(int iSub, DistrGeomState *geomState, Corotator ***allCorot, double *subD);

  // Helmholtz Fluid functions
  void distribBC(int iSub, gsl::span<GenSubDomain<Scalar> *> sd, Domain *domain,
                 int *somToSub, int *scaToSub, int *neumToSub, int (*wetToSub)[2],
                 int *sBoundFlag);
  void buildLocalFFP(int iSub, GenDistrVector<Scalar> *u,
					 Scalar **ffp, int *numSample, double (*dir)[3], bool direction);
  void getWError(int iSub, GenDistrVector<Scalar> *u,
				 double *l2err, double *h1err, double *l2, double *h1);
  // coupled_dph functions
  void preProcessFSIs();
  void distributeWetInterfaceNodes();
  void setSubWetInterface(int iSub, int *nWetInterfaceNodesPerSub, int **subWetInterfaceNodes);
  void getSharedFSIs();
  void markSubWetInterface(int iSub, int *nWetInterfaceNodesPerSub, int **subWetInterfaceNodes);
  void addFsiElements();
  
  GenBasicAssembler<Scalar> * solVecAssemblerNew();
};

#endif
