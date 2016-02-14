#ifndef _SUB_DOMAIN_H_
#define _SUB_DOMAIN_H_

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/ControlLawInfo.h>
#include <Feti.d/DistrVector.h>
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

class BaseSub : virtual public Domain 
//class BaseSub : public Domain
{
 protected:
  SComm *scomm;
  int subNumber;
  int localSubNumber; // relevant when running in distributed
// RT: 030813
//  int *glToLocalNode;
//  int *glToLocalElem;
  GlobalToLocalMap glToLocalNode;
  GlobalToLocalMap glToLocalElem;
  int *glNums;
  int *glElems;
  int glNumNodes;
  int *weight; // DOF weights (i.e. number of subd sharing that dof)
  int *weightPlus;
  double *bcx;
  DComplex *bcxC; // FETI-H
  double *vcx, *acx;
  int *locToGlSensorMap;
  int *locToGlActuatorMap;
  int *locToGlUserDispMap;
  int *locToGlUserForceMap;
  int boundLen;
  int *boundMap;
  int *dualToBoundary;
  int internalLen;
  int *internalMap;
  int crnDofSize;
  int *crnPerNeighb;
  long memK;       // memory necessary to store K(s)
  long memPrec;    // memory necessary to store Preconditioner
#ifdef DISTRIBUTED
  // for distributed output of single nodes
  int numNodalOutput;
  int *outputNodes;
  int *outIndex;
#endif
  int globalNMax;  // highest global node number of all the nodes in this subdomain
  int globalEMax;  // highest global element number of all the elements in this subdomain
  int totalInterfSize;
  int *allBoundDofs;
  Rbm *rigidBodyModes;
  Rbm *rigidBodyModesG;

 public:
  BaseSub();
  BaseSub(Domain &dom, int sn, Connectivity &con, Connectivity &nds, int gn);
  BaseSub(Domain &dom, int sn, int nNodes, int *nds,
          int nElems, int *elems, int gn);
  BaseSub(Domain &dom, int sn, CoordSet* nodes, Elemset* elems, int *glNodeNums,
          int *glElemNums, int gn); // PJSA: for new sower
  virtual ~BaseSub();

  // Multiple Point Constraint (MPC) Data
  int numMPC;             // number of local Multi-Point Constraints
  int *localToGlobalMPC;  // local to global MPC numbering
  GlobalToLocalMap globalToLocalMPC; // alternative data structure for global to local MPC numbering
                                     // not a pointer so don't have to de-reference before using [] operator

  int numMPC_primal;
  int *localToGlobalMPC_primal;
  GlobalToLocalMap globalToLocalMPC_primal;

  int *cornerMap;
  int *cornerEqNums; // unique equation numbers for subdomain corner dofs
  int *localCornerList;
  DofSet *cornerDofs;
  int *cornerNodes;   // corner node in local numbering
  bool *isCornerNode;   // true for node which is a corner node; false otherwise
  int *glCornerNodes; // corner nodes in global numbering
  int *glCrnGroup;    // group of each corner node (global numbering)
  int nGrbm;
  int numCRN;
  int numCRNdof;
  ConstrainedDSA *cc_dsa;
  int *ccToC; // from cc_dsa to c_dsa
  int *cToCC; // from c_dsa to cc_dsa
  DofSet **boundaryDOFs;	
  int nCDofs;
  int *neighbNumGRBMs;
  int *edgeDofSize;       // number of edge dof per neighbor
  int *edgeDofSizeTmp;  // XXXX
  double k_f, k_p, k_s, k_s2;  // wave numbers for FETI-DPH for this subdomain
  double *neighbK_p, *neighbK_s, *neighbK_s2, *neighbK_f;  // neighbors' wave numbers
  double Ymod, Prat, Dens, Thih, Sspe;  // Young's modulus, Poisson ration, density, thickness, speed of sound
  double *neighbYmod, *neighbPrat, *neighbDens, *neighbThih, *neighbSspe;  // neighbor's values

  int dofWeight(int i) { return weight[i]; }
  int crnDofLen()            { return crnDofSize; }
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
  int getNumNodalOutput()     { return numNodalOutput; }
#endif
  int getBC(BCond *, int, int *, BCond *&);
  void setGlNodes(int *globalNodeNums) { glNums = globalNodeNums; }
  int *getGlNodes()           { return glNums; }
  int *getGlElems()           { return glElems; }
  int *getGlMPCs()            { return localToGlobalMPC; }
  int glToPackElem(int e)     { return (geoSource->glToPackElem(e) > globalEMax) ? -1 : glToLocalElem[geoSource->glToPackElem(e)]; }
  int *getSensorDataMap()     { return locToGlSensorMap; }
  int *getActuatorDataMap()   { return locToGlActuatorMap; }
  int *getUserDispDataMap()   { return locToGlUserDispMap; }
  int *getUserForceDataMap()  { return locToGlUserForceMap; }
  int countElemNodes();
  int numMPCs()               { return numMPC; }
  int numMPCs_primal()        { return numMPC_primal; }
  int globalNumNodes(); 
  int numNodes()              { return numnodes; }
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
  int* makeBMaps(DofSetArray *dofsetarray=0);
  int* makeIMaps(DofSetArray *dofsetarray=0);
  int subNum()                { return subNumber; }
  int localSubNum()           { return localSubNumber; }
  int localLen()              { return (cc_dsa) ? cc_dsa->size() : c_dsa->size(); }
  ConstrainedDSA * getCCDSA()  { return (cc_dsa) ? cc_dsa : c_dsa; }
  int localRLen()             { return cc_dsa->size(); }
  void sendNumNeighbGrbm(FSCommPattern<int> *pat);
  void recvNumNeighbGrbm(FSCommPattern<int> *pat);
  void deleteLocalRBMs() { if(rigidBodyModes) delete rigidBodyModesG; rigidBodyModesG = 0; }
  void putNumMPC(int *ptr) { ptr[subNumber] = numMPC; }
  void putLocalToGlobalMPC(int *ptr, int *tg) { for(int i=0; i<numMPC; ++i) tg[ptr[subNumber]+i] = localToGlobalMPC[i]; }
  void putNumMPC_primal(int *ptr) { ptr[subNumber] = numMPC_primal; }
  void putLocalToGlobalMPC_primal(int *ptr, int *tg) { for(int i=0; i<numMPC_primal; ++i) tg[ptr[subNumber]+i] = localToGlobalMPC_primal[i]; }

  // Multiple Point Constraint (MPC) functions
  int getNumMpc() const       { return numMPC; }
  int zColDim() { return rigidBodyModesG->Zstar->numCol(); }
  int zRowDim() { return rigidBodyModesG->Zstar->numRow(); }
  void addMPCsToGlobalZstar(FullM *globalZstar, int startRow, int startCol, int numCol);
  void makeLocalToGroupMPC(Connectivity *groupToMPC);
  void findEdgeNeighbors();
  void makeMpcInterface(Connectivity *subToMpc, Connectivity *lmpcToSub,
                        Connectivity *subToSub_mpc);
  void makeFsiInterface(Connectivity *subToFsi, Connectivity *fsiToSub,
                        Connectivity *subToSub_fsi);

  bool checkForColinearCrossPoints(int numCornerPoints, int *localCornerPoints);
  void addCornerPoints(int *glCornerList);
  int *getLocalCornerNodes()  { return cornerNodes; }
  int numCorners()	      { return numCRN; }
  int *getCornerNodes()       { return glCornerNodes; }
  int numCornerDofs()	      { return numCRNdof; }
  int numCoarseDofs();
  int nCoarseDofs()           { return nCDofs; }
  int numEdgeDofs(int i)      { return edgeDofSize[i]; }
  
  // variables and routines for parallel GRBM algorithm and floating bodies projection
  // and MPCs (rixen method)
 public:
  int group;
 protected:
  int numGroupRBM, groupRBMoffset;
  int *neighbGroup;
  int *neighbNumGroupGrbm;
  int *neighbGroupGrbmOffset;
  int numGlobalRBMs;

 protected:
  Connectivity *localMpcToMpc;
  Connectivity *localMpcToGlobalMpc;
  bool *faceIsSafe;
  int *localToGroupMPC, *localToBlockMPC;
  int *boundDofFlag;  // boundDofFlag[i] = 0 -> perfect interface dof  (not contact or mpc)
                      // boundDofFlag[i] = 1 -> node-to-node contact interface dof
                      // boundDofFlag[i] = 2 -> mpc dof, only used for rixen method, domain->mpcflag = 1
                      // note: boundDofFlag[i] > 2 can be used for future extensions, eg mortar contact
  bool *masterFlag; // masterFlag[i] = true if this sub is the "master" of allBoundDofs[i]
  bool *internalMasterFlag;
  int masterFlagCount;
  int factorial(int n); // computes n!
  bool *mpcMaster;  // mpcMaster[i] = true if this subd contains masterdof for mpc i
  Connectivity *mpcToDof;
  Connectivity *localMpcToBlock;
  Connectivity *blockToLocalMpc;
  Connectivity *blockToBlockMpc;
  Connectivity *localMpcToBlockMpc;
  Connectivity *mpcToBoundDof;
  double *localLambda;  // used for contact pressure output
  int *invBoundMap;
  int *mpclast;

 public:
  void makeLocalMpcToGlobalMpc(Connectivity *mpcToMpc);
  void setLocalMpcToBlock(Connectivity *mpcToBlock, Connectivity *blockToMpc);
  void setGlCrnGroup(int *_gcg) { glCrnGroup = _gcg; }
  void setGroup(Connectivity *subToGroup) { group = (*subToGroup)[subNumber][0]; }
  void setNumGroupRBM(int *ngrbmGr); 
  void getNumGroupRBM(int *ngrbmGr);
  void addNodeXYZ(double *centroid, double* nNodes);
  void sendNeighbGrbmInfo(FSCommPattern<int> *pat);
  void receiveNeighbGrbmInfo(FSCommPattern<int> *pat);
  void setCommSize(FSCommPattern<int> *pat, int size);
  void setCommSize(FSCommPattern<double> *, int size);
  void setMpcNeighbCommSize(FSCommPattern<int> *pt, int size);
  int numSPCs() { return c_dsa->getInvRCNmax(); }
  void addSPCsToGlobalZstar(FullM *globalZstar, int &zRow, int zColOffset);
  void initializeFaceSafety();
  void locateUnsafeFaces();
  void sendFaceSafetyInfo(FSCommPattern<int> *sPat);
  void receiveFaceSafetyInfo(FSCommPattern<int> *sPat);
  void printUnsafeFaces();

 protected:
  int bodyRBMoffset;
  void computeInterfaceCardinality(int *nodeweight);

 public:
  void setCorners(int nCorners, int *crnList);
  SComm *getSComm() { return scomm; }
  void setSComm(SComm *sc);
  void setSendData(int neighb, void *data)
     { scomm->setExchangeData(neighb, data); }
  void *getExchangeData(int neighbID) // Communication mechanism
     { return scomm->getExchangeData(neighbID); }
  void *getExchangePointer(int neighbID) // Communication mechanism
     { return scomm->exchangeData[neighbID]; }
  int numNeighbors() { return scomm->numNeighb;}
  int numEdgeNeighbors() { return scomm->numEdgeNeighb; }
  bool isEdgeNeighbor(int neighb) { return scomm->isEdgeNeighb[neighb]; }
  void setPairsNbTotal(int pairsNbTotal);
  int interfLen(); // Total length for the local interface
  int halfInterfLen(); // Length of the "half interface"
  void computeMasterFlag(Connectivity *mpcToSub);
  bool* getMasterFlag() { return masterFlag; }
  const bool* getInternalMasterFlag();
 
 protected:
  void computeInternalMasterFlag();

 public:
  void setNodeCommSize(FSCommPattern<int> *, int d = 1);
  void setNodeCommSize(FSCommPattern<double> *, int d = 1);
  void setNodeCommSize(FSCommPattern<DComplex> *, int d = 1);
  void setDofCommSize(FSCommPattern<int> *);
  void setDofCommSize(FSCommPattern<double> *);
  void setDofCommSize(FSCommPattern<DComplex> *);
  void setDofPlusCommSize(FSCommPattern<double> *);
  void setDofPlusCommSize(FSCommPattern<DComplex> *);
  void setRbmCommSize(int numRBM, FSCommPattern<double> *);
  void setRbmCommSize(int numRBM, FSCommPattern<DComplex> *);

  // for timing file
  double getSharedDofCount();
  int getTotalDofCount();
  void setBodyRBMoffset(int _boff) { bodyRBMoffset = _boff; }
  SubCornerHandler *getCornerHandler();

  void initialize();
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
  void GramSchmidt(double *Q, bool *isUsed, int numdofperNode, int nQPerNeighb);
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
//  void updateKappa() { kappa = domain->getWaveNumber(); }
//  void updateKappaAndScalings() { kappa = domain->getWaveNumber();
//                       coupledScaling = domain->coupledScaling;
//                       cscale_factor  = domain->cscale_factor;
//                       cscale_factor2 = domain->cscale_factor2;
//  }

  // coupled_dph
 protected:
  bool *wetInterfaceMark;
  bool *wetInterfaceFluidMark;
  bool *wetInterfaceStructureMark;
  bool *wetInterfaceCornerMark;
  int numWIdof;  // number of dofs on the wet interface (both fluid and structure)
  int numWInodes;  // number of nodes on the wet interface (both fluid and structure)
  int *wetInterfaceMap;  // dof map
  int *wetInterfaceNodeMap; 
  int *wetInterfaceNodes;
  int *wDofToNode; //HB
  DofSet *wetInterfaceDofs;
  int *numNeighbWIdof;
  Connectivity *drySharedNodes;
  bool *wiMaster;
  GlobalToLocalMap glToLocalWImap;
  GlobalToLocalMap *neighbGlToLocalWImap;
  int numFsiNeighb;
  int *fsiNeighb;
  bool haveAverageMatProps;
  int *wiInternalMap;
  bool isMixedSub;
  int edgeQindex[2];
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
  void setWIoneCommSize(FSCommPattern<int> *pat);
  void sendNumWIdof(FSCommPattern<int> *sPat);
  void recvNumWIdof(FSCommPattern<int> *sPat);
  void setWImapCommSize(FSCommPattern<int> *pat);
  void sendWImap(FSCommPattern<int> *pat);
  void recvWImap(FSCommPattern<int> *pat);
  void makeDSA();
  void makeCDSA();
  void makeCCDSA();
  int numWetInterfaceDofs() { return numWIdof; }
  GlobalToLocalMap& getGlToLocalWImap() { return glToLocalWImap; }
  GlobalToLocalMap& getNeighbGlToLocalWImap(int i) { return neighbGlToLocalWImap[i]; }
  void zeroEdgeDofSize();
  void mergeInterfaces();
  void markCornerDofs(int *glCornerDofs);

#ifdef HB_COUPLED_PRECOND
  Connectivity* precNodeToNode;
#endif
};

template<class Scalar>
class GenSubDomain : public BaseSub 
{
 private:
  Scalar *kweight; // stiffness weights (i.e. sum of Kii for all subd sharing that dof)
  Scalar *deltaFmpc;
  int *cornerWeight;
  void applyBtransposeAndScaling(Scalar *u, Scalar *v, Scalar *deltaU = 0, Scalar *localw = 0);
  void applyScalingAndB(Scalar *res, Scalar *Pu, Scalar *localw = 0);
  void initialize();

 protected:
  Scalar *scaling;
  void sendDOFList(FSCommPattern<int> *pat); // Send to neighbors the list of DOFs on the shared nodes
  GenSkyMatrix<Scalar> * makeSkyK(Connectivity &nton, Scalar trbm);

 public:
  GenSparseSet<Scalar>      *Src;
  GenSparseSet<Scalar>      *Qrc;
  GenSolver<Scalar>         *Krr;
  GenSparseMatrix<Scalar>   *KrrSparse;
  Scalar                    **BKrrKrc;
  GenAssembledFullM<Scalar> *Kcc;
  GenCuCSparse<Scalar>      *Krc;
  GenCuCSparse<Scalar>      *Grc;
  Scalar                    *rbms;
  Scalar                    *interfaceRBMs;
  GenFullM<Scalar>          *qtkq;
  GenSparseMatrix<Scalar>   *KiiSparse;
  GenSolver<Scalar>         *KiiSolver;
  GenCuCSparse<Scalar>      *Kib;
  GenSparseMatrix<Scalar>   *MPCsparse;
  GenDBSparseMatrix<Scalar> *Kbb;    // for preconditioning
  Corotator           	    **corotators;
  Scalar 		    *fcstar;
  Scalar                    *QtKpBt;
  Scalar                    *locKpQ;

  int *glBoundMap;
  int *glInternalMap;
 
  SubLMPCons<Scalar> **mpc; // multiple point constraints
  SubLMPCons<Scalar> **mpc_primal;

 private:
  GenSolver<Scalar> *localCCtsolver;
  Scalar *diagCCt;
  int lengthCCtData;
  int *CCtrow, *CCtcol;
  Scalar* CCtval;
  Scalar *bcx_scalar;
  int *mpcStatus;
  bool *mpcStatus1, *mpcStatus2;

  // templated RBMs
  GenFullM<Scalar> Rstar;
  GenFullM<Scalar> Rstar_g;
  GenFullM<Scalar> *sharedRstar_g;
  GenFullM<Scalar> *tmpRstar_g;
  GenFullM<Scalar> **G;
  GenFullM<Scalar> **neighbG;

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
  void fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec);
  void fetiBaseOp(Scalar *uc,GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec);
  void fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec, Scalar *beta);
  void interfaceJump(Scalar *iterfData, FSCommPattern<Scalar> *vPat);
  void sendInterf(Scalar *interfvec, FSCommPattern<Scalar> *vPat);
  void extractAndSendInterf(Scalar *subvec, FSCommPattern<Scalar> *pat);
  void assembleInterf(Scalar *subvec, FSCommPattern<Scalar> *pat);
  void splitInterf(Scalar *subvec);
  void assembleInterfInvert(Scalar *subvec, FSCommPattern<Scalar> *pat);
  void getHalfInterf(Scalar *s, Scalar *t);
  void getHalfInterf(Scalar *s, Scalar *t, Scalar *ss, Scalar *tt);
  void scatterHalfInterf(Scalar *s, Scalar *loc);
  void rebuildInterf(Scalar *v, FSCommPattern<Scalar> *vPat);
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

  void expandRBM(Scalar *localR, VectorSet &globalR);
  void getSRMult(Scalar *lvec, Scalar *lbvec, int nRBM, double *locRBMs, Scalar *alpha);
  void sendInterfaceGrbm(FSCommPattern<Scalar> *rbmPat);
  void receiveInterfaceGrbm(FSCommPattern<Scalar> *rbmPat);
  void sendDeltaF(Scalar *deltaF, FSCommPattern<Scalar> *vPat);
  double collectAndDotDeltaF(Scalar *deltaF, FSCommPattern<Scalar> *vPat);
  void makeKbbMpc();
  void rebuildKbb();
  void makeKbb(DofSetArray *dofsetarray=0);
  void factorKii();
  void factorKrr();
  void multKbb(Scalar *, Scalar *, Scalar * = 0, Scalar * = 0, bool errorFlag = true);
  void multDiagKbb(Scalar *, Scalar *);
  void multFi(GenSolver<Scalar> *s, Scalar *, Scalar *);
  void multMFi(GenSolver<Scalar> *s, Scalar *, Scalar *, int numRHS);
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
  void sendDiag(GenSparseMatrix<Scalar> *s, FSCommPattern<Scalar> *vPat);
  void collectScaling(FSCommPattern<Scalar> *vPat);
  void fSend(Scalar *locF, FSCommPattern<Scalar> *vPat, Scalar *locFw = 0);
  void fScale(Scalar *locF, FSCommPattern<Scalar> *vPat, Scalar *locFw = 0);
  void fSplit(Scalar *locF);
  void updatePrescribedDisp(GeomState *geomState, Scalar deltaLambda);
  void updatePrescribedDisp(GeomState *geomState);
  Scalar displacementNorm(Scalar *displacement);
  void firstAssemble(GenSparseMatrix<Scalar> *K);
  void clearTemporaries() { delete [] glToLocalNode; glToLocalNode = 0; }
  void initMpcScaling();
  void makeZstarAndR(double *centroid);  // makes Zstar and R
  void makeKccDofs(DofSetArray *cornerEqs, int augOffset, Connectivity *subToEdge, int mpcOffset = 0);
  void assembleKccStar(GenSparseMatrix<Scalar> *KccStar);
  void deleteKcc();
  void multKbbMpc(Scalar *u, Scalar *Pu, Scalar *deltaU, Scalar *deltaF, bool errorFlag = true);
  void normalizeCstep1(Scalar *cnorm);
  void normalizeCstep2(Scalar *cnorm);
  void getQtKQ(GenSolver<Scalar> *s);
  void getQtKQ(int iMPC, Scalar *QtKQ);
  void multQt(int glMPCnum, Scalar *V, int numV, Scalar *QtV);
  void multQt(int glMPCnum, Scalar *x, Scalar *result);
  void multQtKBt(int glNumMPC, Scalar *G, Scalar *QtKBtG, Scalar alpha=1.0, Scalar beta=1.0);
  void gatherDOFList(FSCommPattern<int> *pat);
  void gatherDOFListPlus(FSCommPattern<int> *pat);

  friend class GenDistrDomain<Scalar>;
  friend class GenDecDomain<Scalar>;
  friend class GenFetiOp<Scalar>;
  friend class GenFetiSolver<Scalar>;

  GenAssembledFullM<Scalar> *getKcc() { return Kcc; }
  void makeQ();
  void precondGrbm();
  void orthoWithOrigR(Scalar *origV, Scalar *V, int numV, int length);
  void ortho(Scalar *vectors, int numVectors, int length);
  DofSet* getCornerDofs() { return cornerDofs; }
  void setMpcSparseMatrix();
  void assembleMpcIntoKcc();
  void multKcc();
  void reMultKcc();
  void multKrc(Scalar *fr, Scalar *uc);
  void multfc(Scalar *fr, Scalar *bf);
  void multFcB(Scalar *bf);
  Scalar *getfc() { return fcstar; }
  void getFr(Scalar *f, Scalar *fr);
  void getFc(Scalar *f, Scalar *fc);
  void getFw(Scalar *f, Scalar *fw);
  void mergeUr(Scalar *ur, Scalar *uc, Scalar *u, Scalar *lambda = 0);
  int numRBM() { return nGrbm; }
  void makeEdgeVectorsPlus(bool isFluidSub = false);
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
  Scalar getMpcRhs(int iMPC);
  Scalar getMpcRhs_primal(int iMPC);
  void constraintProduct(int num_vect, const double* R[], Scalar** V, int trans);
  void addConstraintForces(std::map<std::pair<int,int>, double> &mu, std::vector<double> &lambda, GenVector<Scalar> &f);
  void addCConstraintForces(std::map<std::pair<int,int>, double> &mu, std::vector<double> &lambda, GenVector<Scalar> &fc, double s);
  void locateMpcDofs();
  void makeLocalMpcToDof(); //HB: create the LocalMpcToDof connectivity for a given DofSetArray 
  void makeLocalMpcToMpc();
  void deleteMPCs();

  void projectActiveIneq(Scalar *v);
  void split(Scalar *v, Scalar *v_f, Scalar *v_c);
  void bmpcQualify(std::vector<LMPCons *> *bmpcs, int *pstatus, int *nstatus);
  void updateActiveSet(Scalar *v, double tol, int flag, bool &statusChange);
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
  void setMpcDiagCommSize(FSCommPattern<Scalar> *mpcDiagPat);
  void sendMpcDiag(FSCommPattern<Scalar> *mpcDiagPat);
  void collectMpcDiag(FSCommPattern<Scalar> *mpcDiagPat);

  void extractMpcResidual(Scalar *subv, GenVector<Scalar> &mpcv, SimpleNumberer *mpcEqNums);
  void insertMpcResidual(Scalar *subv, GenVector<Scalar> &mpcv, SimpleNumberer *mpcEqNums);
  void solveLocalCCt(Scalar *subv);
  void extractBlockMpcResidual(int block, Scalar *subv, GenVector<Scalar> *mpcv, 
                               SimpleNumberer *blockMpcEqNums);
  void insertBlockMpcResidual(Scalar *subv, GenVector<Scalar> **mpcv, Connectivity *mpcToBlock,
                              SimpleNumberer **blockMpcEqNums);
  void setMpcCommSize(FSCommPattern<Scalar> *mpcPat);
  void sendMpcInterfaceVec(FSCommPattern<Scalar> *mpcPat, Scalar *interfvec);
  void combineMpcInterfaceVec(FSCommPattern<Scalar> *mpcPat, Scalar *interfvec);
  void sendMpcScaling(FSCommPattern<Scalar> *mpcPat);
  void collectMpcScaling(FSCommPattern<Scalar> *mpcPat);
  void setMpcCommSize(FSCommPattern<int> *mpcPat);
  void sendMpcStatus(FSCommPattern<int> *mpcPat, int flag);
  void recvMpcStatus(FSCommPattern<int> *mpcPat, int flag, bool &statusChange);
  void printMpcStatus();
  void initMpcStatus();
  void saveMpcStatus();
  void restoreMpcStatus();
  void saveMpcStatus1();
  void saveMpcStatus2();
  void cleanMpcData();
  void subtractMpcRhs(Scalar *interfvec);
  void setLocalLambda(Scalar *localLambda);
  void computeContactPressure(Scalar *globStress, Scalar *globWeight);
  void computeLocalContactPressure(Scalar *stress, Scalar *weight);
  void getLocalMpcForces(double *mpcLambda, DofSetArray *cornerEqs,
                         int mpcOffset, GenVector<Scalar> &uc);
  void getConstraintMultipliers(std::map<std::pair<int,int>,double> &mu, std::vector<double> &lambda);
  void getLocalMultipliers(std::vector<double> &lambda);
  void setMpcRhs(Scalar *interfvec, double t, int flag);
  void updateMpcRhs(Scalar *interfvec);
  double getMpcError();
  bool* getMpcMaster() { return mpcMaster; }

 protected:
  double *mpcForces;

  // coupled_dph
  GenDBSparseMatrix<Scalar> *Kww;
  GenCuCSparse<Scalar>      *Kcw;
  GenMpcSparse<Scalar> *Kcw_mpc;
  GenCuCSparse<Scalar> *Krw;
  GenFsiSparse<Scalar> *neighbKww;   
  Scalar *localw;
  Scalar *localw_copy;
  Scalar Bcx(int i);
  void makeBcx_scalar();
  Scalar *deltaFwi;
  bool isWetInterfaceNode(int n) { return (wetInterfaceNodeMap) ? (wetInterfaceNodeMap[n] > -1) : false; }
  bool isWetInterfaceDof(int d) { return (wetInterfaceMap) ? (wetInterfaceMap[d] > -1) : false; }
  Scalar *wweight;
#ifdef HB_COUPLED_PRECOND
  Scalar* kSumWI;
#endif
  
 public:
  void addSingleFsi(LMPCons *localFsi); 
  void constructKrw();
  void constructKww();
  void constructKcw();
  void setWICommSize(FSCommPattern<Scalar> *wiPat);
  void setCSCommSize(FSCommPattern<Scalar> *csPat);
  void fetiBaseOpCoupled1(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec, 
                          FSCommPattern<Scalar> *wiPat);
  void fetiBaseOpCoupled2(Scalar *uc, Scalar *localvec, Scalar *interfvec, 
                          FSCommPattern<Scalar> *wiPat, Scalar *fw = 0);
  void multKbbCoupled(Scalar *u, Scalar *Pu, Scalar *deltaF, bool errorFlag = true);
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
  void pade(GenStackVector<Scalar> *sol,  GenStackVector<Scalar> **u, double *h, double x);
  void setRebuildPade(bool _rebuildPade) { rebuildPade = _rebuildPade; }

  // new B operators
  void multAddBrT(Scalar *interfvec, Scalar *localvec, Scalar *uw = 0);
  void multBr(Scalar *localvec, Scalar *interfvec, Scalar *uc = 0, Scalar *uw = 0);
  void multAddCT(Scalar *interfvec, Scalar *localvec);
  void multC(Scalar *localvec, Scalar *interfvec);

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

  // R matrix construction and access
  void makeLocalRstar(FullM **Qtranspose); // this is used by decomposed domain GRBM algorithm
  void useKrrNullspace();
  // R matrix-vector multiplication
  void addRalpha(Scalar *u, GenVector<Scalar> &alpha);  // u += R_g*alpha
  void addTrbmRalpha(Scalar *rbms, int nrbms, int glNumCDofs, Scalar *alpha, Scalar *ur); // u += R_g*alpha
  void assembleE(GenVector<Scalar> &e, Scalar *f); // e = R^T*f
  void assembleTrbmE(Scalar *rbms, int nrbms, int glNumCDofs, Scalar *e, Scalar *fr); // e = R^T*f

  // G matrix construction and destruction
  void makeG();
  void makeTrbmG(Scalar *rbms, int nrbms, int glNumCDofs);
  void assembleGlobalG(GenFullM<Scalar> *globalG);
  void setGCommSize(FSCommPattern<Scalar> *pat);
  void sendG(FSCommPattern<Scalar> *rbmPat);
  void receiveG(FSCommPattern<Scalar> *rbmPat);
  void zeroG();
  void deleteG();
  // G matrix-vector multiplication
  void multG(GenVector<Scalar> &x, Scalar *y, Scalar alpha);  // y = alpha*G*x
  void trMultG(Scalar *x, GenVector<Scalar> &y, Scalar alpha); // y = alpha*G^T*x
  // (G^T*G) matrix assembly
  void assembleGtGsolver(GenSparseMatrix<Scalar> *GtGsolver);

  // R_g matrix construction and access
  void buildGlobalRBMs(GenFullM<Scalar> &Xmatrix, Connectivity *cornerToSub); // use null space of (G^T*P_H*G) ... trbm method !!!
  void getGlobalRBM(int iRBM, Scalar *Rvec);
  // R_g matrix-vector multiplication
  void subtractRstar_g(Scalar *u, GenVector<Scalar> &beta); // u -= R_g*beta
  void addRstar_gT(Scalar *u, GenVector<Scalar> &beta); // u += R_g*beta
  // (R_g^T*R_g) matrix assembly
  void assembleRtR(GenFullM<Scalar> &RtRu);

  int *l2g;
  void makeLocalToGlobalDofMap();
  void multAddLT(Scalar *localvec, Scalar *globalvec);
  void multAddLinv(Scalar *localvec, Scalar *globalvec);
  void multLTinv(Scalar *globalvec, Scalar *localvec);
  void multL(Scalar *globalvec, Scalar *localvec);

};

typedef GenSubDomain<double> SubDomain;
#ifdef _TEMPLATE_FIX_
  #include <Driver.d/SubDomain.C>
  #include <Driver.d/BOps.C>
  #include <Driver.d/RbmOps.C>
  #include <Driver.d/LOps.C>
#endif

#endif
