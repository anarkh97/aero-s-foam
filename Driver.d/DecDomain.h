#ifndef _DEC_DOMAIN_H_
#define _DEC_DOMAIN_H_

#include <Feti.d/DistrVector.h>
#include <Feti.d/DistrVectorSet.h>
#include <Paral.d/SubDOp.h>

class Domain;
template <class Scalar> class GenSubDomain;
template <class Scalar> class GenParallelSolver;
template <class Scalar> class GenFetiSolver;
template <class Scalar> class GenDomainGroupTask;
template <class V> class SysState;
template <class Scalar> class DiagParallelSolver;
template <class Scalar> class GenSolver;
template <class Scalar> class GenMDDynamMat;
class DistrGeomState;
class MatrixTimers;
class Connectivity;
class Corotator;
class SubCornerHandler;
class Rbm;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;

template<class Scalar>
class GenDecDomain 
{
 protected:
  Domain *domain;
  MatrixTimers &mt;
  Connectivity *subToNode;
  Connectivity *subToSub;
  Connectivity *mpcToSub_dual;
  Connectivity *mpcToSub_primal;
  Connectivity *subToElem;
  Connectivity *elemToNode;
  Connectivity *nodeToSub;
  Connectivity *mpcToMpc;
  Connectivity *cpuToSub;
  Connectivity *elemToSub;
  GenSubDomain<Scalar> **subDomain;
  DistrInfo internalInfo, internalInfo2;
  DistrInfo nodeInfo;
  Connectivity *grToSub;
  FILE *primalFile; // file to store primal residual

  // only used if user requests stress outputs
  GenDistrVector<Scalar> *stress;
  GenDistrVector<Scalar> *weight;
  Scalar *globalStress;
  Scalar *globalWeight;
 
  int *glSubToLocal;            // local numbering of subdomains
  int *localSubToGl;            // local to global map
  int globalNumSub;             // total number of subdomains on all cpus
  int numSub;                   // local number of subdomains

  FSCommunicator *communicator; // PJSA

  int numCPU, myCPU;            // number of CPUs, my CPU number
  Connectivity *cpuToCPU;       // global problem connectivity
  int *subToCPU;                // subdomain to cpu mapping  
  Connectivity *mpcToCpu;
  int numDualMpc, numPrimalMpc;

  int numWetInterfaceNodes;
  int *wetInterfaceNodes;
  int outFreqCount;
  int outEigCount;

  bool soweredInput;
  int *cornerWeight;

  int sizeSfemStress;
  bool firstOutput; 
 public:
  GenDecDomain(Domain *d);
  virtual ~GenDecDomain();

  GenSubDomain<Scalar>** getAllSubDomains() { return subDomain; }
  GenSubDomain<Scalar>* getSubDomain(int isub) { return subDomain[isub]; }
  Connectivity * getSubToSub() { return subToSub; }
  GenFetiSolver<Scalar> *getFetiSolver();
  GenFetiSolver<Scalar> *getDynamicFetiSolver(GenDomainGroupTask<Scalar> &);
  void buildOps(GenMDDynamMat<Scalar>&, double, double, double, Rbm **rbm = 0, FullSquareMatrix **kelArray = 0, bool make_feti = true);
  DiagParallelSolver<Scalar> *getDiagSolver(int nSub, GenSubDomain<Scalar> **, GenSolver<Scalar> **);
  void rebuildOps(GenMDDynamMat<Scalar>&, double, double, double, FullSquareMatrix** = 0, FullSquareMatrix** = 0);
  void subRebuildOps(int iSub, GenMDDynamMat<Scalar>&, double, double, double, FullSquareMatrix**, FullSquareMatrix**);
  int getNumSub() { return numSub; }
  Connectivity *getMpcToSub() { return mpcToSub_dual; }
  virtual void preProcess();
  virtual void postProcessing(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &f, double eigV = 0.0,
                              GenDistrVector<Scalar> *aeroF = 0, int x = 0, GenMDDynamMat<Scalar> *dynOps = 0,
                              SysState<GenDistrVector<Scalar> > *distState = 0, int ndflag = 0); 
  //virtual void postProcessing(DistrGeomState *u, Corotator ***, double x = 0);
  virtual void postProcessing(DistrGeomState *u, Corotator ***, double x = 0, SysState<GenDistrVector<Scalar> > *distState = 0);
  void setUserDefBC(double *, double *); 
  DistrInfo &solVecInfo() { return internalInfo; } // unconstrained dofs
  DistrInfo &sysVecInfo() { return internalInfo2; } // all dofs
  // user defined control functions
  void extractControlData(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &,
                          GenDistrVector<Scalar> &,
                          Scalar *, Scalar *, Scalar *);
  void initUserDefBC();
  void addUserForce(GenDistrVector<Scalar> &, Scalar*);
  void addCtrl(GenDistrVector<Scalar> &, Scalar*);
  // Non linear functions
  DistrInfo* elementVectorInfo();
  DistrInfo* pbcVectorInfo();
  void scaleDisp(GenDistrVector<Scalar> &u);
  void setNewProperties(int);
  void assignRandMat();
  void retrieveElemset();
  virtual void setsizeSfemStress(int fileNumber);  // YYY DG implementation incomplete: do the element stresses
  virtual int getsizeSfemStress() { return sizeSfemStress; }
  virtual Scalar* getSfemStress(int fileNumber) { return globalStress; } // YYY DG ok for nodal stress, modify for element stresses
  virtual void updateSfemStress(Scalar* str, int fileNumber);
  double computeStabilityTimeStep(GenMDDynamMat<Scalar>&);
  void extractSubDomainMPCs(int iSub);

 protected:
  void makeSubDomains();
  void renumberElements(int iSub);
  void createElemToNode();
  void getSharedNodes();
  void addBMPCs();
  void makeSubToSubEtc();
  void preProcessBCsEtc();
  void distributeBCs();
  void distributeControlLawData();
  void distributeDiscreteMass();
  void renumberBC();
  void preProcessMPCs();
  void getSharedDOFs();
  void getSharedMPCs();
  void makeCorners();
  void makeInternalInfo();
  void makeInternalInfo2();
  void makeNodeInfo();
  void getCPUMap();
  void makeSubDMaps();
  void constructSubDomains(int iSub);
  void makeMpcToMpc();
  void makeGlobalMpcToMpc(Connectivity *procMpcToMpc);
  void makeMpcToSub();
  void buildFFP(GenDistrVector<Scalar> &u, FILE *fffp);
  void makeCornerHandler(int iSub, SubCornerHandler **cornerHandler);
  void setLocalCorners(int iSub, SubCornerHandler **cornerHandler);

 private:
  void initialize();
  // output functions
  void getElementForce(GenDistrVector<Scalar>&, int, int, double);
 public:
  virtual void getStressStrain(GenDistrVector<Scalar>&, int, int, double, int printFlag=0);
 private:
  void getElementStressStrain(GenDistrVector<Scalar>&, int, int, double, int printFlag=0); // YYY DG Implement printFlag
  void getElementAttr(int, int, double);
  void getPrincipalStress(GenDistrVector<Scalar>&, int, int, double);
  void getElementPrincipalStress(GenDistrVector<Scalar> &u, int, int, double);
  void getStressStrain(DistrGeomState *u, Corotator ***, int, int, double);
  void getPrincipalStress(DistrGeomState *u, Corotator ***, int, int, double);
  void getElementPrincipalStress(DistrGeomState *u, Corotator ***, int, int, double);
  void computeSubdElemForce(int iSub, Scalar *globForce,
                            GenDistrVector<Scalar> *u, int Findex);
  void computeSubdStress(int, GenDistrVector<Scalar>*, GenDistrVector<Scalar>*,
                         GenDistrVector<Scalar>*, int, int);
  void computeSubdElemStress(int, Scalar *, GenDistrVector<Scalar> *, int, int);
  void computeSubdStress(int iSub, GenDistrVector<Scalar> *globStress,
                         GenDistrVector<Scalar> *globWeight, DistrGeomState *u,
                         Corotator ***allCorot, int *, int *Findex);
  void getElementStressStrain(DistrGeomState *gs, Corotator ***allCorot,
                              int fileNumber, int Findex, double time);
  void computeSubdElemStress(int iSub, Scalar *glElemStress,
                             DistrGeomState *u, Corotator ***allCorot,
                             int fileNumber, int Findex);
  void outputPrimal(GenDistrVector<Scalar>& primal, int iter);
  void getPrimalVector(int fileNumber, Scalar (*xyz)[11], int numNodes,
                       int ndof, double time);//DofSet::max_known_nonL_dof
  void getPrimalScalar(int fileNumber, Scalar (*xyz)[11], int numNodes,
                       int dof, double time);//DofSet::max_known_nonL_dof
  void getAeroForceScalar(int fileNumber, Scalar (*mergedAeroF)[6],
                          int numNodes, int dof, double time);
  void scaleSubDisp(int iSub, GenDistrVector<Scalar> &u);
  // Helmholtz Fluid functions
  void distribBC(int iSub, GenSubDomain<Scalar> **sd, Domain *domain,
     int *somToSub, int *scaToSub, int *neumToSub, int (*wetToSub)[2],
     int *sBoundFlag);
//  void checkSommerTypeBC(int iSub);
  void buildLocalFFP(int iSub, GenDistrVector<Scalar> *u,
                     Scalar **ffp, int *numSample, double (*dir)[3]);
  void getWError(int iSub, GenDistrVector<Scalar> *u,
                 double *l2err, double *h1err, double *l2, double *h1);
  // coupled_dph functions
  void preProcessFSIs();
  void distributeWetInterfaceNodes();
  void setSubWetInterface(int iSub, int *nWetInterfaceNodesPerSub, int **subWetInterfaceNodes);
  void getSharedFSIs();

  // JLchange coupled_dph functions
  void markSubWetInterface(int iSub, int *nWetInterfaceNodesPerSub, int **subWetInterfaceNodes);
  void addFsiElements();
};

#ifdef _TEMPLATE_FIX_
  #include <Driver.d/DecDomain.C>
#endif

#endif
