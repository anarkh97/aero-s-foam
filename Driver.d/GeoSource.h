#ifndef _GEO_SOURCE_H_
#define _GEO_SOURCE_H_

#include <cstring>
#include <list>
#include <vector>

#include <Element.d/Element.h>
#include <Utils.d/OutputInfo.h>
#include <Math.d/DistVector.h>
#include <Driver.d/StructProp.h>
using namespace NewVec;
#define PI 3.14159265358979

class Connectivity;
class Communicator;
template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
class Domain;
class ControlInterface;
class BinFileHandler;
class MatrixTimer;
class LMPCons;
class Decomposition;
class NLMaterial;
class ControlInfo;
class CoordSet;
class OutputInfo;
class LayMat;
class Attrib;
class EFrameData;
class CoefData;
class LayInfo;
struct Group;

enum {SXX=0,SYY=1,SZZ=2,SXY= 3,SYZ= 4,SXZ= 5,VON=6,
      EXX=7,EYY=8,EZZ=9,EXY=10,EYZ=11,EXZ=12,STRAINVON=13,
      VONTOP=14,VONBOT=15,CONPRESS=16,DAMAGE=17,EFFPSTRN=18};
enum {INX,INY,INZ,AXM,AYM,AZM};
enum {YOUNG,MDENS,THICK};
enum {PSTRESS1=0,PSTRESS2=1,PSTRESS3=2,
      PSTRAIN1=3,PSTRAIN2=4,PSTRAIN3=5,
      PSTRESS1DIREC=6,PSTRESS2DIREC=7,PSTRESS3DIREC=8,
      PSTRAIN1DIREC=9,PSTRAIN2DIREC=10,PSTRAIN3DIREC=11};
// Controls the output of heat fluxes or grad(Temp)
enum {HFLX=0, HFLY=1, HFLZ=2, GRTX=3, GRTY=4, GRTZ=5};
enum {SLDX=0, SLDY=1, SLDZ=2};

struct MatchData
{
  int glPos;
  int elemNum;
  double xi, eta;  // natural coordinates
  bool operator < (const MatchData &v) const { return glPos < v.glPos; }
};

#include <Driver.d/Access.h> // TG OffsetData definition moved inside here because of TEMPLATE_FIX

// Control Law Information class, see Control.d/control.C
// to implement user defined control forces, user defined forces or
// user defined displacements

struct ControlLawInfo
{
  char *fileName;
  char *routineName;
  int numSensor;        // number of sensors
  BCond *sensor;
  int numActuator;      // number of actuators
  BCond *actuator;
  int numUserDisp;      // number of user defined displacements
  BCond *userDisp;
  int numUserForce;     // number of user defined forces
  BCond *userForce;

  ControlLawInfo();
  ~ControlLawInfo();
  void print();
  void makeGlobalClaw(ControlLawInfo *subClaw);
};

typedef NLMaterial *(*MatLoader)(int n, double *params);

struct DoubleList {
   int nval;
   double v[32];
};

struct ltstr
{
  inline bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

class GeoSource {
 protected:
  struct OrderPair {
    int node, pos;
    double xyz[3];
    double &operator[] (int i) { return xyz[i]; }
    bool operator<(const OrderPair &op) const { return node < op.node; }
  };

  int curCluster;

  bool decJustCalled; // if true, no need for an external DECOMPOSITION FILE
  bool exitAfterDec; // if true no need to save Subdomain info for FEM in memory

  // input file names
  char *conName;
  char *geoName;
  char *decName;
  char *mapName;
  char *matchName;

  // output file info
  ControlInfo *cinfo;    // contains nodeset, elemset & timer file name
  ResizeArray<OutputInfo> oinfo; // all output files
  int numOutInfo;                // number of Output requests
  int outLimit; // maximum number of frequencies, eigenvectors or timesteps per file
  int numNodalOutput;   // number of disp output for single nodes
  int *outputNodes;     // list of node for singular output
  int *outNodeIndex;
  int *headLen;  	// header description length for binary output

  int numClusters;
  int maxGlobNode;
  int maxClusNode;
  int numClusNodes;
  int numClusElems;
  int nGlobNodes;
  int numNodes, numInternalNodes;
  CoordSet nodes;
  Elemset elemSet;
  Elemset *packedEsetFluid; //ADDED FOR HEV PROBLEM, EC, 20070820
  int nElem;
  int nElemFluid;      //ADDED FOR HEV PROBLEM, EC, 20070820
  int nElemAllClusters;
  int *allNumClusElems;
  int phantomFlag;

  int *elemTypeNumNodesMap;
  map<int, int> glToPckElems;  // glElemNum  -> packedElem Num

  // Match Data
  int *numMatchData;
  MatchData **matchData;
  int *numGapVecs;
  double (**gapVec)[3];      // gap vectors at struct/fluid interface

  ControlLawInfo *claw;      // user defined control law info

  // Material data
  int numLayInfo;
  ResizeArray<LayInfo *> layInfo;
  int numCoefData;
  ResizeArray<CoefData *> coefData;
  int numLayMat;
  ResizeArray<LayMat *> layMat;

  std::map<const char *, MatLoader, ltstr> userDefMat;
  std::map<int, NLMaterial *> materials;
  std::map<int, int> matUsage;

  int numProps;
  SPropContainer sProps;
  int na;  			// number of attributes
  int namax;
  map<int, Attrib> attrib;
  int maxattrib;

  int numEframes;
  ResizeArray<EFrameData> efd;
  vector<OffsetData> offsets;

  int numCframes;
  ResizeArray<double *> cframes;

  int prsflg;
  int prlflg;
  int constpflg, constqflg; // consistent pressure and gravity

  // Connectivities
  Connectivity *clusToSub;
  Connectivity *subToClus;
  Connectivity *subToSub;
  Connectivity *subToNode;
  Connectivity *subToElem;

  int *subToCPU;
  Connectivity *cpuToSub;
  Connectivity *cpuToCPU;

  double mratio; // consistent-lumped matrix ratio; 1==consistent, 0==lumped

 public:
  int fixedEndM; // 0: don't include fixed end moments in gravity force vector for beams and shells
  int *gl2ClSubMap;

  // BC Data

  int numTextDirichlet;
  BCond *textDBC;
  int numDirichlet;           // number of dirichlet bc
  int numDirichletFluid;           // number of dirichlet bc in fluid, ADDED FOR HEV PROBLEM, EC, 20070820
  BCond *dbc;   // set of those dirichlet bc
  BCond *dbcFluid;   // set of those dirichlet bc in fluid, ADDED FOR HEV PROBLEM, EC, 20070820

  int numTextNeuman;
  BCond *textNBC;
  int numNeuman;              // number of Neuman bc
  BCond *nbc;   // set of Neuman bc

  int numIDis;                // number of Initial displacements
  BCond *iDis;  // set of those initial displacements

  int numIDisModal;
  BCond *iDisModal;

  int numIDis6;               // number of Initial displacements (6 column)
  BCond *iDis6; // set of those intitial displacements

  // PITA
  map<int, std::pair<int, int> > timeSliceOutputFiles;
  BCond *PitaIDis6;    // Array of initial seed displacement array
  int numPitaIDis6;    // # of initial seed initial displacement conditions
  int numTSPitaIDis6;  // # of initial seed displacement vectors
  BCond *PitaIVel6;    // Array of initial seed velocity array
  int numPitaIVel6;    // # of initial seed initial displacement conditions
  int numTSPitaIVel6;  // # of initial seed initial velocity conditions

  int numIVel;                // number of initial velocities
  BCond *iVel;  // set of those initial velocities
  int numIVelModal;
  BCond *iVelModal;

  int numComplexDirichlet;
  ComplexBCond *cdbc;

  int numComplexNeuman;
  ComplexBCond *cnbc;

  bool isShift;
  double shiftV;

  int numDampedModes;   // number of modes that have damping
  BCond *modalDamping;  // the value of damping for those modes with damping reuses BCond class

  Decomposition *optDec;

  map<int, Group> group;
  map<int, list<int> > nodeGroup;

  int numSurfaceDirichlet;
  BCond *surface_dbc;
  int numSurfaceNeuman;
  BCond *surface_nbc;
  int numSurfacePressure;
  BCond *surface_pres;

  bool mpcDirect;

public:
  bool binaryInput, binaryOutput;
  bool binaryInputControlLeft;

  GeoSource(int iniSize = 16);
  virtual ~GeoSource();

  vector<int> localToGlobalElementsNumber;

  // Input Read Functions
  void readCpuToSub();
  int readRanges(BinFileHandler &, int &, int (*&r)[2]);
  void readMatchInfo(BinFileHandler &, int (*)[2], int, int, int *);
  template<class Scalar>
    GenSubDomain<Scalar> *getSubDomain(int glSub, Domain *, int locSub = -1);
  template<class Scalar>
    GenSubDomain<Scalar> *getSubDomain(Domain *, int, Connectivity *);
  void setElemTypeMap();
  void getClusterData(BinFileHandler &);
  void applyAuxData(int *, int *, int, int);
  template<class Scalar>
    void distributeBCs(GenSubDomain<Scalar> *&, int *, int *gl2clNodeMap = 0);
  int getBC(BCond *, int, int *, BCond *&, int *gl2clNodeMap = 0);
  void augmentBC(int, BCond *, BCond *&, int &);
  int getCPUMap(FILE *f, int);
  void createSingleCpuToSub(int numSub);
  int getSubCtrl(BCond *, int, BCond *&, int *, int *, int *&); //bin geo
  int getSubCtrl(BCond *, int, BCond *&, int, int *&); // text geo input
  template<class Scalar>
    void distributeCtrlLaw(GenSubDomain<Scalar> *&, int *, int *);  // binary geo input
  template<class Scalar>
    void distributeCtrlLaw(GenSubDomain<Scalar> *&, int); // text geo input
  template<class Scalar>
    void distributeOutputNodes(GenSubDomain<Scalar> *&, int *, int *);
  template<class Scalar>
    void distributeOutputNodesX(GenSubDomain<Scalar> *, Connectivity *nodeToSub);

  void setGeo(char *file) { geoName = file; }
  void setDecomp(char *file) { decName = file; }
  void setMatch(char *file) { matchName = file; }
  void setCpuMap(char *file) { mapName = file; }
  void setGlob(char *file) { conName = file; }
  void setExitAfterDec(bool exit);
  void setNumLocSub(int);
  void setMatchArrays(int);
  void setTextBC();
  void computeGlobalNumElements();

  // duplicate files for the time parallel method
  void duplicateFilesForPita(int, const int*);

  // decomp read functions
  Connectivity* getDecomposition();
  void getTextDecomp(bool sowering = false);
  void getBinaryDecomp();

  // Parser support Functions
  void setControl(char *checkfile, char *nodeset, char *elemset, char *bcondset = 0);

  // Parser support Functions - Geometry
  int  addNode(int nd, double xyz[3]);
  int  addElem(int en, int type, int nn, int *nodeNumbers);
  int  addMat(int, const StructProp &);
  int  addLay(int, LayInfo *);
  int  addCoefInfo(int, CoefData &);
  CoefData* getCoefData(int i) { assert(i >= 0 && i < numCoefData); return coefData[i]; }
  int  addLayMat(int m, double *);
  int  setAttrib(int n, int a, int ca = -1, int cfrm = -1, double ctheta = 0.0);
  int  setFrame(int,double *);
  int  addCFrame(int,double *);
  void setElementPressure(int, double);
  void setElementPreLoad(int, double);
  void setConsistentPFlag(int);
  void setConsistentQFlag(int, int=1);
  void addOffset(OffsetData &od) { offsets.push_back(od); }

  // Parser support Functions - Boundary Conditions
  int  setDirichlet(int, BCond *);
  void convertHEVDirToHelmDir(); // added to use HEFRS and Helmholtz per Charbel's request
  int  setDirichletFluid(int, BCond *); //ADDED FOR HEV PROBLEM, EC, 20070820
  int  setNeuman(int, BCond *);
  int  setIDis(int, BCond *);
  int  setIDisModal(int, BCond *);
  int  setIDis6(int, BCond *);
  int  setIVel(int, BCond *);
  int  setIVelModal(int, BCond *);
  int  addSurfaceDirichlet(int, BCond *);
  int  addSurfaceNeuman(int, BCond *);
  int  addSurfacePressure(int, BCond *);

  // PITA
  int getLocalTimeSliceCount() { return timeSliceOutputFiles.size(); }
  std::pair<int, int> getTimeSliceOutputFileIndices(int timeSliceRank);
  // Set initial seed conditions
  int setPitaIDis6(int, BCond *, int);
  int setPitaIVel6(int, BCond *, int);

  // Parser support Function - Damping
  int setModalDamping(int, BCond *);

  // Parser support Functions - "New" Material Specifications
  void addMaterial(int i, NLMaterial *m) { materials[i] = m; }
  void addMaterial(int i, const char *, DoubleList &d);
  void loadMaterial(const char *, const char *);
  void setMatUsage(int i, int t) { matUsage[i] = t; }
  // Parser support Functions - Output Functions
  void addOutput(OutputInfo &);

  // Parser support Functions - User Defined Control Functions
  void setControlFile(char *_filename);
  void setControlRoutine(char *_routinename);
  int  setSensorLocations(int, BCond *);
  int  setActuatorLocations(int, BCond *);
  int  setUsddLocation(int, BCond *);
  int  setUsdfLocation(int, BCond *);

  void setUpData();
  void setUpData(CoordSet &, Domain *, int);

  Elemset* getPackedEsetFluid() { return packedEsetFluid; }  //ADDED FOR HEV PROBLEM, EC, 20070820

  // PITA: Access to initial seed conditions
  int getNumPitaIDis6() const  { return numPitaIDis6; }     // # of initial seed initial displacement conditions
  int getNumPitaIVel6() const  { return numPitaIVel6; }     // # of initial seed initial velocity conditions
  int getNumTSPitaIDis6() const { return numTSPitaIDis6; }   // # of initial seed displacement vectors
  int getNumTSPitaIVel6() const { return numTSPitaIVel6; }   // # of initial seed velocity vectors
  int getUserProvidedSeedCount() const { return std::min(numTSPitaIDis6, numTSPitaIVel6); }
  const BCond &getPitaIDis6(int i) const { return PitaIDis6[i]; }  // Initial seed displacement for each time-slice
  BCond &getPitaIDis6(int i) { return PitaIDis6[i]; } 
  const BCond &getPitaIVel6(int i) const { return PitaIVel6[i]; }  // Initial seed velocity for each time-slice
  BCond &getPitaIVel6(int i) { return PitaIVel6[i]; } 

  ControlInfo *getCheckFileInfo()  { return cinfo; }
  int  getNumClusNodes()  { return numClusNodes; }
  int  getNodes(CoordSet &);
  CoordSet&  GetNodes(); // HB: return the node CoordSet of GeoSource
  int  getMaxNodeNum()  { return maxGlobNode; }
  int  getNumGlobNodes()  { return nGlobNodes; }
  void  setNumGlobNodes(int n) { nGlobNodes = n; }
  int  getElems(Elemset &, int = 0, int * = 0);
  int getNonMpcElems(Elemset &eset);
  int  getNumClusters() { return numClusters; }
  char *getCpuMapFile()  { return mapName; }
  char *getMatchFileName()  { return matchName; }
  int  numElem() { return nElem; }
  int  numElemFluid() { return nElemFluid; }   //ADDED FOR HEV PROBLEM, EC 20070820
  int  numNode() { return numNodes; }
  void setNumNodes(int n) { numNodes = n; }
  int  totalNumNodes() { return numNodes + numInternalNodes; }
  //int  getPhantomFlag()  { return phantomFlag; }
  //int  glToPack(int i) { return glToPck[i]; }
  int  glToPackElem(int i);
  Connectivity *getClusToSub()  { return clusToSub; }
  Connectivity *getSubToClus()  { return subToClus; }
  Connectivity *getSubToSub()  { return subToSub; }
  Connectivity *getSubToElem()  { return subToElem; }
  void setSubToElem(Connectivity *ste) { subToElem = ste; }
  Connectivity *getSubToNode()  { return subToNode; }

  LayInfo *getLayerInfo(int num)  { return layInfo[num]; }
  SPropContainer &getStructProps()  { return sProps; }
  EFrameData *getEframes()  { return efd+0; }
  double **getCframes()  { return cframes+0; }
  LayInfo **getLayInfo() { return layInfo+0; }
  int getNumLayInfo() { return numLayInfo; }
  int getNumEframes() { return numEframes; }
  int getNumCframes() { return numCframes; }
  int pressureFlag() { return prsflg; }
  int preloadFlag() { return prlflg; }
  int consistentPFlag() { return constpflg; }
  int consistentQFlag() { return constqflg; }
  map<int, Attrib> &getAttributes()  { return attrib; }
  int getNumAttributes() { return na; }

  int getNumDirichlet()  { return numDirichlet; }
  int getNumDirichletFluid()  { return numDirichletFluid; } //ADDED FOR HEV PROBLEM, EC, 20070820
  int getNumNeuman()  { return numNeuman; }
  int getNumIDisModal() { return numIDisModal; }
  int getDirichletBC(BCond *&);
  int getDirichletBCFluid(BCond *&);
  int getTextDirichletBC(BCond *&);
  int getNeumanBC(BCond *&);
  int getTextNeumanBC(BCond *&);
  int getIDis(BCond *&);
  int getIDisModal(BCond *&);
  int getIDis6(BCond *&);
  int getIVel(BCond *&);
  int getIVelModal(BCond *&);
  int getITemp(BCond *&);
  int getNumProps() { return numProps; }
  void setNumProps(int n) { numProps = n; }

  int getSurfaceDirichletBC(BCond *&);
  int getSurfaceNeumanBC(BCond *&);
  int getSurfacePressure(BCond *&);

  int getModalDamping(BCond *&);

  int *getNumMatchData()  { return numMatchData; }
  MatchData *getMatchData(int iSub)  { return matchData[iSub]; }
  double (*getGapVecs(int iSub))[3]  { return gapVec[iSub]; }

  int getNumOutInfo()  { return numOutInfo; }
  void setOutLimit(int _outLimit) { outLimit = _outLimit; }
  int getOutLimit() { return outLimit; }
  void setNumNodalOutput();
  int getNumNodalOutput() { return numNodalOutput; }
  OutputInfo *getOutputInfo()  { return oinfo+0; }
  bool elemOutput();

  int *getSubToCPU()  { return subToCPU; }
  void setCpuToSub(Connectivity *c) { cpuToSub=c; }
  Connectivity *getCpuToSub() { return cpuToSub; }
  Connectivity *getCpuTOCPU()  { return cpuToCPU; }

  ControlLawInfo *getControlLaw()  { return claw; }
//  void makeGlobalClaw(ControlLawInfo *subClaw);
  ControlInterface *getUserSuppliedFunction();

  Elemset* getElemSet(void){return(&elemSet);}

  void simpleDecomposition(int numSubdomains, bool estFlag, bool weightOutFlag); // dec

  // Output Functions
  template<int bound>
    void outputNodeVectors(int, double (*)[bound], int, double time = -1.0);//DofSet::max_known_nonL_dof
  template<int bound>
    void outputNodeVectors(int, DComplex (*)[bound], int, double time = -1.0);
  template<int bound>
    void outputNodeVectors6(int, double (*)[bound], int, double time = -1.0);
  template<int bound>
    void outputNodeVectors6(int, DComplex (*)[bound], int, double time = -1.0);
  //void outputNodeVectors6(int, double (*)[11], int, double time = -1.0) {};
  //void outputNodeVectors6(int, DComplex (*)[11], int, double time = -1.0) {};
  void outputNodeScalars(int, double *, int, double time = -1.0);
  void outputNodeScalars(int, DComplex *, int, double time = -1.0);
  void outputEnergies(int, double, double, double, double, double, double, double);
  void outputEnergies(int, double, DComplex, DComplex, DComplex, DComplex, DComplex, DComplex);
  void outputElemVectors(int, double *, int, double time = -1.0);
  void outputElemVectors(int, DComplex *, int, double time = -1.0);
  void outputElemStress(int, double *, int, int *, double = -1.0);
  void outputElemStress(int, DComplex *, int, int *, double = -1.0);
  void openOutputFiles(int *outNodes = 0, int *outIndex = 0, int num = 0);
  void openOutputFilesForPita(int sliceRank);
  void closeOutputFiles();
  void closeOutputFilesForPita(int sliceRank);
  void createBinaryOutputFile(int, int, int iter = 0);
  void outputHeader(BinFileHandler&, int, int);
  void setHeaderLen(int);
  void outputHeader(int);
  void outputRange(int, int *, int, int, int , int iter = 0);//CBM
  int getHeaderDescription(char *, int);

  void readGlobalBinaryData();
  void computeClusterInfo(int glSub);

  void writeDistributedInputFiles(int nCluster, Domain*);
#ifdef SOWER_SURFS
  template<class Scalar>
    void readDistributedSurfs(int subNum);
#endif
  template<class Scalar>
    GenSubDomain<Scalar> * readDistributedInputFiles(int localSubNum, int subNum);

  // PJSA: new output functions (implemented in Driver.d/BinaryOutput.C and Driver.d/BinaryOutputInclude.C)
  BinFileHandler* openBinaryOutputFile(int glSub, int fileNumber, int iter = 0);
  void writeNodeScalarToFile(double *data, int numData, int glSub, int offset, int fileNumber, int iter,
                             int numRes, double time, int numComponents, int *glNodeNums);
  void writeNodeScalarToFile(DComplex *complexData, int numData, int glSub, int offset, int fileNumber, int iter,
                             int numRes, double time, int numComponents, int *glNodeNums);
  template<class Scalar, int dim>
    void writeNodeVectorToFile(SVec<Scalar, dim> &, int glSub, int offset, int fileNumber, int iter,
                               int numRes, double time, int numComponents, int startComponent, int *glNodeNums);
  void writeElemScalarToFile(double *data, int numData, int glSub, int offset, int fileNumber, int iter,
                             int numRes, double time, int totData, int *glElemNums);
  void writeElemScalarToFile(DComplex *complexData, int numData, int glSub, int offset, int fileNumber, int iter,
                             int numRes, double time, int totData, int *glElemNums);


  // Shifting functions
  bool isShifted() { return isShift; }
  void initShift() {
    if(!isShift) { isShift = true; shiftV = 0.0; }
  }
  void setShift(double w2) { isShift = true; shiftV = w2; }
  void setImpe(double f) { isShift = true; shiftV = 4.0*PI*PI*f*f; }
  void setOmega(double w) { isShift = true; shiftV = w*w; }
  void resetShift(double w2) { isShift = false; shiftV = w2; } //CBM--exp
  double shiftVal() { return shiftV; }
  double freq() { return sqrt(shiftV)/(2.0*PI); }
  double omega() { return sqrt(shiftV); }
  double kappa() { /*if(numProps > 1) cerr << "warning: assuming homogenous fluid (attr #1), k = " << sProps[0].kappaHelm << endl;*/ return sProps[0].kappaHelm; }

  void setMRatio(double _mratio) { assert(_mratio >= 0.0 && _mratio <= 1.0); mratio = _mratio; }
  double getMRatio() const { return mratio; }

  // Housekeeping functions
  void cleanUp();
  void cleanAuxData();

  void makeDirectMPCs(int numLMPC, ResizeArray<LMPCons *> &lmpc);
  int reduceMPCs(int numLMPC, ResizeArray<LMPCons *> &lmpc);
  void addMpcElements(int numLMPC, ResizeArray<LMPCons *> &lmpc);
  void addFsiElements(int numFSI, ResizeArray<LMPCons *> &fsi);
  bool setDirectMPC(bool mode) { return mpcDirect = mode; }
  /// Whether we are doing direct elimination for MPCs
  bool getDirectMPC() { return mpcDirect; }
  Element* getElem(int topid) { return elemSet[topid]; }

  double global_average_E, global_average_nu, global_average_rhof;

  void makeEframe(int ele, int refnode, double *d);

// Group stuff
   void setGroupAttribute(int a, int g);
   void setNodeGroup(int nn, int id);

// Sfem stuff
  enum Rprop { A, E, NU, RHO, T, KX, KY, KZ }; // sfem
  void setGroupRandomProperty(int g, Rprop prop_type, double mean, double std_dev);
  void printGroups();

  // POD-ROM sample nodes
  int sampleNodeCount() const { return sampleNode_.size(); }
  void sampleNodeAdd(int id) { sampleNode_.push_back(id); }

  typedef std::vector<int> SampleNodeList; 
  SampleNodeList::const_iterator sampleNodeBegin() const { return sampleNode_.begin(); }
  SampleNodeList::const_iterator sampleNodeEnd()   const { return sampleNode_.end();   }

private:
  SampleNodeList sampleNode_;

protected:
  void closeOutputFileImpl(int fileIndex);
};


struct RandomProperty
{
  GeoSource::Rprop rprop;
  double mean;
  double std_dev;
  RandomProperty(GeoSource::Rprop rp, double m, double sd) { rprop = rp; mean = m; std_dev = sd; }
  RandomProperty(const RandomProperty &rp) { rprop = rp.rprop; mean = rp.mean; std_dev = rp.std_dev; }
};

struct Group
{
  vector<int> attributes;
  vector<RandomProperty> randomProperties;
};

#ifdef _TEMPLATE_FIX_
 #include <Driver.d/GeoSource.C>
 #ifdef DISTRIBUTED
  #include <Driver.d/BinaryOutputInclude.C>
 #endif
#endif

#endif
