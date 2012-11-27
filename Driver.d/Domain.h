#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include <cassert>
#include <set>
#include <map>

#include <Utils.d/resize_array.h>
#include <Utils.d/SolverInfo.h>
#include <Utils.d/CompositeInfo.h>
#include <Utils.d/MyComplex.h>
#include <Timers.d/Timing.h>
#include <Timers.d/MatrixTimers.h>
#include <Driver.d/HData.h>
#include <Parser.d/AuxDefs.h>
#include <Sfem.d/Sfem.h>
#include <Utils.d/OutputInfo.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>

#include <map>

class MortarHandler;
class MFTTData;
class ControlInterface;
class ControlInfo;
class DofSetArray;
class ConstrainedDSA;
class DOFMap;
template <class Scalar> class GenNBSparseMatrix;
typedef GenNBSparseMatrix<double> NBSparseMatrix;
template <class Scalar> class GenDBSparseMatrix;
typedef GenDBSparseMatrix<double> DBSparseMatrix;
typedef GenDBSparseMatrix<DComplex> DBComplexSparseMatrix;
template <typename Scalar, typename SolverClass> class GenEiSparseMatrix;
template <class Scalar> class GenSkyMatrix;
typedef GenSkyMatrix<double> SkyMatrix;
typedef GenSkyMatrix<DComplex> SkyMatrixC;
template <class Scalar> class GenBlockSky;
typedef GenBlockSky<double> BlockSky;
class Connectivity;
template <class Scalar> class GenCuCSparse;
typedef GenCuCSparse<double> CuCSparse;
template <class Scalar> class GenBLKSparseMatrix;
typedef GenBLKSparseMatrix<double> BLKSparseMatrix;
template <class BaseSolver, class Scalar> class GoldfarbIdnaniQpSolver;
template <class Scalar> class GenDynamMat;
typedef GenDynamMat<double> DynamMat;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
template <class Scalar> class GenDecDomain;
template <class Scalar> class GenDistrDomain;
template <class Scalar> class GenSandiaDomain;
class FlExchanger;
class Rbm;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
typedef GenSolver<DComplex> ComplexSolver;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
typedef GenSparseMatrix<DComplex> ComplexSparseMatrix;
template <class Scalar, class AnyVector> class KrylovProjector;
template <class AnyVector> class Preconditioner;
template <class Scalar, class AnyVector, class AnyOperator> class GenPCGSolver;
template <class Scalar> class GenSpoolesSolver;
typedef GenSpoolesSolver<double> SpoolesSolver;
template <class Scalar> class GenMumpsSolver; 	//Axel
typedef GenMumpsSolver<double> MumpsSolver;   	//Axel
class GeomState;
class DistrGeomState;
class IntFullM;
class ControlLawInfo;
class StaticTimers;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
template <class Scalar> class GenSubDOp;
typedef GenSubDOp<double> SubDOp;
template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
class FSCommunicator;

namespace Rom {
template <typename Scalar> class GenGalerkinProjectionSolver;
template <typename Scalar> class GenEiSparseGalerkinProjectionSolver;
} /* end namespace Rom */

// HB
class SurfaceEntity;

extern Sfem *sfem;

const double defaultTemp = -10000000.0;

// ... Structure used to store problem Operators buildSkyOps
// ... i.e. Only a Solver is needed for a static problem
template<class Scalar>
struct AllOps
{
  GenSolver<Scalar> *sysSolver;  // system solver: to solve (coeM*M+coeC*C+coeK*K)x = b
  GenSparseMatrix<Scalar> *spm;
  GenSolver<Scalar> *prec;       // preconditioner
  GenSparseMatrix<Scalar> *spp;
  
  GenSparseMatrix<Scalar> *Msolver;  // for assembling mass solver: to solve Mx = b
  GenSparseMatrix<Scalar> *K;    // stiffness matrix
  GenSparseMatrix<Scalar> *M;    // mass matrix
  GenSparseMatrix<Scalar> *C;	 // damping matrix
  GenSparseMatrix<Scalar> *Kuc;	 // constrained to unconstrained stiffness
  GenSparseMatrix<Scalar> *Muc;	 // constrained to unconstrained mass matrix
  GenSparseMatrix<Scalar> *Cuc;	 // constrained to unconstrained damping matrix
  GenSparseMatrix<Scalar> *Kcc;  // constrained to constrained stiffness matrix
  GenSparseMatrix<Scalar> *Mcc;	 // constrained to constrained mass matrix
  GenSparseMatrix<Scalar> *Ccc;  // constrained to constrained damping matrix
  GenSparseMatrix<Scalar> **C_deriv;    // derivatives of damping matrix for higher order sommerfeld
  GenSparseMatrix<Scalar> **Cuc_deriv;    // derivatives of constrained to unconstrained damping matrix for higher order sommerfeld

  GenVector<Scalar> *rhs_inpc;
  // Constructor
  AllOps() { sysSolver = 0; spm = 0; prec = 0; spp = 0; Msolver = 0; K = 0; M = 0; C = 0; Kuc = 0; Muc = 0; Cuc = 0; Kcc = 0; Mcc = 0; Ccc = 0; C_deriv = 0; Cuc_deriv = 0; rhs_inpc = 0;}

  void zero() {if(K) K->zeroAll();
               if(M) M->zeroAll();
               if(C) C->zeroAll();
               if(Kuc) Kuc->zeroAll();
               if(Muc) Muc->zeroAll();
               if(Cuc) Cuc->zeroAll();
               if(Kcc) Kcc->zeroAll();
               if(Mcc) Mcc->zeroAll();
               if(Ccc) Ccc->zeroAll();
             }
};

// Structure to store discrete masses
#include <Driver.d/DMassData.h> // TG : struct DMassData moved here

// Structure to store dof renumbering information
struct Renumber {
  int *order;
  int *renumb;
  Renumber() { order = renumb = 0; }
};

struct PrevFrc {
   int lastTIndex;
   double lastFluidTime;
   Vector lastFluidLoad;

   PrevFrc(int neq) : lastFluidLoad(neq, 0.0) { lastTIndex = -1; }
};

/** Class representing a structure and containing all auxiliary data-structures
 *
 */
class Domain : public HData {
  protected:
     MatrixTimers *matrixTimers;// timers to time factoring and assembly
     Solver *solver;            // this domains solver
     int numnodes;		// number of nodes
     int numnodesFluid;		// number of nodes for fluid, ADDED FOR HEV PROBLEM, EC, 20070820
     CoordSet &nodes;   	// All the nodes
     int numele;		// number of elements
     Elemset packedEset;	// The set of compressed elements
     int numdofs; 		// the total number of degrees of freedom

     // BC related data members
     int numDirichlet;		// number of dirichlet bc
     int numDirichletFluid;	// number of dirichlet bc in fluid, ADDED FOR HEV PROBLEM, 20070820
     int numDispDirichlet;      // number of displacement dirichlet bc
     BCond* dbc;		// set of those dirichlet bc
     BCond* dbcFluid;		// set of those dirichlet bc in fluid, ADDED FOR HEV PROBLEM, EC, 20070820
     int numLMPC; 		// total number of Linear Multi-Point Constraints (both real & complex)
     ResizeArray<LMPCons *> lmpc;  // set of LMPCs
     int numFSI; 		// total number of Fluid-Structure interaction Constraints
     ResizeArray<LMPCons *> fsi ;  // set of FSIs
     int numCTC;                // total number of contact constraints
     int numNeuman;		// number of Neuman bc
     BCond* nbc;		// set of Neuman bc
     int numIDis;		// number of Initial displacements
     BCond *iDis;		// set of those initial displacements
     int numIDisModal;
     BCond *iDisModal;
     int numIDis6;		// number of Initial displacements (6 column)
     BCond* iDis6;              // set of those intitial displacements
     int numIVel;		// number of initial velocities
     BCond *iVel;		// set of those initial velocities
     int numIVelModal;
     BCond *iVelModal;

     DofSetArray *dsa;		// Dof set array
     DofSetArray *dsaFluid;	// Dof set array for fluid, ADDED FOR HEV PROBLEM, EC 20070820
     ConstrainedDSA *c_dsa;	// Constrained dof set array
     ConstrainedDSA *c_dsaFluid;// Constrained dof set array for fluid, ADDED FOR HEV PROBLEM, EC 20070820
     ConstrainedDSA *MpcDSA;
     ConstrainedDSA *g_dsa;
     Connectivity *allDOFs;     // all dof arrays for each node
     Connectivity *allDOFsFluid;     // all dof arrays for each node
     int maxNumDOFs; 		// maximum number of dofs an element has
     int maxNumDOFsFluid; 	// maximum number of dofs a fluid element has, ADDED FOR HEV PROBLEM, EC, 20070820
     int maxNumNodes;           // maximum number of nodes an element has
     int maxNumNodesFluid;           // maximum number of nodes a fluid element has, ADDED FOR HEV PRBLEM, EC, 20070820
     Connectivity *elemToNode, *nodeToElem, *nodeToNode, *fsiToNode, *nodeToFsi, *nodeToNodeDirect;
     Connectivity *elemToNodeFluid, *nodeToElemFluid, *nodeToNodeFluid;
                                // ADDED FOR HEV PROBLEM, EC, 20070820
     Connectivity *nodeToNode_sommer; // for higher order sommerfeld
     SolverInfo sinfo;		// solver information structure
     ControlLawInfo *claw;      // contains user defined routine information
     DMassData *firstDiMass;	// Discrete Mass
     int nDimass;               // number of DMASS // TG: for sower
     double *gravityAcceleration;   // (gx,gy,gz)
     double gravitySloshing;    // g, added for sloshing problem, EC, 20070724
     compStruct renumb;         // renumbered nodes per structural component
     compStruct renumbFluid;    // renumbered nodes per fluid component, ADDED FOR HEV PROBLEM, 20070820
     compStruct renumb_nompc;
     MFTTData *mftval;		// Mechanics Force Time Table
     MFTTData *mptval;		// Mechanics Pressure Time Table
     MFTTData *hftval;		// Heat Fluxes Time Table
     ResizeArray<MFTTData *> ymtt;         // PJSA: Young's modulus vs. temperature table
     int numYMTT;                          // number of YM Temp tables
     ResizeArray<MFTTData *> ctett;        // PJSA: Coefficient of thermal expansion vs. temperatur table
     int numCTETT;                         // number of CTE Temp tables
     FlExchanger *flExchanger;  // Fluid Exchanger
     FILE *outFile;

     // for computing stress results in function getStressStrain(...)
     Vector *stress;
     Vector *weight;
     Vector *elDisp;
     Vector *elTemp;
     Vector *elstress;
     Vector *elweight;
     // for computing stress results in function getPrincipalStress(...)
     FullM *p_stress;
     FullM *p_elstress;
     Vector *stressAllElems; // stores stresses of all the elements : used Sfem
     int sizeSfemStress;

     // for compute energies
     double Wext;
     double Waero;
     double Wdmp;
     double pWela;
     double pWkin;
     Vector *previousExtForce;
     Vector *previousAeroForce;
     Vector *previousDisp;
     Vector *previousCq;

     Vector *heatflux;
     Vector *elheatflux;

     //ADDED FOR SLOSHING PROBLEM, EC, 20070723
     Vector *fluidDispSlosh;
     Vector *elFluidDispSlosh;
     Vector *elPotSlosh;
     Vector *elFluidDispSloshAll;

    double savedFreq;
    bool sowering;
    bool output_match_in_top;

    void writeTopFileElementSets(ControlInfo *cinfo, int * nodeTable, int* nodeNumber, int topFlag);

    Elemset elems_copy; // Needed for SFEM
    Elemset elems_fullcopy; // copy of the full elementset, needed for SFEM

    StructProp *p; // property for new constraints

  public:
    // Implements nonlinear dynamics postprocessing for file # fileId
    void postProcessingImpl(int fileId, GeomState*, Vector&, Vector&,
                            double, int, double *, double *,
                            Corotator **, FullSquareMatrix *, double *acceleration = 0,
                            double *acx = 0, GeomState *refState=0, Vector *reactions=0);

     Domain(int iniSize = 16);
     Domain(Domain &, int nele, int *ele, int nnodes, int *nodes);
     Domain(Domain &, Elemset *elems, CoordSet *nodes);  // PJSA: for new sower
     virtual ~Domain();

     double getSavedFreq()  { return savedFreq; }
     void setSavedFreq(double freq)  { savedFreq = freq; }
     double *temprcvd;          // temperature received by structure from
                                // heat solution
     int numContactPairs;  // used for timing file
     char * optinputfile;

     int* glWetNodeMap; //HB

     int*  umap_add;            // mapping for coupling matrix assembly, ADDED FOR HEV PROBLEM, EC, 20070820
     int*  umap;                // mapping for coupling matrix, ADDED FOR HEV PROBLEM, EC, 20070820
     int** pmap;                // mapping for coupling matrix, ADDED FOR HEV PROBLEM, EC, 20070820
     int   nuNonZero;            // ADDED FOR HEV PROBLEM, EC, 20070820
     int*  npNonZero;            // ADDED FOR HEV PROBLEM, EC, 20070820
     double ** C_condensed;

     // functions for controlling printing to screen
     void setVerbose() { outFile = stderr; }
     void setSilent()  { outFile = 0;      }
     void setOutputMatchInTop(bool b) {output_match_in_top = b;};
     void readInModes(char* modesFileName);
     void setSowering(bool b) { sowering = b;}
     bool getSowering() { return sowering;}
     void make_bc(int *bc, double *bcx);
     void make_bc(int *bc, DComplex *bcxC);
     void makeAllDOFs();
     void makeAllDOFsFluid();  //ADDED FOR HEV PROBLEM, EC, 20070820

     void createKelArray(FullSquareMatrix *& kel);
     void createKelArray(FullSquareMatrix *& kel,FullSquareMatrix *& mel);
     void createKelArray(FullSquareMatrix *&kArray, FullSquareMatrix *&mArray, FullSquareMatrix *&cArray);
     void reBuild(int iteration, FullSquareMatrix *kel);
     void getElasticStiff(FullSquareMatrix *kel);
     void getElementMass(FullSquareMatrix *kel);
     void getElementForces( GeomState &geomState, Corotator **allCorot,
                            int fileNumber, int forceIndex, double time);
     void getStressStrain(GeomState &geomState, Corotator **allCorot,
                          int fileNumber, int stressIndex, double time, 
                          GeomState *refState = NULL);
     void getPrincipalStress(GeomState &geomState, Corotator **allCorot,
                          int fileNumber, int stressIndex, double time);
     void updateStates(GeomState *refState, GeomState& geomState, Corotator **allCorot);
     void getElemStiffAndForce(const GeomState &geomState, double time, 
                               const GeomState *refState, const Corotator &elemCorot,
                               double *elemForce, FullSquareMatrix &elemStiff);
     void getElemStiffAndForce(const GeomState &geomState, double time, 
                               const Corotator &elemCorot,
                               double *elemForce, FullSquareMatrix &elemStiff);
     void getStiffAndForce(GeomState &u, Vector &elementInternalForce,
                           Corotator **allCorot, FullSquareMatrix *kel,
                           Vector &residual, double lambda = 1.0, double time = 0.0,
                           GeomState *refState = NULL, Vector *reactions = NULL);
     void getFollowerForce(GeomState &u, Vector &elementInternalForce,
                           Corotator **allCorot, FullSquareMatrix *kel,
                           Vector &residual, double lambda = 1.0, double time = 0.0,
                           GeomState *refState = NULL, Vector *reactions = NULL);
     void getWeightedStiffAndForceOnly(const std::map<int, double> &weights,
                                       GeomState &u, Vector &elementInternalForce,
                                       Corotator **allCorot, FullSquareMatrix *kel,
                                       Vector &residual, double lambda, double time,
                                       GeomState *refState);
     void getElemInternalForce(const GeomState &geomState, double time,
                               const GeomState *refState, const Corotator &elemCorot,
                               double *elemForce, FullSquareMatrix &elemStiff);
     void getElemInternalForce(const GeomState &geomState, double time,
                               const Corotator &elemCorot,
                               double *elemForce, FullSquareMatrix &elemStiff);
     void getInternalForce(GeomState &u, Vector &elementInternalForce,
                           Corotator **allCorot, FullSquareMatrix *kel,
                           Vector &residual, double lambda = 1.0, double time = 0.0,
                           GeomState *refState = NULL, Vector *reactions = NULL);
     void getWeightedInternalForceOnly(const std::map<int, double> &weights,
                                       GeomState &u, Vector &elementInternalForce,
                                       Corotator **allCorot, FullSquareMatrix *kel,
                                       Vector &residual, double lambda, double time,
                                       GeomState *refState);

     void applyResidualCorrection(GeomState &geomState, Corotator **corotators, Vector &residual, double rcoef = 1.0);
     void initializeParameters(GeomState &geomState, Corotator **corotators);
     void updateParameters(GeomState &geomState, Corotator **corotators);
     double getError(Corotator **corotators);

     void getGeometricStiffness(GeomState &u, Vector &elementInternalForce,
        			Corotator **allCorot, FullSquareMatrix *&kel);
     void computeGeometricPreStress(Corotator **&allCorot, GeomState *&geomState,
                             FullSquareMatrix *&kelArray, StaticTimers *times,
                             FullSquareMatrix *&geomKelArray, FullSquareMatrix *&melArray,
                             bool melFlag = false);
     ControlInterface* getUserSuppliedFunction();
     void setsizeSfemStress(int fileNumber);
     int getsizeSfemStress() { return sizeSfemStress; }
     double * getSfemStress(int fileNumber, double* dummystress); // dummystress is an arbitrary vector that tells that the stress is double*
     DComplex * getSfemStress(int fileNumber, DComplex* dummystress) {cerr <<"DComplex * Domain::getSfemStress not implemented" << endl; return 0;};
     void updateSfemStress(double* str, int fileNumber);
     void updateSfemStress(DComplex* str, int fileNumber) {cerr <<"Domain::updateSfemStress(DComplex*,.) not implemented" << endl;};
     // Functions to suport the parsing:
     int  addDMass(int, int, double, int jdof = -1);
     DMassData* getFirstDMassData() { return firstDiMass; }
     int getDMassNumber() { return nDimass; }
     int  setDirichlet(int,BCond *);
     int  setDirichletFluid(int,BCond *);
     int  addLMPC(LMPCons *, bool checkflag=true);
     void printLMPC();
     void printLMPC2();
     void normalizeLMPC();
     void setPrimalLMPCs(int &numDual, int &numPrimal);
     void checkLMPCs(Connectivity *nodeToSub);
     void printFSI(FILE* file=stderr);
     Connectivity * makeLmpcToNode();
     Connectivity * makeLmpcToNode_primal();
     void makeFsiToNode();
     Connectivity *getFsiToNode() { return fsiToNode; }
     int getNumFSI() { return numFSI; }
     void setNumFSI(int n) { numFSI = n; }
     ResizeArray<LMPCons *> &getFSI() { return fsi; }
     virtual int  setNeuman(int,BCond *);
     int  setIDis6(int, BCond *);
     int  setIDisModal(int, BCond *);
     int  setIDis(int, BCond *);
     int  setIVel(int, BCond *);
     int  setIVelModal(int, BCond *);
     int  setIAcc(int, BCond *);
     int  setMFTT(MFTTData *);
     MFTTData * getMFTT() const { return mftval; }
     int  setMPTT(MFTTData *);
     int  setHFTT(MFTTData *);
     int  addYMTT(MFTTData *);
     void printYMTT();
     int  addCTETT(MFTTData *);
     std::pair<int, ResizeArray<MFTTData*>* >* getCTETT() { return new std::pair<int, ResizeArray<MFTTData*>* >(numCTETT,&ctett); };
     std::pair<int, ResizeArray<MFTTData*>* >* getYMTT() { return new std::pair<int, ResizeArray<MFTTData*>* >(numYMTT,&ymtt); };
     void printCTETT();
     void computeTDProps();

     virtual void setUpData();
     SolverInfo  &solInfo() { return sinfo; }
     MatrixTimers &getTimers() { return *matrixTimers; }
     void setGravity(double ax, double ay, double az);

     void setGravitySloshing(double gg);

     virtual int glToPackElem(int e);
     void ProcessSurfaceBCs();
     void setNewProperties(int s);
     void assignRandMat();
     void retrieveElemset();

     /** Abstract method to assemble any type of operator
      *
      * the OpList takes care of distributing the contribution to the various components.
      */
     template<class Scalar, class OpList>
     void assembleSparseOps(OpList &ops);
/** ... General build functions to replace the specialized build
  * ... functions and allow us to reuse the code in each problem
  * ... type (i.e. use makeSparseOps in statics, dynamics, eigen, etc.)
  */
     template<class Scalar>
       void buildOps(AllOps<Scalar> &ops, double Kcoef, double Mcoef, double Ccoef,
                     Rbm *rbm = 0, FullSquareMatrix *kelArray = 0, FullSquareMatrix *melArray = 0,
                     FullSquareMatrix *celArray = 0, bool factor = true);

     template<class Scalar>
       void makeStaticOpsAndSolver(AllOps<Scalar> &ops, double Kcoef, double Mcoef,
                 double Ccoef, GenSolver<Scalar> *&systemSolver, GenSparseMatrix<Scalar> *&spm,
                 Rbm *rbm = 0, FullSquareMatrix *kelArray = 0, FullSquareMatrix *melArray = 0,
                 FullSquareMatrix *celArray = 0);

     template<class Scalar>
       void makeDynamicOpsAndSolver(AllOps<Scalar> &ops, double Kcoef, double Mcoef,
                 double Ccoef, GenSolver<Scalar> *&systemSolver, GenSparseMatrix<Scalar> *&spm,
                 Rbm *rbm = 0, FullSquareMatrix *kelArray = 0, FullSquareMatrix *mel = 0,
                 FullSquareMatrix *celArray = 0);

     template<class Scalar>
       void rebuildOps(AllOps<Scalar> &ops, double Kcoef, double Mcoef, double Ccoef,
	  	       Rbm* rbm = 0, FullSquareMatrix *kelArray = 0, FullSquareMatrix *mel = 0,
                       FullSquareMatrix *celArray = 0, bool factor = true);

     template<class Scalar>
       void makeSparseOps(AllOps<Scalar> &ops, double Kcoef, double Mcoef,
	 		  double Ccoef, GenSparseMatrix<Scalar> *mat = 0,
                          FullSquareMatrix *kelArray = 0, FullSquareMatrix *melArray = 0,
                          FullSquareMatrix *celArray = 0);

     template<class Scalar>
       GenDBSparseMatrix<Scalar> *constructDBSparseMatrix(DofSetArray *dof_set_array=0,
                           Connectivity *cn=0);

     template<typename Scalar, typename SolverClass>
       GenEiSparseMatrix<Scalar,SolverClass> *constructEiSparseMatrix(DofSetArray *dof_set_array=0,
                           Connectivity *cn=0, bool flag=true);

     template<typename Scalar, typename SolverClass>
       GenEiSparseMatrix<Scalar,SolverClass> *constructGoldfarb(DofSetArray *dof_set_array=0,
                           Connectivity *cn=0);

     template<class Scalar>
       GenCuCSparse<Scalar> *constructCuCSparse(DofSetArray *dof_set_array=0);

     template<class Scalar>
       GenCuCSparse<Scalar> *constructCCSparse(DofSetArray *dof_set_array=0);

     template<class Scalar>
       GenSkyMatrix<Scalar> *constructSkyMatrix(DofSetArray*, Rbm *rbm = 0);

     template<class Scalar>
       GenBlockSky<Scalar> *constructBlockSky(DofSetArray *DSA);

     template<class Scalar>
       GenBLKSparseMatrix<Scalar> *constructBLKSparseMatrix(DofSetArray*, Rbm *rbm = 0);

     template<class Scalar>
       GenNBSparseMatrix<Scalar> *constructNBSparseMatrix();

     template<class Scalar>
       GenPCGSolver<Scalar, GenVector<Scalar>, GenSparseMatrix<Scalar> >
          *constructPCGSolver(GenSparseMatrix<Scalar> *K, Rbm *rbm=0);

     template<class Scalar>
       GenSpoolesSolver<Scalar> *constructSpooles(ConstrainedDSA *CDSA = 0);

     template<class Scalar>
       GenMumpsSolver<Scalar> *constructMumps(ConstrainedDSA *CDSA = 0, Rbm *rbm=0, FSCommunicator *com = 0);

     template<class Scalar>
       Rom::GenGalerkinProjectionSolver<Scalar> *constructGalerkinProjectionSolver();
     
     template<class Scalar>
      Rom::GenEiSparseGalerkinProjectionSolver<Scalar> *constructEiSparseGalerkinProjectionSolver();

     Rbm              *constructRbm(bool printFlag = true);
     Rbm              *constructHzem(bool printFlag = true);
     Rbm              *constructSlzem(bool printFlag = true);

     template<class Scalar>
       void addGravityForce(GenVector<Scalar>& force);

     template<class Scalar>
       void addPressureForce(GenVector<Scalar>& force, double lambda = 1.0);

     template<class Scalar>
       void addAtddnbForce(GenVector<Scalar>& force, double lambda = 1.0);

     template<class Scalar>
       void addAtdrobForce(GenVector<Scalar>& force, double lambda = 1.0);

     template<class Scalar>
       void addThermalForce(GenVector<Scalar>& force);

     template<class Scalar>
       void addMpcRhs(GenVector<Scalar>& force, double t = 0);

     void buildPrescDisp(Vector &v, double lambda);
     void buildPrescDisp(Vector &v, double t, double dt);

     template<class Scalar>
       void buildRHSForce(GenVector<Scalar> &force, GenSparseMatrix<Scalar> *kuc = 0);

     template<class Scalar>
       void computeReactionForce(GenVector<Scalar> &fc, GenVector<Scalar> &Vu,
                                 GenSparseMatrix<Scalar> *kuc, GenSparseMatrix<Scalar> *kcc = 0);
     void computeReactionForce(Vector &fc, Vector &Du, Vector &Vu, Vector &Au,
                               double *bcx, double *vcx, double *acx,
                               SparseMatrix *_kuc, SparseMatrix *_kcc,
                               SparseMatrix *_cuc, SparseMatrix *_ccc,
                               SparseMatrix *_muc, SparseMatrix *_mcc);
     void computeReactionForce(Vector &fc, GeomState *geomState, Corotator **corotators,
                               FullSquareMatrix *kel, double lambda, GeomState *refState);
     void computeReactionForce(Vector &fc, GeomState *geomState, Corotator **corotators,
                               FullSquareMatrix *kel, double time, GeomState *refState,
                               Vector &Vu, Vector &Au, double *vcx, double *acx,
                               SparseMatrix *_cuc, SparseMatrix *_ccc,
                               SparseMatrix *_muc, SparseMatrix *_mcc);
     bool reactionsReqd(double time, int step);

     template<class Scalar>
       void buildFreqSweepRHSForce(GenVector<Scalar> &force, GenSparseMatrix<Scalar> *muc,
                                   GenSparseMatrix<Scalar> **cuc_deriv, int iRHS, double omega);
     template<class Scalar>
       void buildRHSForce(GenVector<Scalar> &force,GenVector<Scalar> &tmp,
                         GenSparseMatrix<Scalar> *kuc,
                         GenSparseMatrix<Scalar> *muc,
                         GenSparseMatrix<Scalar> **cuc_deriv,
                         double omega, double delta_omega,
                         GeomState *gs=0);

     void initNodalTemperatures();
     double * getNodalTemperatures();

     // Main program control functions.
     void arcLength();
     void createCorotators(Corotator **allCorot);
     void preProcessing();
     FILE * openFile(char *fileName, const char *extension);
     void printStatistics();

     // static & freq response post processing function
     template<class Scalar>
     void postProcessing(GenVector<Scalar> &sol, Scalar *bcx, GenVector<Scalar> &force, 
                         int ndflag = 0, int index = 0, double time = 0, double eigV = 0.0,
                         GenSparseMatrix<Scalar> *kuc = NULL, GenSparseMatrix<Scalar> *kcc = NULL);

     void resProcessing(Vector &, int index=0, double t=0);

     // Nonlinear post processing function
     void postProcessing(GeomState *geomState, Vector &force, Vector &aeroForce, double time=0.0,
                         int step=0, double *velocity=0, double *vcx=0,
                         Corotator **allCorot=0, FullSquareMatrix *mArray=0, double *acceleration=0,
                         double *acx=0, GeomState *refState=0, Vector *reactions=0);

     // Pita Nonlinear post processing function
     void pitaPostProcessing(int timeSliceRank, GeomState *geomState, Vector &force, Vector &aeroForce, double time = 0.0,
                             int step = 0, double *velocity = 0, double *vcx = 0,
                             Corotator **allCorot = 0, FullSquareMatrix *mArray = 0, double * acceleration = 0);

     // Dynamic functions and thermal functions
     void aeroSolve();
     void aeroheatSolve();
     void thermoeSolve();
     void getHeatFlux(Vector &tsol, double *bcx, int fileNumber, int hgIndex,
                      double time=0);
     void getHeatFlux(ComplexVector &tsol, DComplex *bcx, int fileNumber, int hgIndex, double time=0)
       { cerr << " *** WARNING: Domain::getHeatFlux(Complex) is not implemented \n"; }
     void getTrussHeatFlux(Vector &tsol, double *bcx, int fileNumber,
                           int hgIndex, double time=0);
     void getTrussHeatFlux(ComplexVector &tsol, DComplex *bcx, int fileNumber, int hgIndex, double time=0)
       { cerr << " *** WARNING: Domain::getTrussHeatFlux(Complex) is not implemented \n"; }
     template <class Scalar>
       void computeConstantForce(GenVector<Scalar>& constantForce, GenSparseMatrix<Scalar>* kuc = 0);
     template <class Scalar> 
       void computeExtForce4(GenVector<Scalar>& force, const GenVector<Scalar>& constantForce, double t,
                             GenSparseMatrix<Scalar> *kuc = 0, ControlInterface *userSupFunc = 0,
                             GenSparseMatrix<Scalar> *cuc = 0, double tm = 0, GenSparseMatrix<Scalar> *muc = 0);
     template <class Scalar> 
       void computeExtForce(GenVector<Scalar>& force, double t, GenSparseMatrix<Scalar> *kuc = 0, ControlInterface *userSupFunc = 0,
                            GenSparseMatrix<Scalar> *cuc = 0, double tm = 0, GenSparseMatrix<Scalar> *muc = 0);
     void computeExtForce(Vector &f, double t, int tIndex,
                          SparseMatrix *kuc, Vector &prev_f);


     int  probType() { return sinfo.probType; }

     double computeStabilityTimeStep(DynamMat&);
     double computeStabilityTimeStep(FullSquareMatrix *kelArray, FullSquareMatrix *melArray, GeomState *geomState);

     void initDispVeloc(Vector& d_n, Vector& v_n, Vector& a_n, Vector &v_p, const char* = "");
     void initDispVelocOnTimeSlice (Vector& d_n, Vector& v_n, int sliceRank); // PITA: Use user-provided initial seeds
     void initTempVector(Vector& d_n, Vector& v_n, Vector& v_p);
     void writeRestartFile(double time, int timeIndex,
                           Vector &d_n, Vector &v_n, Vector &v_p, double Fref = 0.0, const char* = "");
     void writeRestartFile(double time, int timeIndex, Vector &v_n,
                           GeomState *geomState, const char* = "");
     void readRestartFile(Vector &d_n, Vector &v_n,
                          Vector &a_n, Vector &v_p, double *bcx,
                          double *vcx, GeomState &geomState, const char* = "");
     void getOrAddDofForPrint(bool ad, Vector& d_n, double* bcx, int iNode,
             double *xdata, int *dofx, double *ydata=0, int *dofy=0, double *zdata=0, int *dofz=0);
     void addVariationOfShape_StructOpt(int iNode, CoordSet *nodescopy, double &x, double &y, double &z);
     void aeroSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState = 0);
     void aeroheatSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState = 0);
     void thermohSend(Vector& d_n, Vector& v_n, Vector& a_n, Vector& v_p, double* bcx, double* vcx, GeomState* geomState = 0);
     void buildAeroelasticForce(Vector &f, PrevFrc& prevFrc, int tIndex, double t, double gamma, double alphaf, GeomState* geomState = 0);
     void buildAeroheatFlux(Vector &f, Vector &prev_f, int tIndex, double t);
     void thermoeComm();
     void dynamOutput(int, double, double*, DynamMat&, Vector&, Vector &, Vector&, Vector&, Vector&, Vector &, double*, double* = 0);
     void pitaDynamOutput(int, double* bcx, DynamMat&, Vector&, Vector &, Vector&, Vector&, Vector&, Vector &,
                          double* vcx, double* acx, int sliceRank, double time);

  protected:
     void dynamOutputImpl(int, double* bcx, DynamMat&, Vector&, Vector &, Vector&, Vector&, Vector&, Vector &, double*, double*, double, int, int);

  public:
     void tempdynamOutput(int, double*, DynamMat&, Vector&, Vector&, Vector&,
                          Vector&);

     double computeStructureMass(bool printFlag = true);
     double computeFluidMass();
     double getStructureMass();

     // Condition number estimated routines
     double computeConditionNumber(DynamMat&);

     //   Output Related functions
     template<class Scalar>
     int processOutput(OutputInfo::Type &type, GenVector<Scalar> &sol, Scalar *bcx, int iInfo,
                       double time, double freq = 0, int printFlag = 0);

     template<class Scalar>
     int processDispTypeOutputs(OutputInfo &oinfo, Scalar (*sol)[11], int iInfo,
                                 int numNodes, double time, double freq = 0, int printFlag = 0);
     void getStressStrain(Vector &sol, double *bcx, int fileNumber,
                          int strInd, double time = 0, int printFlag =0);
     void getStressStrain(ComplexVector &sol, DComplex *bcx, int fileNumber,
                          int strInd, double time = 0, int printFlag =0) { cerr << "Domain::getStressStrain is not implemented for complex\n"; }
     void getPrincipalStress(Vector &sol, double *bcx, int fileNumber,
                             int strInd, double time = 0);
     void getPrincipalStress(ComplexVector &sol, DComplex *bcx, int fileNumber,
                             int strInd, double time = 0) { cerr << "Domain::getPrincipalStress is not implemented for complex\n"; }
     void getElementForces(Vector &sol, double *bcx, int fileNumber,
                           int forceIndex, double time = 0);
     void getElementForces(ComplexVector &sol, DComplex *bcx, int fileNumber,
                           int forceIndex, double time = 0) { cerr << "Domain::getElementForces is not implemented for complex\n"; }
     void getElementAttr(int fileNumber, int typ, double time=0.0);
     void getElasticForces(Vector &dsp, double *bcx, Vector &ext_f, double eta,
                           FullSquareMatrix *kelArray=0);
     void getKtimesU(Vector &dsp, double *bcx, Vector &ext_f, double eta,
                     FullSquareMatrix *kelArray=0);
     //ADDED FOR SLOSHING PROBLEM, EC, 20070723
     void getSloshDispAll(Vector &tsol, double *bcx, int fileNumber, double time);
     void getSloshDispAll(ComplexVector &tsol, complex<double> *bcx, int fileNumber, double time) { cerr << "getSloshDispAll(complex) not implemented\n"; }
     //ADDED FOR SLOSHING PROBLEM, EC, 20070723
     void getSloshDisp(Vector &tsol, double *bcx, int fileNumber, int hgIndex, double time);
     void getSloshDisp(ComplexVector &tsol, complex<double> *bcx, int fileNumber, int hgIndex, double time) { cerr << "getSloshDisp(complex) not implemented\n"; }
     double getStrainEnergy(Vector &sol, double *bcx, SparseMatrix *gStiff=0,
                            FullSquareMatrix *kelArray=0);
     double getNodalStressStrain(Vector &sol ,double *bcx,
                                 int inode, int stressIndex, int surface);
     double getNodalDisp(Vector &sol ,double *bcx,
                         int node, int dispTyp);
     double getKineticEnergy(Vector &sol, SparseMatrix *gMass);
     double getDampingEnergy(Vector &sol, SparseMatrix *gDamp);

     template<class Scalar>
       void scaleDisp(Scalar *u);
     template<class Scalar>
       void scaleInvDisp(Scalar *u);

     template<class Scalar>
     int mergeDistributedDisp(Scalar (*xyz)[11], Scalar *u, Scalar *bcx = 0);//DofSet::max_known_nonL_dof

     Connectivity *makeSommerToNode();
     Connectivity *prepDirectMPC();
     // renumbering functions
     Renumber  getRenumbering();
     Renumber* getRenumberingFluid(); //ADDED FOR HEV PROBLEM, EC, 20070820

     void makeNodeToNode_sommer();

     // Eigen solver
     void eigenOutput(Vector& eigenValues, VectorSet& eigenVectors, double* bcx = 0, int convEig = 0); // modified for SLOSHING PROBLEM, EC, 20070723
     void setEigenValue(double _lbound, int _nshifts, int _maxArnItr = 0); //CBM
     void setEigenValues(double _lbound, double _ubound, int _neigps = 0, int _maxArnItr = 0);

     Solver *getSolver() { return solver; }

     template<class Scalar>
       void getSolverAndKuc(GenSolver<Scalar> *&solver, GenSparseMatrix<Scalar> *&kuc,
                            FullSquareMatrix *kelArray = 0, bool factorize=true);
     template<class Scalar>
       void getSolverAndKuc(AllOps<Scalar> &allOps, FullSquareMatrix *kelArray, bool factorize=true);

     void make_constrainedDSA();
     void make_constrainedDSA(int *bc);
     void make_constrainedDSA(int fake);

     // returns the number of dof
     int numdof()   { return dsa ? dsa->size() : -1; }

     // returns the number of unconstrained dof
     int numUncon() {
       return c_dsa ? c_dsa->size() : dsa ? dsa->size() : 0;
     }

     // returns the number of unconstrained Fluid dof
     //ADDED FOR HEV PROBLEM, EC, 20070820
     int numUnconFluid() {
       return c_dsaFluid ? c_dsaFluid->size() : dsaFluid->size();
     }
     // returns a pointer to the dirichlet boundary condtions
     BCond* getDBC() { return dbc; }

     // returns the number of dirichlet bc
     int  nDirichlet() { return numDirichlet; }
     int  nDirichletFluid() { return numDirichletFluid; } //ADDED FOR HEV PROBLEM, EC, 20070820
     int  nDispDirichlet() { return numDispDirichlet; }
     void setNumDispDirichlet(int n) { numDispDirichlet = n; }

     // returns the number of initial displacements
     int numInitDisp()  { return numIDis;  }
     int numInitDispModal() { return numIDisModal; }
     int numInitDisp6() { return numIDis6; }

     // returns a pointer to the initial displacement boundary condtions
     BCond* getInitDisp()  { return iDis;  }
     BCond* getInitDispModal() { return iDisModal; }
     BCond* getInitDisp6() { return iDis6; }

     // returns a pointer to the initial velocities
     int    numInitVelocity() { return numIVel; }
     BCond* getInitVelocity() { return iVel; }
     int    numInitVelocityModal() { return numIVelModal; }
     BCond* getInitVelocityModal() { return iVelModal; }

     // returns the number of neumann bc
     int  nNeumann() { return numNeuman; }

     // returns a pointer to the neumann boundary condtions
     BCond* getNBC() { return nbc; }

     // returns the number of nodes
     int  numNodes() { return (nodeToNode) ? nodeToNode->csize() : numnodes; }
     int  numNode()  { return numnodes; }
     void setNumNodes(int n)  { numnodes = n; }  // includes virtual nodes
     int  numGlobalNodes() { return numnodes; }

     // returns the number of elements
     int  numElements() { return numele; }
     void setNumElements(int n) { numele = n; }
     int addElem(int ele, int type, int nnd, int *nd)
        { packedEset.elemadd(ele, type, nnd, nd); return ele; }

     // returns the number of dofs
     void setNumDofs(int n) { numdofs = n; }
     int  numDofs() { return numdofs; }

     // returns the packed element set (allows element numbering gaps)
     Elemset& getElementSet() { return packedEset; }

     // returns the value of the gravity Force flag
     int  gravityFlag() { return gravityAcceleration ? 1: 0; }

     // returns the value of the pressure force flag
     int  pressureFlag();

     // returns the value of the contact force flag
     int  tdenforceFlag() { return int(nMortarCond > 0 && sinfo.newmarkBeta == 0.0 && sinfo.penalty == 0.0); } // TD enforcement (contact/tied surfaces with ACME) used for explicit dynamics

     int  thermalFlag() { return sinfo.thermalLoadFlag || sinfo.thermoeFlag >= 0; }

     // returns the maximum number of dofs per element
     int  maxNumDOF() { return maxNumDOFs; }

     // returns a pointer to the Connectivity allDOFs
     Connectivity *getAllDOFs() { return allDOFs; }
     void deleteAllDOFs() { delete allDOFs; allDOFs = 0; }

     CoordSet& getNodes() { return nodes; }

     ConstrainedDSA * getCDSA() { return c_dsa; }
     DofSetArray *    getDSA()  { return dsa;   }
     ConstrainedDSA * getCDSAFluid() { return c_dsaFluid; }
     DofSetArray *    getDSAFluid()  { return dsaFluid;   }

     // function that returns composite layer info
     LayInfo *getLayerInfo(int num);

     // function to get fluid exchanger
     FlExchanger * getFileExchanger() { return flExchanger; }
     void makeNodeTable(int topFlag);
     void makeTopFile(int topFlag);
     void makeAxiTopFile(int topFlag,int numSlices);

     template<class Scalar> friend class GenDecDomain;
     template<class Scalar> friend class GenDistrDomain;
     template<class Scalar> friend class GenSandiaDomain;
     friend class HData;

     // for output of tensor data

     double crossScale;

     double *** CPoint;
     double **  MidPoint;

     char *getBasename( char * fname);

     void getCompositeData(int iInfo,double time);

     void aeroPreProcess(Vector&, Vector&, Vector&, Vector &,double*,double*);
     void thermoePreProcess();

     void aeroHeatPreProcess(Vector&, Vector&, Vector&, double *bcx );
     void thermohPreProcess(Vector&, Vector&, Vector&, double *bcx );

     Connectivity *getNodeToNode();
     Connectivity *getNodeToElem() { return nodeToElem; }
     Connectivity* getNodeToNode_sommer() { return nodeToNode_sommer;}
     ConstrainedDSA *makeCDSA(int nbc, BCond *bcs);
     void getNormal2D(int node1, int node2, double &nx, double &ny);
     Elemset* getEset() { return &packedEset; }
     void setNumnodes(int n) { numnodes = n; }

     void addNode(int nd, double *xyz) { nodes.nodeadd(nd,xyz); }

     // PJSA: mpc stuff
  protected:
     bool haveNodes;
  public:
     int getNumLMPC() { return numLMPC; }
     void setNumLMPC(int n) { numLMPC = n; }
     LMPCons* getLMPC(int i) { return lmpc[i]; }
     LMPCons* getFsi(int i) { return fsi[i]; }
     ResizeArray<LMPCons *>* getLMPC() { return &lmpc; }
     int addNodalCTC(int n1, int n2, double nx, double ny, double nz,
                     double normalGap = 0.0, int _mode = -1, int lagrangeMult = -1, double penalty = 0.0);
     int getNumCTC() { return numCTC; }
     void addNodeToNodeLMPCs(int lmpcnum, int n1, int n2, double face_normal[3], double gap_vector[3], int itype);
     void addDirichletLMPCs(int _numDirichlet, BCond *_dbc);
     void deleteAllLMPCs();
     void deleteSomeLMPCs(mpc::ConstraintSource s);
     void UpdateContactSurfaceElements(GeomState *);

     // HB: mortar stuff (EXPERIMENTAL)
  protected:
     int nSurfEntity;                          // I know, this should be in GeoSource !!!
     ResizeArray<SurfaceEntity*> SurfEntities; //
     int nMortarLMPCs;                         // total number of Mortar LMPCs generated
     Connectivity* mortarToMPC;                //
     vector<int> contactSurfElems;
     std::set<int> aeroEmbeddedSurfaceId;  //KW: Ids of wet surfaces
  public:
     int AddSurfaceEntity(SurfaceEntity*);
     int AddSurfaceEntity(SurfaceEntity*, int isurf);
     void PrintSurfaceEntities();

     int AddAeroEmbedSurfaceId(int Id);
     std::set<int> & GetAeroEmbedSurfaceId() { return aeroEmbeddedSurfaceId; }

     int nMortarCond;
     int nContactSurfacePairs; 
     ResizeArray<MortarHandler*> MortarConds;

     int AddMortarCond(MortarHandler*);
     void PrintMortarConds();
     void DeleteMortarConds();

     void SetMortarPairing();
     void SetUpSurfaces(CoordSet* cs = 0);
     void UpdateSurfaces(GeomState *, int config_type = 1);
     void UpdateSurfaces(DistrGeomState *geomState, int config_type, SubDomain **sd);
     void MakeNodalMass(SparseMatrix *M, SparseMatrix *Mcc);
     void MakeNodalMass(SubDOp *M, SubDomain **sd);

     void InitializeDynamicContactSearch(int numSub = 0, SubDomain **sd = 0);
     void RemoveGap(Vector &g);
     void PerformDynamicContactSearch(double dt_old, double dt);
     void AddContactForces(double dt_old, double dt, Vector &f);
     void AddContactForces(double dt_old, double dt, DistrVector &f);

     void InitializeStaticContactSearch(MortarHandler::Interaction_Type t, int numSub = 0, SubDomain **sd = 0);
     void UpdateSurfaces(MortarHandler::Interaction_Type t, GeomState *geomState);
     void UpdateSurfaces(MortarHandler::Interaction_Type t, DistrGeomState *geomState, SubDomain **sd);
     void PerformStaticContactSearch(MortarHandler::Interaction_Type t);
     void ExpComputeMortarLMPC(MortarHandler::Interaction_Type t, int nDofs = 0, int* Dofs = 0);

     void ComputeMortarLMPC(int nDofs = 0, int* Dofs = 0);
     void computeMatchingWetInterfaceLMPC();

     void CreateMortarToMPC();

     Connectivity* GetMortarToMPC();
     int GetnMortarConds() { return nMortarCond; }
     int GetnContactSurfacePairs() { return nContactSurfacePairs; }
     int GetnMortarLMPCs();
     MortarHandler* GetMortarCond(int i) { return MortarConds[i]; }
     ResizeArray<SurfaceEntity*>* viewSurfEntities() { return(&SurfEntities); }
     void setNumSurfs(int nSurfs) { nSurfEntity = nSurfs; }
     int getNumSurfs() { return(nSurfEntity); }

#ifdef HB_ACME_FFI_DEBUG
     void WriteFFITopFile(FILE* file);
#endif
  protected:
     void initialize();

     // FETI-DPH acoustics
     template<class Scalar>
       void assembleSommer(GenSparseMatrix<Scalar> *K, AllOps<Scalar> *ops = 0);
     template<class Scalar>
       void computeSommerDerivatives(double HH, double KK, int curvatureFlag, int *dofs, FullSquareMatrix &ms,
                                     DComplex **bt2nMatrix, double kappa, double ss, AllOps<Scalar> *ops);
     template<class Scalar>
       void assembleATDROB(GenSparseMatrix<Scalar> *K, AllOps<Scalar> *ops = 0, double Kcoef = 0.0);
     template<class Scalar>
       void updateMatrices(AllOps<Scalar> *ops, GenSparseMatrix<Scalar> *K, int *dofs,
                           FullSquareMatrix *reEl, FullSquareMatrix *imEl,double Kcoef = 0.0);
     template<class Scalar>
       void updateDampingMatrices(AllOps<Scalar> *ops, int *dofs, FullSquareMatrix *reEl,
                                  FullSquareMatrix *imEl, double ss, int n);
     struct WetInterface {
       int fluidSurfaceID;
       int structureSurfaceID;
     };
     ResizeArray<WetInterface *> *wetInterfaces;
     int nWetInterface;
     bool firstOutput;

  public:
     int isFluidElement(int i);
     int isStructureElement(int i);
     bool isComplex() {
       return ((numComplexDirichlet > 0) || (numComplexNeuman > 0)
               || (numSommer > 0) || (numComplexLMPC > 0) ||
                   sinfo.hasDamping() || PMLFlag || packedEset.hasDamping() );
     }
     bool isHomogeneous();
     void addWetInterface(int fluidSurfaceID, int structureSurfaceID, double normal_tol = 0.1, double tangential_tol = 0.001);
     int* getAllWetInterfaceNodes(int &count);
     double getFrequencyOrWavenumber();
     void computeCoupledScaleFactors();
     void getInterestingDofs(DofSet &ret, int glNode);

     double** getCMatrix(); //ADDED FOR HEV PROBLEM, EC, 20070820
     void multC(const Vector&, Vector&);
     void trMultC(const Vector&, Vector&);

     ControlLawInfo* getClaw() { return  claw;}

     virtual double densProjCoeff(int dof) { return 1.0; }
     virtual void densProjectStiffness(GenFullSquareMatrix<double>& kel, int num) { /* do nothing */ }
     virtual void densProjectStiffnessC(GenFullSquareMatrix<DComplex>& kel, int num) { /* do nothing */ }

     int getMaxNumNodes() { return maxNumNodes; }

     void initSfem();

     ConstrainedDSA *makeMaps(DofSetArray *dsa, ConstrainedDSA *cdsa, DOFMap *baseMap, DOFMap *eqMap);

     /** Replaces all 6-DOF rigid elements that share a node by a single element */
     void collapseRigid6();

     // new controls
     void updateUsddInDbc(double* userDefineDisp, int* map = 0);
     void updateUsdfInNbc(double* userDefineForce, int* map = 0, double* weight = 0);
     void updateActuatorsInNbc(double* actuatorsForce, int* map = 0, double* weight = 0);

     int outFlag; // if this is set to 1 then the output file should use the compressed numbering (i.e. with gaps removed)
                  // compatible with top files generated using -T command line argument
                  // the default is 0, for compatibility with top files generated using -t command line argument
     int *nodeTable;
     int exactNumNodes;
};

#ifdef _TEMPLATE_FIX_
  #include <Driver.d/OpMake.C>
#endif

#endif
