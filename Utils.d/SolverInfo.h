#ifndef _SOLVER_INFO_
#define _SOLVER_INFO_

#include <Utils.d/NonlinearInfo.h>
#include <Utils.d/SensitivityInfo.h>
#include <Feti.d/FetiInfo.h>
#include <Utils.d/Conwep.d/BlastLoading.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <algorithm>
#include <list>
#include <vector>
#include <limits>
#include <fstream>
#include <set>

#if defined(WINDOWS) || defined(MACOSX)
 #include <cfloat>
#else
 #include <limits>
#endif

const double defaultTemp = -10000000.0;

struct AdaptiveSweepParams {
public:
     int maxP,minRHS,maxRHS,deltaRHS,numS;
     double w1,w2,atol;
     bool dgp_flag;
};

struct SweepParams {
public:
   SweepParams() { 
                  nFreqSweepRHS = 8;
                  freqSweepMethod = Taylor;
                  isAdaptSweep = false;
                  padeL = 9;
                  padeM = 10;
                  padeN = 2;
                  pade_pivot = false;
                  pade_tol = 1.0e-16;
                  pade_poles = false;
                  pade_poles_sigmaL = 0.0;
                  pade_poles_sigmaU = std::numeric_limits<double>::max();
                  alphaD = 0.0; betaD = 0.0;
   }
   int nFreqSweepRHS;
   enum { Taylor, Pade1, Pade, Fourier, PadeLanczos, GalProjection, KrylovGalProjection, QRGalProjection };
   AdaptiveSweepParams adaptSweep;
   bool isAdaptSweep;
   int freqSweepMethod;
   int padeL, padeM, padeN;
   bool pade_pivot;
   double pade_tol;
   bool pade_poles;
   double pade_poles_sigmaL, pade_poles_sigmaU;
// RT: Does not belong here but added to allow different damping factors for
// different impedance sections
   double alphaD;
   double betaD;
};

struct SolverInfo {

 private:
   double dt;           // time step value
   double dtemp;        // thermal time step value

 public:
   // Problem Type parameters
   enum { Static, Dynamic, Modal, NonLinStatic, NonLinDynam, 
          ArcLength, ConditionNumber, TempDynamic, Top,
          AxiHelm, MatNonLinStatic, MatNonLinDynam,
          Helmholtz, HelmholtzFreqSweep, HelmholtzDirSweep, HelmholtzMF, HelmholtzSO,
          Decomp, NonLinTempDynam, DisEnrM, PodRomOffline,
          None }; // note to developers: if you add a new entry in ths enum then
                  // you should also modify problemTypeMessage in Driver.d/Static.C
   enum SensitivityMethod { Direct, Adjoint }; 

   SensitivityMethod sensitivityMethod;
   int probType;
   int soltyp; // from CONTROL statement: 1 = statics, 2 = heat conduction, etc...

   float ATDARBFlag;
   float ATDDNBVal;
   float ATDROBVal;
   float ATDROBalpha;
   float ATDROBbeta;
   int aeroFlag;
   bool dyna3d_compat;
   int aeroheatFlag;
   int thermoeFlag;
   int thermohFlag;
   int modeDecompFlag;
   int isCollocated;
   double alphas[2];
   double alphat[2];

   // Added by Alex Main.  A predictor on the velocity of the structure
   double alphasv;

   int structoptFlag;
   int thermalLoadFlag; // whether there is a thermal load applied
   int radiationFlag; // whether there are radiation lhs/rhs terms for linear analysis
   int hzemFlag; // Zero energy mode (for thermal Problems),
                 // Equivalent to grbm for structures
   int hzemFilterFlag;
   int slzemFlag; // Flag for zero energy mode
   int slzemFilterFlag; // Flag for zero energy mode filter

   double mppFactor;  // modal amplification factor for mpp command

   // Solver parameters
   int type;     // 0 = direct, 1 = iterative, 2 = FETI, 3 = Block Diag
   int subtype;  // subtype ... 9 is mumps  10 is diag
   int iterType; // 0 = CG, 1 = GMRES, 2 = GCR, 4 = BCG, 5 = CR
   int iterSubtype; // matrix storage
   int precond;  // preconditioner 0 = none, 1 = jacobi
   int maxit;    // maximum number of iterations
   double tol;   // tolerance for convergence
   int renum;    // renumbering scheme for skyline: 0 = none (default), 1 = sloan, 2 = RCM
   int maxvecsize;  // for pcg # of krylov vectors to store default = 0
   bool lastIt;

   // Time Integration Algorithms
   enum { Newmark, Qstatic };

   // Parameters for PITA
   bool activatePita;                // Enables PITA (both linear and nonlinear problems)
   bool mdPita;                // Enables time-space parallelism (linear problem only)
   bool pitaNoForce;           // Enables NoForce optimization (linear problem only)
   int pitaTimeGridRatio;      // Coarse/fine time-grid ratio (always required)  
   int pitaMainIterMax;        // Maximum number of main iterations (always required)
   int pitaProcessWorkloadMax; // Maximum number of active time-slices on each CPU (always required)
   //int numSpaceMPIProc;      // Number of CPUs in the space domain (only when time-space parallelism is enabled)
   int pitaGlobalBasisImprovement;    // 1 (default) = Seeds only, 2 = Seeds and propagated seeds (nonlinear problem only)
   int pitaLocalBasisImprovement;     // 0 (default) = Global only, 1 = Local increments only (nonlinear problem only)
   bool pitaRemoteCoarse;      // Coarse grid integrator on dedicated CPU 
   double pitaProjTol;         // Tolerance used to build the projector
   bool pitaTimeReversible;    // true if PITA exploits the time-reversibility of the problem
   bool pitaReadInitSeed;      // true if PITA uses provided initial seed information
   double pitaJumpCvgRatio;    // Accumulated jump-based convergence criterion
                               // (0.0 = Use default value: pitaTimeGridRatio^2, -1.0 = Deactivated)
   bool pitaJumpMagnOutput;    // Enables the output of the relative magnitude of the jumps

   bool modal;          // true iff system is to be solved in modal coordinates
   bool acoustic;       // true iff system is to be solved for acoustic time domain
   bool modifiedWaveEquation; // true if solving using the modified wave equation
   double modifiedWaveEquationCoef; // value for the coefficient (theoretically 12.0)
   int order;           // used for dynamics... 1: first order (heat), 2: second order (mech, acoustics)
   int timeIntegration; // type of time integration scheme
   int initialTimeIndex;// initial time index (either 0 or from restart)
   double initialTime;  // initial time (either 0.0 or from restart)
   double initExtForceNorm;  // initial Force Norm (for restarting qstatics)
   double tmax;         // maximum time
   double t_AeroF;
   bool stop_AeroF, stop_AeroS;
/* these are now private, use getTimeStep and setTimeStep functions to access
   double dt;           // time step value
   double dtemp;        // thermal time step value
*/
   double alphaDamp;    // Rayleigh Mass damping coefficient 
   double betaDamp;     // Rayleigh Stiffness damping coefficient
   double etaDamp;      // Structural stiffness damping coefficient
   int mtypeDamp;       // type of damping moments (0 = axial, 2 = follower)
   double alphaTemp;
   double newmarkBeta;  // Newmark algorithm parameter (beta)
   double newmarkGamma; // Newmark algorithm parameter (gamma)
   double newmarkAlphaF;  // Newmark algorithm parameter (alphaf)
   double newmarkAlphaM; // Newmark algorithm parameter (alpham)
   int stable;          // 0: do not compute the stability timestep for explicit dynamics
                        // 1: compute stability timestep for explicit dynamics but only use the computed value if 
                        //    it is less than the previously computed (or initially specified) value
                        // 2: compute stability timestep for explicit and use it regardless of the previously computed
                        //    or initially specified value (i.e. the timestep can either increase or decrease)
   double stable_tol;   // convergence tolerance for computionation of stability timestep
   int stable_maxit;    // stable_maxit*n is the maximum number of iterations for computation of stability timestep
   double stable_cfl;   // factor by which computed stability time step is multiplied
   int stable_freq;     // number of timesteps between stability timestep updates (nonlinear only)
   bool check_energy_balance;
   double epsilon1;
   double epsilon2;
   double qsMaxvel;     // final relaxation parameter in quasi-static alg.
   double delta;        // translation factor from time index to real time for quasistatics
   int nRestart;        // how many time steps per restart
   int steadyFlag;      // quasi-transient steady state simulation
   int steadyMin;       // minimum number of steps in steady state simulation
   int steadyMax;       // maximum number of steps in steady state simulation
   double steadyTol;    // steady state criteria
   double sensitivityTol;    // Sensitivity analysis tolerance
   double ratioSensitivityTol;    // Sensitivity analysis tolerance
   double qsBeta;       // relaxation parameter for computing alphaR in case of
                        // thermal quasistatic
   bool no_secondary;
   // Eigenvalue Problem parameters
   enum {SubSpace, LobPcg, Arpack};

   int eigenSolverType; // type of eigensolver
   int nEig;            // number of eigen pairs to compute 
   int eigenSolverSubType; // subtype of eigensolver
   
   int subspaceSize;    // number of vectors in the subspace
   double tolEig;       // subspace iteration tolerance
   double tolJac;       // jacobi iteration tolerance
   bool   explicitK;    // determines whether or not to explicitly form K during eigprob
   int maxitEig;       
   bool qrfactorization; // determines if qr factorization of eigen modes to be conducted or not   
                         // if true, qmatrix and rmatrix in OUTPUT has to be specified.
                         // this is mainly used for interpolation of ROM.
   const char* xmatrixname;
   const char* qmatrixname;
   const char* rmatrixname;
   const char* eigenvaluename;

   // Sloshing problem flag
   int sloshing;
   // Hydroelastic vibration problem flag
   int HEV;
 
   const char *which;   /* for ARPACK: specifies "which" of the eigenvalues to compute, where:
                          "LA" - compute the eigenvalues just to the right of sigma (sigma is the shift)
                          "SA" - compute the eigenvalues just to the left of sigma
                          "LM" - undefined for shift-invert on generalized eigenvalue problem
                          "SM" - undefined for shift-invert on generalzied eigenvalue problem
                          "BE" - compute eigenvalues to either side of sigma */
   int arpack_mode;     // 3 is shift-invert, 4 is buckling mode
   double lbound;       // lower bound of a set or range of eigenvalues to be computed by ARPACK
   double ubound;       // upper bound of a set or range of eigenvalues to be computed by ARPACK
   int nshifts;         // number of shifts to be used by ARPACK when computing neigenpa eigenvalues
                        // that are greater then lbound
   int neigps;          // number of consecutive eigenvalues computed within the range [lbound,ubound]
                        // for every shift starting with lbound
   bool filtereig;      // filtering very small eigenvalues (if necessary)
   int maxArnItr;       // EXPERIMENTAL: specifies the maximum number of Arnoldi iterations (in ARPACK)


   // Rigid Body Mode parameters
   double trbm;         // algebraic rbm tolerance 
   double trbm2;        // algebraic rbm tolerance used for sparse/skyline when GRBM is activated
   double tolsvd;       // singular value decomposition tolerance
   int rbmflg;          // 0 = algebraic rbm, 1 = geometric rbm
   std::set<int> rbmFilters; // selective rbm filtering for linear dynamics problems
   bool grbm_use_lmpc;  // true (default) = lmpcs treated algebraically in GRBM method
                        // false = lmpcs are assumed to not introduce any mechanisms (i.e. like beams)
   std::vector<double> grbm_ref; // coordinates of reference point for rotational modes; if empty then 1st node is used.

   double condNumTolerance; // Condition number tolerance
   int condNumMaxit;

   int massFlag;
   int filterFlags;
   int filterQ; // 0 = Q=M for statics/quasistatics, 1 = Q=I for statics/quasistatics

   int zeroInitialDisp; // flag to set initial disp to zero

   bool pivot;  // true if pivoting is to be used in spooles solver
   int spooles_scale; // true if scaling is to be used in spooles solver
   double spooles_tau;  // used when pivoting is enabled, all entries in L and U have magnitude
                        // less than or equal to tau, default is 100.
   double spooles_maxzeros; // see Solvers.d/Spooles.C for description
   int spooles_maxsize, spooles_maxdomainsize, spooles_seed, spooles_msglvl; // see Solvers.d/Spooles.C for description
   int spooles_renum; // renumbering scheme for spooles: 0 = best of ND and MS, 1 = MMD, 2 = MS, 3 = ND
   int sparse_renum;  // renumbering scheme for BLKSparseMatrix: 0 = esmond MMD (default), 1 = metis ND
   int sparse_maxsup, sparse_defblk;

   double goldfarb_tol;

   // YC : sensitivity and optimization parameters
   bool sensitivity;

   // map object for Mumps control CNTL and ICNTL matrices
   std::map<int, int> mumps_icntl;
   std::map<int, double> mumps_cntl;

   bool localScaled, coarseScaled;

   int curSweepParam;
   std::map<int,SweepParams> sweepParams;
   SweepParams* getSweepParams() { return &(sweepParams[curSweepParam]); }
   bool doFreqSweep,doEigSweep;
/*
   int nFreqSweepRHS;
   enum { Taylor, Pade1, Pade, Fourier, PadeLanczos, GalProjection, KrylovGalProjection, QRGalProjection };
   AdaptiveSweepParams adaptSweep;
   bool isAdaptSweep;
   int freqSweepMethod;
   int padeL, padeM, padeN;
   bool pade_pivot;
   double pade_tol;
   bool pade_poles;
   double pade_poles_sigmaL, pade_poles_sigmaU;
*/

   bool test_ulrich;
   int modeFilterFlag;


   int addedMass;
   bool isCoupled;
   double coupled_scale;
   bool isMatching;  // true if only one wet interface is given, otherwise false (default) 
   bool farfield; // true if farfield output requested

   bool dbccheck;
   int contact_mode;
   int contactsurface_mode;
   bool trivial_detection;

   bool noninpc;
   bool inpc;
   int nsample;

   std::map<int, int> debug_icntl;   // used for debugging
   std::map<int, float> debug_cntl;  // used for debugging

   bool iacc_switch; // mech/acou: true --> compute consistent initial second time derivative ie, a^0 = M^{-1}(fext^0 - fint^0 - Cv^0) for a second order differential equation (ie mech/acou)
                     //            false --> a^0 = 0
                     // heat: true --> compute consistent initial first time derivative ie, v^0 = M^{-1}(fext^0 - fint^0) for a first order differential equation (ie heat)
                     //       false --> v^0 = 0
   bool zeroRot;
   

   int dist_acme; // 0: sequential, 1: parallel with centralized input on host (cpu with id 0), 2: parallel with distributed input by subdomain
                  // NOTE: currently only dist_acme == 0 is supported for Mortar method (statics and implicit dynamics) ... see main.C 
                  //       dist_acme == 2 requires a special decomposition, master and slave must be in the same subdomain
   bool allproc_acme; // true: use all available processors, false: use only processors with subdomain/s having surface interactions
   bool ffi_debug;
   double mortar_scaling;
   int mortar_integration_rule;
   bool tdenforceFlag;

   bool lagrangeMult;
   double penalty;
   int mpcDirect;
   bool usePrescribedThreshold;
   double mpcDirectTol; // threshold for definition of a null pivot is defined as mpcDirectTol*epsilon
   double coefFilterTol, rhsZeroTol, inconsistentTol;
   int constraint_hess;
   double constraint_hess_eps;

   // options used for the augmented Lagrangian constraint enforcement method
   int num_penalty_its;
   double penalty_tol;
   double penalty_beta;
   bool reinit_lm; // true : the Lagrange multipliers are re-initialized to zero at the start of each load/time-step
                   // false (default) : the Lagrange multipliers are initialized to zero at the start of the first load/time-step
                   //                   subsequently, their values are retained.
   int lm_update_flag; // 1 (default) : the Lagrange multipliers are updated after each Newton solve in the penalty iteration loop.
                       // 0 : the Lagrange multipliers are updated after each Newton solve except for the final one (i.e. the one
                       //     done before terminating the penalty iteration loop due to maximum number of iterations or satisfaction
                       //     of the minimum constraint violation condition).

   std::vector<std::string> RODConversionFiles;
   std::vector<std::string> PODerrornorm;
   int numSnap;
   std::vector<std::string> snapfiPodRom;
   std::vector<std::string> robfi;
   std::vector<double> snapshotWeights;
   std::vector<std::string> readInROBorModes;
   std::vector<std::string> readInDualROB; 
   std::map<std::pair<int,int>,std::string> readInLocalBasesAuxi;
   std::vector<std::string> readInLocalBasesCent;
   const char * readInModes;
   const char * SVDoutput;
   const char * reducedMeshFile;
   const char * readInShapeSen;
   std::vector<BlastLoading::BlastData> conwepConfigurations;
   std::vector<std::string> statePodRomFile;
   std::vector<std::string> velocPodRomFile;
   std::vector<std::string> accelPodRomFile;
   std::vector<std::string> dsvPodRomFile;
   const char * isvPodRomFile;
//   const char * dsvPodRomFile;
   const char * forcePodRomFile;
   const char * constraintPodRomFile;
   const char * constraintSnapshotFile;
   const char * constraintViolationFile;
   const char * residualPodRomFile;
   const char * jacobianPodRomFile;
   bool ROMPostProcess;
   bool statevectPodRom;
   bool velocvectPodRom;
   bool accelvectPodRom;
   bool isvPodRom;
   bool dsvPodRom;
   bool forcevectPodRom;
   bool residvectPodRom;
   bool jacobvectPodRom;
   bool readmodeCalled;
   bool modalCalled;
   bool modalLMPC;
   bool readShapeSen;
   double ksParameter;
   double ksMax;
   bool activatePodRom;
   bool snapshotsPodRom;
   bool checkPodRom;
   bool svdPodRom;
   int  svdBlockSize;
   bool clusterSubspaceAngle;
   int clustering;
   int  solverTypeCluster; // 0: Random, 1: K-means 2: Sparse Supspace Clustering
   int use_nmf;
   int nmfNumROBDim;
   int nmfDelROBDim;
   int nmfRandInit;
   int nmfMaxIter;
   int nmfNumSub;
   double nmfTol;
   int nmfPqnNumInnerIter;
   double nmfPqnAlpha;
   bool DEIMBasisPod;
   bool UDEIMBasisPod;
   bool ConstraintBasisPod;
   bool ReducedStiffness;
   bool computeForceSnap;
   bool computeConstraintSnap;
   bool filterSnapshotRows;
   bool orthogForceSnap;
   bool orthogConstraintSnap;
   bool computeDEIMIndices;
   bool DEIMPodRom;
   bool UDEIMPodRom;
   bool samplingPodRom;
   bool snapProjPodRom;
   bool galerkinPodRom;
   bool elemLumpPodRom;
   bool onlineSvdPodRom;
   int  maxSizePodRom;
   double romEnergy;
   std::vector<int> localBasisSize;
   std::vector<int> localDualBasisSize;
   int  maxSizeDualBasis;
   int  maxDeimBasisSize;
   bool selectFullNode;
   bool selectFullElem;
   int  forcePodSize;
   int  constraintPodSize;
   int  normalize;
   bool subtractRefPodRom;
   bool useScalingSpnnls;
   bool useCenterSpnnls;
   bool useReverseOrder;
   bool projectSolution;
   bool positiveElements;
   int  solverTypeSpnnls; // 0: Lawson & Hanson, 1: Conjugate Gradient Pursuit
   double maxSizeSpnnls;
   int maxElemSpnnls;
   double maxIterSpnnls;
   bool reduceFollower;
   bool randomVecSampling;
   int  skipPodRom;
   int  skipOffSet;
   int  skipState;
   int  skipVeloc;
   int  skipAccel;
   int  skipInternalStateVar;
   int  skipDualStateVar;
   int  skipForce;
   int  skipResidual;
   int  skipJacobian;
   int  orthogPodRom;
   int  numRODFile;
   int  romresidType;
   double tolPodRom;
   bool useMassNormalizedBasis;
   bool useMassOrthogonalProjection;
   bool performMassNormalization;
   bool ConwepOnOff;
   std::list<int> loadcases;
   bool basicDofCoords; // if this is true then all of the nodes use the basic coordinate frame 0 for DOF_FRM
   bool basicPosCoords; // if this is true then all of the nodes use the basic coordinate frame 0 for POS_FRM
   bool scalePosCoords;
   double xScaleFactor, yScaleFactor, zScaleFactor;
   double xLMPCFactor, yLMPCFactor, zLMPCFactor;
   std::vector<double> xScaleFactors, yScaleFactors, zScaleFactors;
   int inertiaLumping; // 1: diagonal lumping (default), 2: block-diagonal 6x6 lumping
                       // note #1: this flag is automatically set to 2 when a product of inertia is defined using DIMASS
                       //          or when a discrete mass element (type 131) is defined.
                       // note #2: for diagonal lumping, consistent mass matrices are diagonally lumped and
                       //          off-diagonal entries of lumped mass matrices which are block-diagonal are dropped.
                       // note #3: for block-diagonal lumping, consistent mass matrices are still diagonally lumped
                       //          but lumped mass matrices which are block-diagonal remain so.
   bool printMatLab;
   const char * printMatLabFile;

   bool elementDeletion;
   std::map<int,double> deleteElements; // elements to be deleted at specific time or times (specified in input file)
   bool piecewise_contact;
   bool piecewise;
   bool freeplay;
   double piecewise_dlambda, piecewise_maxLambda;

   int npMax;   // Max number of elements in the reduced mesh for the ScalaPack LH parse solver.
   int scpkMB;  // Scalapack row block size
   int scpkNB;  // Scalapack column block size
   int scpkMP;  // Scalapack row processor grid size
   int scpkNP;  // Scalapack column processor grid size

   bool useMassAugmentation; // adjust lumped mass matrix for shell element type 15 to increase stability time-step

   // Constructor
   SolverInfo() { filterFlags = 0;
                  filterQ = 1;
                  type = 0;     
                  soltyp = -1;
                  subtype = 0; // By default we use direct Skyline
                  renum = 0;
                  probType = SolverInfo::None;
                  alphaDamp = 0.0;
                  betaDamp = 0.0;
                  etaDamp = 0.0;
                  mtypeDamp = 2;
                  modal = false;
                  lastIt = false;
                  mppFactor = 1.0;
                  sensitivityMethod = SolverInfo::Direct;
 
                  // Parameters for sensitivity
                  sensitivity = false;
                 
                  // Parameters for PITA 
                  activatePita = false;
                  mdPita = false;
                  pitaNoForce = false;
                  pitaTimeGridRatio = 1;
                  pitaMainIterMax = 0;
                  pitaProcessWorkloadMax = 1;
                  //numSpaceMPIProc = 1;
                  pitaGlobalBasisImprovement = 1;
                  pitaLocalBasisImprovement = 0;
                  pitaRemoteCoarse = false;
                  pitaProjTol = 1.0e-6;
                  pitaTimeReversible = false;
                  pitaReadInitSeed = false;
                  pitaJumpCvgRatio = -1.0;
                  pitaJumpMagnOutput = false;

                  acoustic = false;
                  modifiedWaveEquation = false;
                  modifiedWaveEquationCoef = 12.0;

                  tmax = 0.0; 
                  t_AeroF = 0.0;
                  stop_AeroF = false;
                  stop_AeroS = false;
                  dt = 0.0;
                  dtemp = 0.0;

                  steadyFlag = 0;
                  steadyMin = 1;
                  steadyMax = 10;
                  steadyTol = 1.0e-3; 
                  sensitivityTol = 1.0e-5; 
                  ratioSensitivityTol = 1.0; 
                  qsMaxvel = 1.0;
                  qsBeta = 1.0;
                  delta = 0.0;
                  no_secondary = false;

                  trbm = 1.0E-16;   // default zero pivot tolerance
                  trbm2 = 1.0E-16;  // default zero pivot tolerance
                  tolsvd = 1.0E-6;  // default singular value tolerance
                  massFlag = 0;     // whether to calculate total structure mass
                                  
                  ATDARBFlag = -2.0;
                  ATDDNBVal = 0.0;
                  ATDROBVal = 0.0;
                  ATDROBalpha = 0.0; //this value can not be 0 when Robin boundary is set, it is the flag!
                  ATDROBbeta = 0.0;

                  aeroFlag = -1;
                  dyna3d_compat = false;
                  aeroheatFlag = -1;
                  thermoeFlag = -1;
                  thermohFlag = -1;
                  thermalLoadFlag = 0;
                  radiationFlag = 0;
                  modeDecompFlag = 0;
                  hzemFlag = 0;
                  hzemFilterFlag = 0;
                  slzemFlag = 0;
                  slzemFilterFlag = 0;
                  condNumTolerance = 1.0E-3;
                  condNumMaxit = 100;
                  nRestart = -1;    // default do not write binary restart file
                  initialTimeIndex = 0;
                  initialTime = 0.0;
                  initExtForceNorm = 0.0;
                  zeroInitialDisp = 0;
                  gepsFlg = 0;
                  alphas[0] = alphas[1] = 0.0;
                  structoptFlag = 0;
                  // 2nd order newmark default: constant average acceleration
                  order = 2; timeIntegration = Newmark; newmarkBeta = 0.25; newmarkGamma = 0.5; newmarkAlphaF = 0.5; newmarkAlphaM = 0.5;
                  stable = 1;
                  stable_tol = 1.0e-3;
                  stable_maxit = 100;
                  stable_cfl = 0.8;
                  stable_freq = 1000;
                  check_energy_balance = false;
                  epsilon1 = 0.01;
                  epsilon2 = 1e-10;
                  rbmflg = 0;
                  grbm_use_lmpc = true;
                  buckling = 0;

                  sparse_renum = 0;
                  sparse_maxsup = 100;
                  sparse_defblk = 10;
                  pivot = false;
                  spooles_scale = 0;
                  spooles_tau = 100.;
                  spooles_seed = 532196;
                  spooles_maxsize = 64;
                  spooles_maxdomainsize = 24;
                  spooles_maxzeros = 0.04;
                  spooles_msglvl = 0;
                  spooles_renum = 0;
                  goldfarb_tol = 1.0;
                  explicitK = false;
                  qrfactorization = false;
                  eigenvaluename = "";
                  xmatrixname = "";
                  qmatrixname = "";
                  rmatrixname = "";
                  localScaled = false;
                  coarseScaled = false;

                  doFreqSweep = false;
                  doEigSweep = false;
/*
                  nFreqSweepRHS = 8;
                  freqSweepMethod = Taylor;
                  isAdaptSweep = false;
                  padeL = 9;
                  padeM = 10;
                  padeN = 2;
                  pade_pivot = false;
                  pade_tol = 1.0e-16;
                  pade_poles = false;
                  pade_poles_sigmaL = 0.0; pade_poles_sigmaU = std::numeric_limits<double>::max();
*/
                  modeFilterFlag = 0;
                  test_ulrich = false;
                  addedMass = 1;
                  isCoupled = false;
                  coupled_scale = 1.0;
                  isMatching = false; 
                  farfield = false;

                  dbccheck = true;
                  contact_mode = 1;
                  contactsurface_mode = 1;
                  trivial_detection = false;

                  nEig = 0;
                  eigenSolverSubType = 0;
                  which = "";
                  arpack_mode = 3;
                  // CBM: new stuff
                  lbound = 0.0;
                  ubound = 0.0;
                  nshifts = 0;
                  neigps = 50;
                  filtereig = false;
                  maxArnItr = 0;

                  noninpc = false;
                  inpc = false;
                  nsample = 1;

                  subspaceSize = -1;
                  tolEig = tolJac = 0.001;
                  nEig = 1;
                  maxitEig = 0;

                  sloshing = 0;
                  HEV = 0;
                  
                  mumps_icntl[3] = 0; // supress diagnostic output
                  //mumps_icntl[7] = 7; // renumbering, default = 7 (auto)
                  //mumps_icntl[8] = 7; // scaling, default = 7 (auto)
                  //mumps_cntl[1] = 0.01; // relative threshold for numerical pivoting (larger value may increase fill in but lead to more accurate factorization)

                  for(int i=0; i<20; ++i) { debug_icntl[i] = 0; debug_cntl[i] = 0.0; }

                  // iterative solver defaults
                  precond = 0; 
                  tol = 1.0e-8; 
                  maxit = 1000;
                  iterType = 0; 
                  iterSubtype = 3; 
                  maxvecsize = 0; 

                  iacc_switch = true;
                  zeroRot = false;

                  dist_acme = 0;
                  allproc_acme = true;
                  ffi_debug = false;
                  mortar_scaling = 1.0;
                  mortar_integration_rule = 6;
                  tdenforceFlag = true;

                  lagrangeMult = true;
                  penalty = 0;
                  num_penalty_its = 1;
                  penalty_tol = 1e-8;
                  penalty_beta = 10;
                  reinit_lm = false;
                  lm_update_flag = 1;
                  mpcDirect = 0;
                  usePrescribedThreshold = false;
                  mpcDirectTol = 10;
                  coefFilterTol = 10;
                  rhsZeroTol = 0;
                  inconsistentTol = 1e-8;
                  constraint_hess = 1;
                  constraint_hess_eps = 0;

                  numSnap            = 1;
                  readInModes        = "";
                  readInShapeSen     = "";
                  SVDoutput          = "pod.rob";
                  reducedMeshFile    = "";
                  isvPodRomFile      = "";
                  //dsvPodRomFile      = "";
                  forcePodRomFile    = "";
                  constraintPodRomFile  = "";
 		  constraintSnapshotFile = "";
                  constraintViolationFile = "";
                  residualPodRomFile = "";
                  jacobianPodRomFile = "";
                  ROMPostProcess     = false;
                  statevectPodRom    = false;
                  velocvectPodRom    = false;
                  accelvectPodRom    = false;
                  isvPodRom          = false;
                  dsvPodRom          = false;
                  forcevectPodRom    = false;
                  residvectPodRom    = false;
                  jacobvectPodRom    = false;
                  readmodeCalled     = false;
                  modalCalled        = false;
                  modalLMPC          = false;
                  readShapeSen       = false;
                  activatePodRom     = false;
                  snapshotsPodRom    = false;
                  checkPodRom        = false;
                  svdPodRom          = false;
                  DEIMBasisPod       = false;
                  UDEIMBasisPod      = false;
                  ConstraintBasisPod = false;
                  ReducedStiffness   = false;
                  computeForceSnap   = false;
                  computeConstraintSnap = false;
                  filterSnapshotRows = false;
                  orthogForceSnap    = false;
                  orthogConstraintSnap = false;
                  computeDEIMIndices = false;
                  DEIMPodRom         = false;
                  UDEIMPodRom        = false;
                  svdBlockSize       = 64;
		  clusterSubspaceAngle = false;
                  clustering         = 0;
                  solverTypeCluster  = 1; // K-means
                  use_nmf            = 0;
                  nmfNumROBDim       = 1;
                  nmfDelROBDim       = 10;
                  nmfRandInit        = 1;
                  nmfMaxIter         = 100;
                  nmfNumSub          = 0;
                  nmfTol             = 1e-6;
                  nmfPqnNumInnerIter = 2;
                  nmfPqnAlpha        = 0.4;
                  ksParameter        = 50;
                  ksMax              = 1.0e5;
                  samplingPodRom     = false;
                  snapProjPodRom     = false;
                  galerkinPodRom     = false;
                  elemLumpPodRom     = false;
                  onlineSvdPodRom    = false;
                  maxSizePodRom      = 0;
                  romEnergy          = 0.0;
                  maxSizeDualBasis   = 0;
                  maxDeimBasisSize   = 0;
                  selectFullNode     = false;
                  selectFullElem     = false;
                  forcePodSize       = 0;
                  constraintPodSize  = 0;
                  normalize          = 1;
                  subtractRefPodRom  = false;
                  useScalingSpnnls   = true;
                  useCenterSpnnls    = false;
                  useReverseOrder    = false;
                  projectSolution    = false;
                  positiveElements   = true;
                  maxSizeSpnnls      = 1.0;
                  maxElemSpnnls      = 0;
                  maxIterSpnnls      = 3.0;
                  solverTypeSpnnls   = 0;
                  reduceFollower     = false;
                  randomVecSampling  = false;
                  skipPodRom         = 1;
                  skipOffSet         = 0;
                  skipState          = 1;
                  skipVeloc          = 1;
                  skipAccel          = 1;
                  skipInternalStateVar = 1;
                  skipDualStateVar   = 1;
                  skipForce          = 1;
                  skipResidual       = 1;
                  skipJacobian       = 1;
                  orthogPodRom       = 1;
                  numRODFile         = 0;
                  romresidType       = 0;
                  tolPodRom          = 1.0e-6;
                  useMassNormalizedBasis = true;
                  useMassOrthogonalProjection = false;
                  performMassNormalization = false;
                  ConwepOnOff        = false;
                  basicDofCoords     = true;
                  basicPosCoords     = true;
                  scalePosCoords     = false;
                  xScaleFactor       = 1.0;
                  yScaleFactor       = 1.0;
                  zScaleFactor       = 1.0;
		  xLMPCFactor        = 1.0;
                  yLMPCFactor        = 1.0;
                  zLMPCFactor        = 1.0;
                  inertiaLumping     = 0;
                  printMatLab        = false;
                  printMatLabFile    = "";
                  elementDeletion    = false;
                  piecewise_contact  = true;
                  npMax              = 0;  // 0 => reduced mesh size is not limited
                  scpkMB             = 0;  // 0 => Scalapack LH solver will use default block size
                  scpkNB             = 0;  // 0 => Scalapack LH solver will use default block size
                  scpkMP             = 0;  // 0 => Scalapack LH solver will use default processor grid 
                  scpkNP             = 0;  // 0 => Scalapack LH solver will use default processor grid
                  useMassAugmentation = true;
                  piecewise = false;
                  freeplay = false;
                  piecewise_dlambda = 1.0;
                  piecewise_maxLambda = 1.0;
                }

   void setDirectMPC(int mode) { mpcDirect = mode; }
   // Whether we are doing direct elimination for MPCs
   int getDirectMPC() { return mpcDirect; }

   // Set RbmFilter level
   void useRbmFilter(int rbmfil) { filterFlags = rbmfil; }

   // Set Aeroelastic Algorithm
   void setAero(int alg)
    { if(alg == 22) { aeroFlag = 20; dyna3d_compat = true; }
      else aeroFlag = alg;
    }

   void setAeroHeat(int alg, double alt0 =0, double alt1 =0)
    { aeroheatFlag = alg;
      alphat[0] = alt0+alt1;
      alphat[1] = -alt1;
    }

   void setThermoe(int alg)
      { thermoeFlag = alg; }

   void setThermoh(int alg)
      { thermohFlag = alg; }

   void setModeDecomp(int alg)
      { modeDecompFlag = alg; }

   // SET CONDITION NUMBER TOLERANCE
   void setCondNumTol(double tolerance, int maxit) { condNumTolerance = tolerance; condNumMaxit = maxit; }

   // ... NON LINEAR SOLVER INFORMATION VARIABLES
   NonlinearInfo NLInfo;
   const NonlinearInfo &getNLInfo() const { return NLInfo; }
   NonlinearInfo &getNLInfo() { return NLInfo; }

   bool unsym() { return NLInfo.unsymmetric; }

   int gepsFlg;         // Geometric pre-stress flag
   int buckling;        // Buckling analysis flag
   void setGEPS() { gepsFlg  = 1; }

   // This could be a pointer to a FetiInfo type
   FetiInfo fetiInfo;

   FetiInfo &getFetiInfo() { return fetiInfo; }

   // KHP: MOVE TO NonlinearInfo
   void setNewton(int n)     { NLInfo.updateK    = n; }
   void setKrylov()          { NLInfo.kryflg     = 1; }
   void setInitialization()  { NLInfo.initflg    = 1; }
   void setReOrtho()         { NLInfo.reorthoflg = 1; }

   // SET DYNAMIC VALUE FUNCTIONS

   void setTimes(double _tmax, double _dt, double _dtemp)
   { tmax = _tmax; dt = _dt; dtemp = _dtemp; }

   void setParallelInTime(int J, int k, int workloadMax)
   { pitaTimeGridRatio = J; pitaMainIterMax = k; pitaProcessWorkloadMax = workloadMax; }

   // Set Rayleigh damping stiffness coefficient beta
   // Set Rayleigh damping mass coefficient alpha
   void setDamping(double beta, double alpha)
   { alphaDamp = alpha; betaDamp = beta; }
   bool hasDamping() { return ((alphaDamp != 0.0) || (betaDamp != 0.0)); }
   void setSDamping(double eta, double beta)
   { etaDamp = eta; betaDamp = beta; }
   bool hasSDamping() { return ((etaDamp != 0.0) ); }

   // SET RENUMBERING SCHEME (ND = nested dissection, MMD = multiple minimum degree, MS = multi-section, RCM = Reverse Cuthill McKee
   void setRenum(int renumberid) { 
     switch(renumberid) { 
       case 1 : renum = 1; break; // sloan
       case 2 : renum = 2; break; // RCM
     }
   }
   void setSparseRenum(int renumberid) { 
     switch(renumberid) {
       case 3 : sparse_renum = 0; break; // esmond MMD
       case 4 : sparse_renum = 1; break; // metis ND
     }
   }
   void setSpoolesRenum(int renumberid) {
     switch(renumberid) {
       case 6 : spooles_renum = 0; break; // spooles best of ND and MS
       case 3 : spooles_renum = 1; break; // spooles MMD
       case 5 : spooles_renum = 2; break; // spooles MS
       case 4 : spooles_renum = 3; break; // spooles ND
     }
   }

   void setProbType(int pbt)
   { 
     if((probType == Top) || (probType == Decomp)) return;

     // PJSA 6-9-04: ignore STATIC keyword if different probType already defined
     if((pbt == Static) && (probType != None)) return; 

     // PJSA 4-2-08
     if((pbt == Dynamic) && (probType != None)) { setDynamicProbType(); return; }

     probType = pbt;
   }

   void setDynamicProbType()
   {
     switch(probType) {
       case(Top) : 
       case(Decomp) :
         break;
       case(NonLinStatic) : 
       case(NonLinDynam) :
         probType = NonLinDynam;
         break;
       case(ArcLength) :
         std::cerr << "ARCLENGTH + DYNAMICS not supported\n"; 
         exit(-1);
         break;
       case(MatNonLinDynam) :
       case(MatNonLinStatic) :
         probType = MatNonLinDynam;
         break;
       default :
         if(probType != None && probType != Static && probType != Dynamic) std::cerr << "WARNING: switching problem type from " << probType << " to Dynamic\n";
         probType = Dynamic;
         break;
     }
   }

   // set parameters for 2nd order dynamics problem using newmark/generalize-alpha time integration scheme
   void setNewmarkSecondOrderInfo(double beta = 0.25, double gamma = 0.50, 
                                  double alphaf = 0.50, double alpham = 0.50, double rho=1.0)
   { 
     timeIntegration = Newmark;
     order = 2;
     if(rho != 1.0) { // compute alpha parameters from spectral radius if given
       alphaf = rho/(rho+1.0);
       alpham = (2.0*rho-1.0)/(rho+1.0);
       beta = 0.25*(1+alphaf-alpham)*(1+alphaf-alpham);
       gamma= 0.5+alphaf-alpham;
     }
     else if ((alphaf == 10.0) && (alpham == 10.0)) {
       alphaf = 0.5;
       alpham = 0.5;
       beta = 0.25;
       gamma= 0.5;
     }
     if(alphaf >= 1.0) alphaf = 0.99;
     newmarkBeta = beta; 
     newmarkGamma = gamma;
     newmarkAlphaF = alphaf; 
     newmarkAlphaM = alpham;
   }

   //set parameters for a 1st order dynamics problem using newmark time integration scheme
   void setNewmarkFirstOrderInfo(double epsiln = 0.5)
   {
     timeIntegration = Newmark;
     order = 1;
     newmarkGamma = alphaTemp = epsiln;
     if(newmarkGamma == 0) newmarkBeta = 0; // this is used to flag explicit scheme
   }

   // set parameters for quasistatic problem
   void setQuasistaticInfo(double _qsMaxvel, double _qsBeta,
                           double _steadyTol, int _steadyMax, double _delta = 0) 
   {
     timeIntegration = Qstatic;
     qsMaxvel = _qsMaxvel;
     qsBeta = _qsBeta;
     steadyTol = _steadyTol;
     steadyMax = _steadyMax;
     delta = _delta;
     if(hzemFilterFlag == 1) {
       fprintf(stderr," ... WARNING: HZEMFILTER IGNORED FOR QSTATICS ...\n");
       hzemFilterFlag = 0;
     } 
     steadyFlag = 1;
     if (tmax == 0.0) { tmax = 1.0e8; dt = 1.0; dtemp = 1.0; }

     // delta is a translation from iterations to time 
     if(delta)  {
       qsMaxvel = 1.0;
       fprintf(stderr, " ... Setting Max Quasi-Static Velocity Parameter to 1.0\n");
     }
   }

   double getTimeStep() const
   {
     if(thermoeFlag > -1 || thermohFlag > -1) return std::min(dt, dtemp);
     else return (order == 1) ? dtemp : dt;
   }

   void setTimeStep(double _dt)
   {
     if(order == 1) dtemp = _dt;
     else dt = _dt;
   }
 
   // set parameters for eigen problem  
   void setSubSpaceInfo(int _subspaceSize, double _tolEig, double _tolJac)
   { 
     eigenSolverType = SubSpace; 
     subspaceSize = _subspaceSize; 
     tolEig = _tolEig; 
     tolJac = _tolJac; 
   }

   // direct solver
   void setSolver(int _substype) { type = 0; subtype = _substype;  }

   // iterative solver
   void setSolver(int _iterType, int _precond, double _tol=1.0e-8, int _maxit=1000,
                  int _iterSubtype=3, int _maxvecsize=0)
    { type = 1; precond = _precond; tol = _tol; maxit = _maxit; 
      iterType = _iterType; iterSubtype = _iterSubtype; maxvecsize = _maxvecsize; 
    }

   void setTrbm(double _tolzpv)
    { trbm = _tolzpv; mumps_cntl[3] = trbm; }

   void setGrbm(double _tolsvd, double _tolzpv)
    { trbm2 = _tolzpv; tolsvd = _tolsvd; rbmflg = 1; }
   void setGrbm(double _tolzpv)
    { trbm2 = _tolzpv; rbmflg = 1; }
   void setGrbm() 
    { rbmflg = 1; }

   bool isAcousticHelm() {
     return ((probType == Helmholtz) || (probType == HelmholtzDirSweep)
             || (probType == HelmholtzMF) || (probType == HelmholtzSO)
             || (probType == HelmholtzFreqSweep));
   }
   bool isAcoustic() {
     return (isAcousticHelm()||acoustic);
   }
   
   bool isATDARB() {
     return (ATDARBFlag != -2.0);
   }

   bool isDynam() {
     return ((timeIntegration != Qstatic) && ((probType == Dynamic) || (probType == NonLinDynam)
             || (probType == TempDynamic) || (probType == MatNonLinDynam) || (probType == PodRomOffline)));
   }

   bool isNonLin() {
     return ((probType == NonLinStatic) || (probType == NonLinDynam)
             || (probType == MatNonLinStatic) || (probType == MatNonLinDynam) || (probType == ArcLength)
             || (probType == PodRomOffline));
   }

   bool isNonLinExtF() {
     return (isNonLin() && NLInfo.linearelastic != 2);
   }

   bool isStatic() {
     return ((probType == Static) || (probType == NonLinStatic)
             || (probType == MatNonLinStatic) || (probType == ArcLength));
   }

   bool keepModalInitialConditions() {
     if(galerkinPodRom) {
       // for nonlinear ROMs, keep the modal ivel and non-modal ivel separate
       // only if the basis used to define the initial conditions is the same
       // as the ROB.
       std::string rob(readInROBorModes[0]);
       if(useMassNormalizedBasis) rob.append(".normalized");
       std::string icb(readInModes);
       return (rob.compare(icb) == 0); 
     }
     else {
       // for linear ROMs, keep the modal ivel and non-modal ivel separate
       return modal;
     }
   }

   int classifySolver();
   void activatePiecewise();
};

#endif
