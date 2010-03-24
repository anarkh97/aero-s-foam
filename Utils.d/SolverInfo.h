#ifndef _SOLVER_INFO_
#define _SOLVER_INFO_

#include <Utils.d/NonlinearInfo.h>
#include <Feti.d/FetiInfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <map>

#if defined(WINDOWS) || defined(MACOSX)
 #include <cfloat>
#else
 #include <limits>
#endif

using namespace std;

struct SolverInfo {

   // Problem Type parameters
   enum { Static, Dynamic, Modal, NonLinStatic, NonLinDynam, 
	  ArcLength, ConditionNumber, TempDynamic, Top,
          AxiHelm, MatNonLinStatic, MatNonLinDynam,
          Helmholtz, HelmholtzFreqSweep, HelmholtzDirSweep, HelmholtzMF, HelmholtzSO,
          Decomp, None, NonLinTempDynam, DisEnrM };
   
   int probType;
   int soltyp; // from CONTROL statement: 1 = statics, 2 = heat conduction, etc...

   float ATDARBFlag;
   float ATDDNBVal;
   float ATDROBVal;
   float ATDROBalpha;
   float ATDROBbeta;
   int aeroFlag;
   int aeroheatFlag;
   int thermoeFlag;
   int thermohFlag;
   int modeDecompFlag;
   int isCollocated;
   double alphas[2];
   double alphat[2];

   int structoptFlag;
   int thermalLoadFlag; // whether there is a thermal load applied
   int hzemFlag; // Zero energy mode (for thermal Problems),
                 // Equivalent to grbm for structures
   int hzemFilterFlag;
   int slzemFlag;// Flag for zero energy mode, ADDED FOR SLOSHING PROBLEMS, EC, 20070723
   int slzemFilterFlag; // Flag for zero energy mode filter, ADDED FOR SLOSHING PROBLEMS, EC, 20070723

   double mppFactor;  // modal amplification factor for mpp command

   // Solver parameters
   int type;     // 0 = direct, 1 = iterative, 2 = FETI, 3 = Block Diag
   int subtype;  // subtype and matrix storage... 9 is mumps  10 is diag
   int iterType; // 0 = CG, 1 = CR, 2 = BCG
   int precond;  // preconditioner 0 = none, 1 = jacobi
   int maxit;    // maximum number of iterations
   double tol;   // tolerance for convergence
   int renum;    // renumbering scheme for skyline: 0 = none (default), 1 = sloan, 2 = RCM
   int maxvecsize;  // for pcg # of krylov vectors to store default = 0
   bool lastIt;

   // Time Integration Algorithms
   enum { Newmark, Qstatic };

   // Parameters for PITA
   bool tiParall; 	    // Enables PITA (both linear and nonlinear problems)
   bool mdPita;             // Enables time-space parallelism (linear problem only)
   bool NoForcePita;        // Enables NoForce optimization (linear problem only)
   bool ConstForcePita;     // Enables ConstantForce optimization (linear problem only)
   bool CkCoarse;           // Enable coarse time-grid correction (linear problem only)
   int Jratio;              // Coarse/fine time-grid ratio (always required)  
   int kiter;               // Maximum number of main iterations (always required)
   int numTSperCycleperCPU; // Maximum number of active time-slices on each CPU (always required)
   int numSpaceMPIProc;     // Number of CPUs in the space domain (only when time-space parallelism is enabled)
   int baseImprovementMethodForPita; // 0 = All seeds, 1 (default) = Local increments (nonlinear problem only) 
   bool remoteCoarse;       // Coarse grid integrator on dedicated CPU 
   double pitaProjTol;      // Tolerance used to build the projector

   bool newPitaImplementation; // JC TODO Remove when new PITA implementation ready

   bool modal;          // true iff system is to be solved in modal coordinates
   bool acoustic;       // true iff system is to be solved for acoustic time domain
   bool modifiedWaveEquation; // true if solving using the modified wave equation
   double modifiedWaveEquationCoef; // value for the coefficient (theoretically 12.0)
   bool modalIDisp;     // true iff initial displacements are given in modal coord.
   int order;           // used for dynamics... 1: first order (heat), 2: second order (mech, acoustics)
   int timeIntegration; // type of time integration scheme
   int initialTimeIndex;// initial time index (either 0 or from restart)
   double initialTime;  // initial time (either 0.0 or from restart)
   double initExtForceNorm;  // initial Force Norm (for restarting qstatics)
   double tmax;         // maximum time
   double dt;           // structural time step value
   double dtemp;        // thermal time step value
   double alphaDamp;    // Rayleigh Mass damping coefficient 
   double betaDamp;     // Rayleigh Stiffness damping coefficient
   double alphaTemp;
   double newmarkBeta;  // Newmark algorithm parameter (beta)
   double newmarkGamma; // Newmark algorithm parameter (gamma)
   double newmarkAlphaF;  // Newmark algorithm parameter (alphaf)
   double newmarkAlphaM; // Newmark algorithm parameter (alpham)
   bool stable;         // compute stability timestep for explicit
   double stable_tol;   // convergence tolerance for computionation of stability timestep
   int stable_maxit;    // stable_maxit*n is the maximum number of iterations for computation of stability timestep
   double qsMaxvel;     // final relaxation parameter in quasi-static alg.
   double delta;        // translation factor from time index to real time for quasistatics
   int nRestart;        // how many time steps per restart
   int steadyFlag;      // quasi-transient steady state simulation
   int steadyMin;       // minimum number of steps in steady state simulation
   int steadyMax;       // maximum number of steps in steady state simulation
   double steadyTol;    // steady state criteria
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
    
   // Sloshing problem flag, ADDED FOR SLOSHING PROBLEM, EC, 20070723
   int sloshing;
   // Hydroelastic vibration problem flag, ADDED FOR HEV PROBLEM, EC, 20070820
   int HEV;
 
   const char *which;   /* for ARPACK: specifies "which" of the eigenvalues to compute, where:
                          "LA" - compute the eigenvalues just to the right of sigma (sigma is the shift)
                          "SA" - compute the eigenvalues just to the left of sigma
                          "LM" - undefined for shift-invert on generalized eigenvalue problem
                          "SM" - undefined for shift-invert on generalzied eigenvalue problem
                          "BE" - compute eigenvalues to either side of sigma */
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
   double tolsvd;       // singular value decomposition tolerance
   int rbmflg;          // 0 = algebraic rbm, 1 = geometric rbm
   int rbmFilters[6];   // rbm filtering for nonlinear modal problems for each
                        //   of the six rigid body modes: 0 = do not filter, 1 = do filter

   double condNumTolerance; // Condition number tolerance
   int condNumMaxit;

   int massFlag;
   int filterFlags;

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

   // KAS :  map object for Mumps control CNTL and ICNTL matrices
   map<int, int> mumps_icntl;
   map<int, double> mumps_cntl;

   bool localScaled, coarseScaled;

   int nFreqSweepRHS;
   bool doFreqSweep,doEigSweep;
   //--- UH --- 05/21/08
   enum { Taylor, Pade1, Pade, Fourier, PadeLanczos, GalProjection };
   //--- UH --- 05/21/08
   int freqSweepMethod;
   int padeL, padeM, padeN;
   bool pade_pivot;
   double pade_tol;
   bool pade_poles;
   double pade_poles_sigmaL, pade_poles_sigmaU;
   int modeFilterFlag;
   bool test_ulrich;
   int addedMass;
   bool isCoupled;
   double coupled_scale;
   bool isMatching;  // true if only one wet interface is given, otherwise false (default) 
   bool farfield; // true if farfield output requested

   bool dbccheck;
   int contact_mode;
   double coulomb_fric_coef;

   bool noninpc;
   bool inpc;
   int nsample;

   map<int, int> debug_icntl;   // used for debugging
   map<int, float> debug_cntl;  // used for debugging

   bool iacc_switch; // mech/acou: true --> compute consistent initial second time derivative ie, a^0 = M^{-1}(fext^0 - fint^0 - Cv^0) for a second order differential equation (ie mech/acou)
                     //            false --> a^0 = 0
                     // heat: true --> compute consistent initial first time derivative ie, v^0 = M^{-1}(fext^0 - fint^0) for a first order differential equation (ie heat)
                     //       false --> v^0 = 0

   // Constructor
   SolverInfo() { filterFlags = 0;
                  NLInfo = 0; 
                  type = 0;     
                  soltyp = -1;
                  subtype = 0; // By default we use direct Skyline
                  renum = 0;
                  probType = SolverInfo::None;
                  alphaDamp = 0.0;
                  betaDamp = 0.0;
                  modal = false;
                  modalIDisp = false;
                  lastIt = false;
                  mppFactor = 1.0;
		 
                  // Parameters for PITA 
		              tiParall       = false;
                  mdPita         = false;
                  NoForcePita    = false;
                  ConstForcePita = false;
                  CkCoarse       = false;
		  Jratio = 1;
		  kiter  = 0;
                  numTSperCycleperCPU = 1;
                  numSpaceMPIProc     = 1;
                  baseImprovementMethodForPita = 0;
                  remoteCoarse = false;
                  pitaProjTol = 1.0e-6;
                  newPitaImplementation = false; // JC TODO Remove when new PITA implementation ready

                  acoustic = false;
                  modifiedWaveEquation = false;
                  modifiedWaveEquationCoef = 12.0;

                  tmax = 0.0; 
                  dt = 0.0;
                  dtemp = 0.0;

                  steadyFlag = 0;
                  steadyMin = 1;
                  steadyMax = 10;
                  steadyTol = 1.0e-3; 
                  no_secondary = false;

                  trbm = 1.0E-16;   // default zero pivot tolerance
                  tolsvd = 1.0E-6;  // default singular value tolerance
                  massFlag = 0;     // whether to calculate total structure mass
				  
	          ATDARBFlag = -2.0;
                  ATDDNBVal = 0.0;
	   	  ATDROBVal = 0.0;
		  ATDROBalpha = 0.0; //this value can not be 0 when Robin boundary is set, it is the flag!
		  ATDROBbeta = 0.0;

		  aeroFlag = -1;
                  aeroheatFlag = -1;
                  thermoeFlag = -1;
                  thermohFlag = -1;
                  thermalLoadFlag = 0;
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
                  stable = true;
                  stable_tol = 1.0e-3;
                  stable_maxit = 100;
                  rbmflg = 0;
                  rbmFilters[0] = rbmFilters[1] = rbmFilters[2] = rbmFilters[3] = rbmFilters[4] = rbmFilters[5] = 0;
                  buckling = 0;

                  sparse_renum = 0;
                  sparse_maxsup = 100;
                  sparse_defblk = 30;
                  pivot = false;
                  spooles_scale = 0;
                  spooles_tau = 100.;
                  spooles_seed = 532196;
                  spooles_maxsize = 64;
                  spooles_maxdomainsize = 24;
                  spooles_maxzeros = 0.04;
                  spooles_msglvl = 0;
                  spooles_renum = 0;
                  explicitK = false;
                  localScaled = false;
                  coarseScaled = false;

                  nFreqSweepRHS = 8;
                  doFreqSweep = false;
                  doEigSweep = false;
                  freqSweepMethod = Taylor;
                  padeL = 9;
                  padeM = 10;
                  padeN = 2;
                  pade_pivot = false;
                  pade_tol = 1.0e-16;
                  pade_poles = false;
                  pade_poles_sigmaL = 0.0; pade_poles_sigmaU = numeric_limits<double>::max();
                  modeFilterFlag = 0;
                  test_ulrich = false;
                  addedMass = 1;
                  isCoupled = false;
                  coupled_scale = 1.0;
                  isMatching = false; 
                  farfield = false;

                  dbccheck = false;
                  contact_mode = 1;
                  coulomb_fric_coef = 0.0;

                  nEig = 0;
                  eigenSolverSubType = 0;
                  which = "";
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

                  //ADDED FOR SLOSHING PROBLEM, EC, 20070723
                  sloshing = 0;
                  HEV = 0;
                  
                  mumps_icntl[3] = 0; // supress diagnostic output
                  //mumps_icntl[7] = 7; // renumbering, default = 7 (auto)
                  //mumps_icntl[8] = 7; // scaling, default = 7 (auto)
                  //mumps_cntl[1] = 0.01; // relative threshold for numerical pivoting (larger value may increase fill in but lead to more accurate factorization)

                  for(int i=0; i<20; ++i) { debug_icntl[i] = 0; debug_cntl[i] = 0.0; }
                  // some temporary defaults for FETI contact
                  debug_icntl[0] = 0; debug_icntl[1] = 0; debug_icntl[3] = 0;

                  // iterative solver defaults
                  precond = 0; 
                  tol = 1.0e-6; 
                  maxit = 100;
                  iterType = 1; 
                  subtype = 3; 
                  maxvecsize = 0; 

                  iacc_switch = true;
                 }

   // Set RbmFilter level
   void useRbmFilter(int rbmfil) { filterFlags = rbmfil; }

   // Set Aeroelastic Algorithm
   void setAero(int alg) { aeroFlag = alg; }

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
   NonlinearInfo *NLInfo;
   void initNLInfo() { if(NLInfo == 0) NLInfo = new NonlinearInfo; }
   NonlinearInfo & getNLInfo() { return *NLInfo; }

   int gepsFlg;         // Geometric pre-stress flag
   int buckling;        // Buckling analysis flag
   void setGEPS() { if(NLInfo == 0) NLInfo = new NonlinearInfo; gepsFlg  = 1; }

   // This could be a pointer to a FetiInfo type
   FetiInfo fetiInfo;

   FetiInfo &getFetiInfo() { return fetiInfo; }

   // KHP: MOVE TO NonlinearInfo
   void setNewton(int n)     { if(NLInfo == 0) NLInfo = new NonlinearInfo; NLInfo->updateK    = n; }
   void setKrylov()          { if(NLInfo == 0) NLInfo = new NonlinearInfo; NLInfo->kryflg     = 1; }
   void setInitialization()  { if(NLInfo == 0) NLInfo = new NonlinearInfo; NLInfo->initflg    = 1; }
   void setReOrtho()         { if(NLInfo == 0) NLInfo = new NonlinearInfo; NLInfo->reorthoflg = 1; }

   // SET DYNAMIC VALUE FUNCTIONS

   void setTimes(double _tmax, double _dt, double _dtemp)
    { tmax = _tmax; dt = _dt; dtemp = _dtemp; }

   void setParallelInTime( int _Jratio, int _kiter, int _numTSperCycleperCPU)    
    { Jratio = _Jratio; kiter=_kiter; numTSperCycleperCPU=_numTSperCycleperCPU; }

   // Set Rayleigh damping stiffness coefficient alpha
   // Set Rayleigh damping mass coefficient beta
   void setDamping(double beta, double alpha)
   { alphaDamp = alpha; betaDamp = beta; }
   bool hasDamping() { return ((alphaDamp > 0.0) || (betaDamp > 0.0)); }

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
         cerr << "ARCLENGTH + DYNAMICS not supported\n"; 
         exit(-1);
         break;
       case(MatNonLinDynam) :
       case(MatNonLinStatic) :
         probType = MatNonLinDynam;
         break;
       default :
         if(probType != None && probType != Static && probType != Dynamic) cerr << "WARNING: switching problem type from " << probType << " to Dynamic\n";
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
     alphaTemp = epsiln;
     //if(alphaTemp < 0.5){
     //   fprintf(stderr," ... WARNING: UNSTABLE EXPLICIT SCHEME ...\n");
     // }
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
 
   // set parameters for eigen problem  
   void setSubSpaceInfo(int _subspaceSize, double _tolEig, double _tolJac)
   { 
     eigenSolverType = SubSpace; 
     subspaceSize = _subspaceSize; 
     tolEig = _tolEig; 
     tolJac = _tolJac; 
   }

   // DEFINE SOLVER
   void setSolver(int _substype) { type = 0; subtype = _substype;  }

   void setSolver(int _precond, double _tol,int _maxit,int _iterType=1,
                  int _subtype=3, int _maxvecsize=0)
    { type = 1; precond = _precond; tol = _tol; maxit = _maxit; 
      iterType = _iterType; subtype = _subtype; maxvecsize = _maxvecsize; }

   void setTrbm(double _tolzpv)
    { trbm = _tolzpv; rbmflg = 0; }

   void setGrbm(double _tolsvd, double _tolzpv)
    { trbm = _tolzpv; tolsvd = _tolsvd; rbmflg = 1; }
   void setGrbm(double _tolzpv)
    { trbm = _tolzpv; rbmflg = 1; }

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
     return ((probType == Dynamic) || (probType == NonLinDynam)
             || (probType == TempDynamic) || (probType == MatNonLinDynam));
   }

   bool isNonLin() {
     return ((probType == NonLinStatic) || (probType == NonLinDynam)
             || (probType == MatNonLinStatic) || (probType == MatNonLinDynam) || (probType == ArcLength));
   }

   bool isStatic() {
     return ((probType == Static) || (probType == NonLinStatic)
             || (probType == MatNonLinStatic) || (probType == ArcLength));
   }

};

#endif
