#ifndef _STRUCTOPT_H_
#define _STRUCTOPT_H_

#ifdef STRUCTOPT

#include <Utils.d/MyComplex.h>

#include <Structopt.d/Optcrit.h>
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optsol.h>

#include <Driver.d/Domain.h>
#include <Structopt.d/Driver_opt.d/Domain_opt.h>

#include <Problems.d/StaticDescr.h>
#include <Structopt.d/Problems_opt.d/StaticDescr_opt.h>

#include <Driver.d/StaticProbType.h>
#include <Structopt.d/Driver_opt.d/StaticProbType_opt.h>

#include <Math.d/Vector.h>
#include <Math.d/VectorSet.h>
#include <Math.d/SparseMatrix.h>

#include <Utils.d/resize_array.h>

#include <Solvers.d/Solver.h>
#include <Element.d/Element.h>

#include <Structopt.d/Structopt_base.h>

#define numBaseOptpar 5


//------------------------------------------------------------------------------
class Structopt_sd : public Structopt
{
 public:
  int numvar;
  int type; //0 - Interface for Design Optimization
            //1 - Interface for Reliability Analysis
  int numAnalysis;      
  int numFunc;
  int anagrdType;
  int optInfset;
  Domain_opt *structdom;
  Optpro *optpro;
  
  int reliabilityFlag;
  Optpro*       reliabilityProb; 
  Structopt_sd* reliabilityStrc;
  
  int failcritOnly;
  
  int designoptFlag;
  Structopt_sd* designoptStrc;
  
  StcAnalysisData *analysisData;
  
  int aeroact;                       //Flag to save Sens. of Aeroforces
  
  int sendInitDisp;
  
  int numElecAttrVar;
  int numGradElecAttrVar;
  int numCurElecAttrVar;
  
  int numThermAttrVar;
  int numGradThermAttrVar;
  int numCurThermAttrVar;
  
  OptActInfo** optInf;
  
  //ResizeArray<Stcvar*>   opv;        //Structural Variables     

  //for single domain static analysis
  int numStaticProb; 
  SingleDomainStatic_opt<double, Vector, Solver> **staticpros;  
  StaticSolver_opt< double, Solver, Vector,
    SingleDomainPostProcessor_opt<double, Vector, Solver>,
    SingleDomainStatic_opt<double, Vector, Solver>,
    ComplexVector > **staticsolvers;

  SingleDomainStatic_opt<DComplex, ComplexVector, ComplexSolver> **c_staticpros;
  StaticSolver_opt<DComplex, ComplexSolver, ComplexVector,
    SingleDomainPostProcessor_opt<DComplex, ComplexVector, ComplexSolver>,
    SingleDomainStatic_opt<DComplex, ComplexVector, ComplexSolver>,
    ComplexVector > **c_staticsolvers;
  
  /*
  //for nonlinear static analysis  
  int numNlnstcProb;  
  NonLinStatic **nonlinpros;  
  NLStaticSolver <Solver,Vector,SingleDomainPostProcessor<double, Vector, Solver>,
    NonLinStatic, GeomState> **nonlinsolvers;
  
  //for single domain eigenvalue analysis  
  int numEigenProb;  
  SingleDomainEigen **eigenpros;  
  EigenSolver<DynamMat, Vector, VectorSet, 
    SDEigenPostProcessor, SingleDomainEigen> 
    **eigensolvers;
  
  //for single linear domain dynamic and multiphysics analysis  
  int numDynamProb;  
  SingleDomainDynamic **dynpros;
  DynamicSolver <DynamMat,Vector,SDDynamPostProcessor,SingleDomainDynamic>
    **dynsolvers;
  
  
  //for single nonlinear domain dynamic and multiphysics analysis  
  int  numNLdynProb;  
  NonLinDynamic **NLdynpros;  
  NLDynamSolver  <Solver, Vector, SDDynamPostProcessor, 
    NonLinDynamic, GeomState>
    **NLdynsolvers;
  
  // forced vibrations
  int numFVibrProb;
  SingleDomainFVibr** fvibrpros;
  FVibrSolver**       fvibrsolvers;
  */
  
  // pointers to various solutions and fem arrays  
  int * analysisMap;
  int   activeAnalysis;
  
  double * bcx;
  
  Vector * sol;
  Vector * vel;
  Vector * acc;
  
  Vector * grad;
  Vector * adj;
  
  GeomState*  gs;
  Corotator** aC;
  double lam;
  double lam_anal;
  double dqdlam;
  
  SparseMatrix * gStiff;
  SparseMatrix * gMass;
  SparseMatrix * gDamp;
  
  FullSquareMatrix * kelArray;
  
  ComplexVector* csol;
  ComplexVector* cgrad;
  ComplexVector* cadj;
  DComplex*      cbcx;
  
 public:  
  // Constructors  
  explicit Structopt_sd(int, Domain_opt* d=0, Optpro* o=0);
  
  void build      (Domain_opt*, Optpro*);
  void buildInopt (Domain_opt*, Optpro*, Structopt_sd*);
  
  void initStructProbSol(Structopt_sd* stcopt=0);
  void initReliabilityAnalysis();
  void initAnalyticalSA();
  void cleanup();
  
  // set pointer to attributes, coordinates, nodal forces
  
  double * getptrnodattr    ( int& loc1, int& loc2 );
  double * getptrnodcord    ( int& loc1, int& loc2 );
  double * getptrnodalforce ( int& loc1, int& loc2 );
  double * getptrfluidattr  ( int& loc1);
  double * getptrelecattr   ( int& num,  int& loc1, int& loc2 );
  double * getptrthermattr  ( int &num, int& loc1, int& loc2 );
  
  // set pointer to attributes and coordinates
  
  double * getptrgradnodattr    ( int& loc1, int& loc2 );
  double * getptrgradnodcord    ( int& loc1, int& loc2 );
  double * getptrgradnodalforce ( int& loc1, int& loc2 );
  double * getptrgradfluidattr  ( int& loc1);
  double * getptrgradelecattr   ( int& num,  int& loc2 );
  double * getptrgradthermattr  ( int &num, int& loc2 );
  
  // arrays to store electrostatic variables
  
  double           *elecattr; 
  double           *gradelecattr; 
  ResizeArray<int>  elecattrList;
  
  // arrays to store thermal variables
  
  double	       *thermattr; 
  double	       *gradthermattr; 
  ResizeArray<int>  thermattrList;
  
  // get design criteria values
  
  void getStressInfo(int,int&,int&);
  
  void func();
  void evaluate();
  
  int getAnalysisData(int, double&);
  int getAnalysisGradData(int, double&);
  int getAnalysisAdjData(int, double&, int=0);
  
  double getstrainenergy(int*,int,double&,int);
  double getmass(int*,int,int);
  double getMomOfInertia(int*,int,int,int);
  double getfrequency(int,int);
  double getcontrolcost(int);
  double getdisp(int,int,int,double&,int);
  double getLambda(double&,int);
  double getinternalforce(int,int,double&,int);
  double getnodalpos(int,int,double&,int);
  double getnodstr(int,int,double&,int);
  double getkineticenergy(double&,int);
  double getdampingenergy(double&,int);
  double getExtVal(int,int,int,double&);
  double getaeroforce(int,double&,int);
  double getstressint(int,double&,double&,int,int,int*,int,double&,int);
  double getdisplevel(int,int,int,int*,int*,double*,double*,double&,double&,
		      int);
  double getstresslevel(int,int,int,int,int*,int*,double*,double*,double&,
			double&,int);
  double getforml2norm2();
  double getFailprob(int,int);
  double getestat(int,int*,int,double&,int);
  double getPullIn(int);
  double getDCmass(double&,int*,int,double&,int);
  
  // sensitivity analysis
  
  void initOptInf (int);
  void setoptInf(int);
  void addOptInf(int,int,int,int);
  
  void getMidPointMass(double*,double**);
  
  void zerograd();
  void graddirect(int,double**,int n=0,int* nn=0); 
  void gradadjoint(int, double**);   
  void gradadjointStress(int, double**);   
  void gradadjointModal (int, double**);   
  void gradadjointGeom  (int, double**);   
  
  double getgradstrainenergy(int*,int,double&,int);
  double getgradmass(int*,int,int);
  double getgradMomOfInertia(int*,int,int,int);
  double getgradfrequency(int,int);
  double getgraddisp(int,int,int,double&,int);
  double getgradLambda(double&,int);
  double getgradinternalforce(int,int,double&,int);
  double getgradnodalpos(int,int,double&,int);
  double getgradnodstr(int,int,double&,int);
  double getgradkineticenergy(double&,int);
  double getgraddampingenergy(double&,int);
  double getgradaeroforce(int idir,double&,int);
  double getgradstressint(int,double&,double&,int,int,int*,int,double&,int);
  double getgraddisplevel(int,int,int,int*,int*,double*,double*,double&,
			  double&,int);
  double getgradstresslevel(int,int,int,int,int*,int*,double*,double*,
			    double&,double&,int);
  double getgradforml2norm2();
  double getgradFailprob(int,int,int=-1);
  double getgradestat(int,int*,int,double&,int);
  double getgradPullIn(int,int);
  double getgradDCmass(double&,int*,int,double&,int);
  
  double getgradpartstrainenergy(int*,int,double&,int);
  double getgradpartinternalforce(int,int,double&,int);
  double getgradpartnodstr(int,int,double&,int); 
  double getgradpartstressint(int,double&,double&,int,int,int*,int,double&,
			      int);    
  double getgradpartstresslevel(int,int,int,int,int*,int*,double*,double*,
				double&,double&,int);
  double getgradpartestat(int,int*,int,double&,int);
  double getgradpartPullIn(int);
  double getgradpartDCmass(double&,int*,int,double&,int);
  
  void getgraddustrainenergy(int*,int,double&,int);
  void getgraddudisp(int,int,int,double&,int);
  void getgraddulambda(double&,int);
  void getgradduinternalforce(int,int,double&,int);
  void getgraddunodstr(int,int,double&,int);
  void getgraddumass(int*,int,int);
  void getgradduMomOfInertia(int*,int,int,int);
  void getgradduaeroforce(int,double&,int);
  void getgradduExtVal(int,int,double&);
  void getgraddustressint(int,double&,double&,int,int,int*,int,double&,int);
  void getgraddudisplevel(int,int,int,int*,int*,double*,double*,double&,
			  double&,int);
  void getgraddustresslevel(int,int,int,int,int*,int*,double*,double*,
			    double&,double&,int);
  void getgradduDCmass(double&,int*,int,double&,int);
  
  // transient related functions (mostly for aeroelastic)
  
  int     maxstep;
  double  stepcoef;
  
  void   addAnalysisData   ( anadata &);


#ifdef AEROELASTIC	  
  void   sndOptpar         ( int, int);
  void   setFluidvariables ( double*, int);
  void   setFluidcriteria  ( int, double*);
  
  void   setElecvariables  (double*,int,int);
  void   setThermvariables (double*,int,int);
  
  int    chkOptInfEL();
  int    chkOptInfTH();
#endif
  
  int    procComputeMaxSteps( double&);
  int    procSetvar         ( double, double);
  int    procTransition     ( int,    int);
  void   procEvaluate       ( double, double);
  void   procComputeStepSize( int); 
  void   procSavevar();
  
  int    initProcess(double&, double&);
  
  // output routines  
  void postProcessing(int);
};

#endif

#endif
