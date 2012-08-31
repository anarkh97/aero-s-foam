#ifndef __STRUCTOPT_BASE_HPP__
#define __STRUCTOPT_BASE_HPP__

#ifdef STRUCTOPT

#include <Structopt.d/Optinp.h>
#include <Structopt.d/Driver_opt.d/Domain_opt.h>
#include <memory>

//------------------------------------------------------------------------------
struct StcAnalysisData 
{
  int dynFlag;
  int dynTyp;
  
  int dampFlag;
  int dampTyp;
  int dampFrq[2];
  
  int    transFlag;
  int    transMin;
  int    transMax;
  int    initBlend;
  double designVel;
  double forceVel;
};

//------------------------------------------------------------------------------
class Structopt 
{
 public:  
  virtual ~Structopt();

  //virtual void build      (Domain_opt*, Optpro*) { assert(0); }
  //virtual void buildInopt (Domain_opt*, Optpro*, Structopt*);
  
  //virtual void initStructProbSol(Structopt* stcopt=0);
  virtual void initReliabilityAnalysis() { assert(0); }
  virtual void initAnalyticalSA() { assert(0); }
  virtual void cleanup() { assert(0); }
  
  // set pointer to attributes, coordinates, nodal forces  
  virtual double * getptreleattr    ( int& loc1, int& loc2 );
  virtual double * getptrnodattr    ( int& loc1, int& loc2 ) { assert(0); }
  virtual double * getptrnodcord    ( int& loc1, int& loc2 ) { assert(0); }
  virtual double * getptrcomposite  ( int& loc1, int& loc2 );
  virtual double * getptrnodalforce ( int& loc1, int& loc2 ) { assert(0); }
  virtual double * getptrfluidattr  ( int& loc1) { assert(0); }
  virtual double * getptrelecattr   ( int& num,  int& loc1, int& loc2 ) { assert(0); }
  virtual double * getptrthermattr  ( int &num, int& loc1, int& loc2 ) { assert(0); }
  
  // set pointer to attributes and coordinates  
  virtual double * getptrgradeleattr    ( int& loc1, int& loc2 );
  virtual double * getptrgradnodattr    ( int& loc1, int& loc2 ) { assert(0); }
  virtual double * getptrgradnodcord    ( int& loc1, int& loc2 ) { assert(0); }
  virtual double * getptrgradcomposite  ( int& loc1, int& loc2 );
  virtual double * getptrgradnodalforce ( int& loc1, int& loc2 ) { assert(0); }
  virtual double * getptrgradfluidattr  ( int& loc1) { assert(0); }
  virtual double * getptrgradelecattr   ( int& num,  int& loc2 ) { assert(0); }
  virtual double * getptrgradthermattr  ( int &num, int& loc2 ) { assert(0); }
  
  
  // get design criteria values  
  virtual void getStressInfo(int,int&,int&) { assert(0); }
  
  virtual void func() { assert(0); }
  virtual void evaluate() { assert(0); }
  
  virtual int getAnalysisData(int, double&) { assert(0); }
  virtual int getAnalysisGradData(int, double&) { assert(0); }
  virtual int getAnalysisAdjData(int, double&, int=0) { assert(0); }
  
  virtual double getstrainenergy(int*,int,double&,int) { assert(0); }
  virtual double getmass(int*,int,int) { assert(0); }
  virtual double getMomOfInertia(int*,int,int,int) { assert(0); }
  virtual double getfrequency(int,int) { assert(0); }
  virtual double getcontrolcost(int) { assert(0); }
  virtual double getdisp(int,int,int,double&,int) { assert(0); }
  virtual double getLambda(double&,int) { assert(0); }
  virtual double getinternalforce(int,int,double&,int) { assert(0); }
  virtual double getnodalpos(int,int,double&,int) { assert(0); }
  virtual double getnodstr(int,int,double&,int) { assert(0); }
  virtual double getkineticenergy(double&,int) { assert(0); }
  virtual double getdampingenergy(double&,int) { assert(0); }
  virtual double getExtVal(int,int,int,double&) { assert(0); }
  virtual double getaeroforce(int,double&,int) { assert(0); }
  virtual double getstressint(int,double&,double&,int,int,int*,int,double&,int) { assert(0); }
  virtual double getdisplevel(int,int,int,int*,int*,double*,double*,double&,double&,
		      int) { assert(0); }
  virtual  double getstresslevel(int,int,int,int,int*,int*,double*,double*,double&,
			double&,int) { assert(0); }
  virtual double getforml2norm2() { assert(0); }
  virtual double getFailprob(int,int) { assert(0); }
  virtual double getestat(int,int*,int,double&,int) { assert(0); }
  virtual double getPullIn(int) { assert(0); }
  virtual double getDCmass(double&,int*,int,double&,int) { assert(0); }
  
  // sensitivity analysis
  
  virtual void initOptInf (int) { assert(0); }
  virtual void setoptInf(int) { assert(0); }
  virtual void addOptInf(int,int,int,int) { assert(0); }
  
  virtual void getMidPointMass(double*,double**) { assert(0); }
  
  virtual void zerograd() { assert(0); }
  virtual void graddirect(int,double**,int n=0,int* nn=0) { assert(0); }
  virtual void gradcrit(int,double**,int,int*);// { assert(0); }
  virtual void gradadjoint(int, double**) { assert(0); }
  virtual void gradadjointStress(int, double**) { assert(0); }
  virtual void gradadjointModal (int, double**) { assert(0); }
  virtual void gradadjointGeom  (int, double**) { assert(0); }   
  
  virtual double getgradstrainenergy(int*,int,double&,int) { assert(0); }
  virtual double getgradmass(int*,int,int) { assert(0); }
  virtual double getgradMomOfInertia(int*,int,int,int) { assert(0); }
  virtual double getgradfrequency(int,int) { assert(0); }
  virtual double getgraddisp(int,int,int,double&,int) { assert(0); }
  virtual double getgradLambda(double&,int) { assert(0); }
  virtual double getgradinternalforce(int,int,double&,int) { assert(0); }
  virtual double getgradnodalpos(int,int,double&,int) { assert(0); }
  virtual double getgradnodstr(int,int,double&,int) { assert(0); }
  virtual double getgradkineticenergy(double&,int) { assert(0); }
  virtual double getgraddampingenergy(double&,int) { assert(0); }
  virtual double getgradaeroforce(int idir,double&,int) { assert(0); }
  virtual double getgradstressint(int,double&,double&,int,int,int*,int,double&,int) { assert(0); }
  virtual double getgraddisplevel(int,int,int,int*,int*,double*,double*,double&,
				  double&,int) { assert(0); }
  virtual double getgradstresslevel(int,int,int,int,int*,int*,double*,double*,
				    double&,double&,int) { assert(0); }
  virtual double getgradforml2norm2() { assert(0); }
  virtual double getgradFailprob(int,int,int=-1) { assert(0); }
  virtual double getgradestat(int,int*,int,double&,int) { assert(0); }
  virtual double getgradPullIn(int,int) { assert(0); }
  virtual double getgradDCmass(double&,int*,int,double&,int ){ assert(0); }
  
  virtual double getgradpartstrainenergy(int*,int,double&,int) { assert(0); }
  virtual double getgradpartinternalforce(int,int,double&,int) { assert(0); }
  virtual double getgradpartnodstr(int,int,double&,int) { assert(0); }
  virtual double getgradpartstressint(int,double&,double&,int,int,int*,int,double&,
				      int) { assert(0); }
  virtual double getgradpartdisplevel(int,int,int,int*,int*,double*,double*,double&,
				      double&,int);
  virtual double getgradpartstresslevel(int,int,int,int,int*,int*,double*,double*,
					double&,double&,int) { assert(0); }
  virtual double getgradpartestat(int,int*,int,double&,int) { assert(0); }
  virtual double getgradpartPullIn(int) { assert(0); }
  virtual double getgradpartDCmass(double&,int*,int,double&,int) { assert(0); }
  
  virtual void getgraddustrainenergy(int*,int,double&,int) { assert(0); }
  virtual void getgraddudisp(int,int,int,double&,int) { assert(0); }
  virtual void getgraddulambda(double&,int) { assert(0); }
  virtual void getgradduinternalforce(int,int,double&,int) { assert(0); }
  virtual void getgraddunodstr(int,int,double&,int) { assert(0); }
  virtual void getgraddumass(int*,int,int) { assert(0); }
  virtual void getgradduMomOfInertia(int*,int,int,int) { assert(0); }
  virtual void getgradduaeroforce(int,double&,int) { assert(0); }
  virtual void getgradduExtVal(int,int,double&) { assert(0); }
  virtual void getgraddustressint(int,double&,double&,int,int,int*,int,double&,int) { assert(0); }
  virtual void getgraddudisplevel(int,int,int,int*,int*,double*,double*,double&,
				  double&,int) { assert(0); }
  virtual void getgraddustresslevel(int,int,int,int,int*,int*,double*,double*,
				    double&,double&,int) { assert(0); }
  virtual void getgradduDCmass(double&,int*,int,double&,int) { assert(0); }
  
  // transient related functions (mostly for aeroelastic)  
  virtual void   addAnalysisData   ( anadata &) { assert(0); }
  virtual void   sndOptpar         ( int, int) { assert(0); }
  virtual void   setFluidvariables ( double*, int) { assert(0); }
  virtual void   setFluidcriteria  ( int, double*) { assert(0); }
  
  virtual void   setElecvariables  (double*,int,int) { assert(0); }
  virtual void   setThermvariables (double*,int,int) { assert(0); }
  
  virtual int    chkOptInfEL() { assert(0); }
  virtual int    chkOptInfTH() { assert(0); }
  
  virtual int    procComputeMaxSteps( double&) { assert(0); }
  virtual int    procSetvar         ( double, double) { assert(0); }
  virtual int    procTransition     ( int,    int) { assert(0); }
  virtual void   procEvaluate       ( double, double) { assert(0); }
  virtual void   procComputeStepSize( int) { assert(0); }
  virtual void   procSavevar() { assert(0); }
  
  virtual int    initProcess(double&, double&) { assert(0); }
  
  // output routines  
  virtual void postProcessing(int) { assert(0); }
};


// global singletons
extern std::auto_ptr<Structopt> structopt;
extern std::auto_ptr<Structopt> structrel;

extern std::auto_ptr<Optpro> optpro;
extern std::auto_ptr<Optpro> relpro;

#endif
#endif
