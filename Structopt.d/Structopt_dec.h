#ifndef _STRUCTOPT_DEC_HPP_
#define _STRUCTOPT_DEC_HPP_

#ifdef STRUCTOPT

#include <Utils.d/MyComplex.h>

#include <Structopt.d/Structopt_base.h>
#include <Structopt.d/Optcrit.h>
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optsol.h>

#include <Driver.d/DecDomain.h>

#include <Paral.d/MDStatic.h>
template<class Scalar> class GenMultiDomainStatic_opt;
template<class Scalar> class GenMultiDomainPostProcessor_opt;

#include <Driver.d/StaticProbType.h>
#include <Structopt.d/Driver_opt.d/StaticProbType_opt.h>

#include <Feti.d/DistrVector.h>
#include <Feti.d/Feti.h>
#include <Corotational.d/Corotator.h>
#include <Math.d/DistVector.h>
#include <Math.d/SparseMatrix.h>

#include <Utils.d/resize_array.h>

#include <Element.d/Element.h>

#define numBaseOptpar 5

template<class Scalar> class GenDecDomain_opt;


//------------------------------------------------------------------------------
template<class Scalar>
class Structopt_dec : public Structopt
{
 public:
  int numvar;
  int type; //0 - Interface for Design Optimization
            //1 - Interface for Reliability Analysis
  int numAnalysis;      
  int numFunc;
  int anagrdType;
  int optInfset;
  GenDecDomain_opt<Scalar> *structdom;
  Optpro *optpro;
  
  StcAnalysisData *analysisData;
  
  int sendInitDisp;
  
  OptActInfo** optInf;
  
  //ResizeArray<Stcvar*>   opv;        //Structural Variables  

  //for single domain static analysis
  int numStaticProb; 
  GenMultiDomainStatic_opt<Scalar> **staticpros;  
  StaticSolver_opt<Scalar, GenFetiSolver<Scalar>, GenDistrVector<Scalar>,
		   GenMultiDomainPostProcessor_opt<Scalar>, GenMultiDomainStatic_opt<Scalar>,
		   GenDistrVector<DComplex > >  **staticsolvers;
  
  // pointers to various solutions and fem arrays  
  int * analysisMap;
  int   activeAnalysis;
  
  //Scalar*  bcx;
  GenDistrVector<Scalar> * sol;
  GenDistrVector<Scalar> * vel;
  GenDistrVector<Scalar> * acc;
  
  GenDistrVector<Scalar> * grad;
  GenDistrVector<Scalar> * adj;
  
 public:  
  // Constructors  
  explicit Structopt_dec(int, GenDecDomain_opt<Scalar>* domain = 0, Optpro *_optpro = 0);
  ~Structopt_dec() {}

  void build(GenDecDomain_opt<Scalar>*, Optpro*);
  void initStructProbSol();
  void initAnalyticalSA();
  void getMidPointMass(double* totmas, double** midPoint);
  void func();
  void evaluate();
  int getAnalysisData(int anaId, double& time);
  int getAnalysisGradData(int anaId, double& time);
  int getAnalysisAdjData(int anaId, double& time, int zeroAll=0);
  void gradadjointGeom(int icrit, double** gc);

  double getstrainenergy(int*,int,double&,int);
  double getmass(int* eleList, int listSize, int anaId);
  void postProcessing(int giter);
  void gradadjoint(int icrit, double** gc);
  void gradadjointStress(int icrit, double** gc);
  void zerograd();
  double getgradstrainenergy(int* eleList, int listSize, double& time, int anaId);
  void getgraddustrainenergy(int* eleList, int listSize, double& time, int anaId);
  double getgradpartstrainenergy(int* eleList, int listSize, double& time, int anaId);
  double getgradmass(int* eleList, int listSize, int anaId);
  void cleanup();
  void initOptInf(int);
  void setoptInf(int);
  void addOptInf(int,int,int,int);
  void graddirect(int,double**,int n=0,int* nn=0); 

  double getdisp(int node, int typ, int dva, double& time, int anaId);
  double getgraddisp(int node, int typ, int dva, double& time, int anaId);
  void getgraddudisp(int,int,int,double&,int);

  double getdisplevel(int aFlag, int quoFlag, int size, 
		      int* nodeid, int* distyp, double* refVal, 
		      double* difVal, double& powFac, double& time, 
		      int anaId);
  void getgraddudisplevel(int aFlag, int quoFlag, int size, 
			  int* nodeid, int* distyp, double* refVal, 
			  double* difVal, double& powFac, double& time,
			  int anaId);
  double getgraddisplevel(int,int,int,int*,int*,double*,double*,double&,
			  double&,int);

  /*  
  void initReliabilityAnalysis();
  void cleanup();
  
  // set pointer to attributes, coordinates, nodal forces  
  double * getptreleattr    ( int& loc1, int& loc2 );
  double * getptrnodattr    ( int& loc1, int& loc2 );
  double * getptrnodcord    ( int& loc1, int& loc2 );
  double * getptrcomposite  ( int& loc1, int& loc2 );
  double * getptrnodalforce ( int& loc1, int& loc2 );
  double * getptrfluidattr  ( int& loc1);
  double * getptrelecattr   ( int& num,  int& loc1, int& loc2 );
  double * getptrthermattr  ( int &num, int& loc1, int& loc2 );
  
  // set pointer to attributes and coordinates  
  double * getptrgradeleattr    ( int& loc1, int& loc2 );
  double * getptrgradnodattr    ( int& loc1, int& loc2 );
  double * getptrgradnodcord    ( int& loc1, int& loc2 );
  double * getptrgradcomposite  ( int& loc1, int& loc2 );
  double * getptrgradnodalforce ( int& loc1, int& loc2 );
  double * getptrgradfluidattr  ( int& loc1);
  double * getptrgradelecattr   ( int& num,  int& loc2 );
  double * getptrgradthermattr  ( int &num, int& loc2 );
  
  
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
  
  void getMidPointMass(double*,double**);
  
  void zerograd();
  void graddirect(int,double**,int n=0,int* nn=0); 
  void gradcrit(int,double**,int,int*);  
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
  double getgradpartdisplevel(int,int,int,int*,int*,double*,double*,double&,
			      double&,int);
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
  
  
  // output routines  
  void postProcessing(int);

  */
};


#ifdef _TEMPLATE_FIX_
#include <Structopt.d/Structopt_dec.C>
#endif

#endif
#endif
