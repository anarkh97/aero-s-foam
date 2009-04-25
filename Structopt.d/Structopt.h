#ifndef _STRUCTOPT_H_
#define _STRUCTOPT_H_

#include <Structopt.d/Optcrit.h>
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optsol.h>

#include <Driver.d/Domain.h>
#include <Driver.d/Dynam.h>
#include <Driver.d/StaticProbType.h>
#include <Driver.d/EigenProbType.h>
#include <Driver.d/DynamProbType.h>

#include <Problems.d/StaticDescr.h>
#include <Problems.d/EigenDescr.h>
#include <Problems.d/DynamDescr.h>

#include <Math.d/Vector.h>
#include <Math.d/VectorSet.h>
#include <Math.d/SparseMatrix.h>

#include <Utils.d/resize_array.h>

#include <Solvers.d/Solver.h>

#include <Element.d/Element.h>


class Domain;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
class Optcrit;
class Stcvar;
class Optpro;
class Optsol;

//------------------------------------------------------------------------------

struct StcAnalysisData {

   int dynFlag;
   int dynTyp;
   
   int dampFlag;
   int dampTyp;
   int dampFrq[2];
   
   int    transFlag;
   int    transMin;
   int    transMax;
   double designVel;
   double forceVel;
};

//------------------------------------------------------------------------------

class Structopt {

   public:

      int numcrit;
      int numvar;
      
      int numfunc;
      
      int anagrdType;

      double  *critval;
      double **critgard;
      
      Domain *structdom;
      Optpro *optpro;
      Optsol *optsol;  
      
      StcAnalysisData *analysisData;
      
      int aeroact;                       //Flag to save Sens. of Aeroforces
      double  *aeroforce;                //Aeroforces
      double **aerosens;                 //Sensitivity of Aeroforces

      ResizeArray<Optcrit*>  opc;        //Optimization Criteria
      ResizeArray<Stcvar*>   opv;        //Structural Variables     

                             //for single domain static analysis

      SingleDomainStatic<double,Vector,Solver> *staticpro;

      StaticSolver<double, Solver, Vector, 
                   SingleDomainPostProcessor<double, Vector, Solver>,
		   SingleDomainStatic<double,Vector,Solver>, ComplexVector > 
                                               *staticsolver;
 
                         //for single domain eigenvalue analysis
			 
      SingleDomainEigen *eigenpro;
	
      EigenSolver<DynamMat, Vector, VectorSet, 
	          SDEigenPostProcessor, SingleDomainEigen> 
		                               *eigensolver;

                         //for single domain eigenvalue analysis

      SingleDomainDynamic<double> *dynpro;
      
      DynamicSolver <DynamMat, Vector,
                     SDDynamPostProcessor, SingleDomainDynamic<double>,double>
		                               *dynsolver;

   public:
   
      // Constructors

      Structopt();

      void build (Domain *_domain, Optpro *_optpro, Optsol *_optsol);
      void cleanup();

      // set pointer to attributes, coordinates, nodal forces

      double * getptreleattr    ( int& loc1, int& loc2 );
      double * getptrnodattr    ( int& loc1, int& loc2 );
      double * getptrnodcord    ( int& loc1, int& loc2 );
      double * getptrcomposite  ( int& loc1, int& loc2 );
      double * getptrnodalforce ( int& loc1, int& loc2 );
      double * getptrfluidattr  ( int& loc1);

      // set pointer to attributes and coordinates

      double * getptrgradeleattr   ( int& loc1, int& loc2 );
      double * getptrgradnodattr   ( int& loc1, int& loc2 );
      double * getptrgradnodcord   ( int& loc1, int& loc2 );
      double * getptrgradcomposite ( int& loc1, int& loc2 );
      double * getptrgradfluidattr ( int& loc1);

      // get design criteria values

      void func();
      void evaluate();
      
      double getstrainenergy();
      double getmass();
      double getfrequency( int ieig );
      double getcontrolcost();
      double getdisp(int node, int typ, int dva);
      double getnodstr(int node, int typ);
      double getkineticenergy();
      double getdampingenergy();
      double getaeroforce(int idir);
      double getaeromoment(int node, int idir);

      // sensitivity analysis
        
      void zerograd();
      void graddirect(int); 
      void gradcrit();  
      void gradadjoint(int, double**);   

      double maxSensAeroforces(int);

      double getgradstrainenergy();
      double getgradmass();
      double getgradfrequency( int ieig );
      double getgraddisp(int node, int typ, int dva);
      double getgradnodstr(int node, int typ);
      double getgradkineticenergy();
      double getgraddampingenergy();
      double getgradaeroforce(int idir);
      double getgradaeromoment(int node, int idir);

      void getgraddustrainenergy();
      void getgraddudisp(int node, int typ, int dva);
      void getgraddunodstr(int node, int typ);
      void getgraddumass();
      void getgradduaeroforce();

      // transient related functions (mostly for aeroelastic)
      
      int     maxstep;
      double  stepcoef;
      
      void   addAnalysisData   ( anadata &);
      void   procEvaluate      ( double, double);
      void   sndOptpar         ( int, int);
      void   setFluidvariables ( double*, int);
      void   setFluidcriteria  ( int, double*);

      int    procSetvar ( double, double);

      double initProcess(int&, int&, double&, double&);

      // output routines
      
      void postProcessing(int);

};

#endif
