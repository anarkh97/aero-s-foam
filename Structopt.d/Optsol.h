#ifndef _OPTSOL_H_
#define _OPTSOL_H_

#ifdef STRUCTOPT

#include <Utils.d/resize_array.h>
#include <Structopt.d/Optinp.h>
#include <Structopt.d/Optgrad.h>

class Optsol;
class Optalg;
class Structopt;
class Optgrad;

struct graddata;

//------------------------------------------------------------------------------
//                              Optimization Solvers 
//                          Base Class and Derived Classes
//                                                  created  9/4/98 by Kurt
//------------------------------------------------------------------------------

class Optsol {

     public:

        FILE * optunitout;               //Outputfile

        char * optprotfile;              //Filename for solution protocol
	int    fsize;

        int iactalg;                     //Id number of current algorithm
        int numalg;                      //Number of Algorithms

        int giter;                       //Total Number of Iterations
	int liter;                       //Last Iterationnumber for output

        ResizeArray<Optalg*> optalg;     //Optimization Algorithms

        Optpro    *optpro;               //Current Optimization Problem
        Structopt *structopt;            //Current Structural Opt.Prob.

	int numvar;                      //Problem Characteristics
	int numcon;
	int numeqc;
	
	int *typcon;                     //Type of Constraint
	int *actcon;                     //Indicates active Constraints

	double *var;                     //Arrays for scaled values
	double *varup; 
	double *varlow;       	

        double obj;
	double *con;	
	double *gradobj;
        double *pgradcon;
	double **gradcon;

        double gnorm;                    // norm of gradients for SA of LSF
        double relSwitch;                // Switching reliability-index
 
     public:
     
        Optsol();

        int  getanagrad();

        double *  getVariables() { return var; }
        Optgrad * getActOptgrad();
	
	void solve(Optpro *optpro, Structopt *structopt);

        void addSolver( int &, int &,nlpdata &, graddata &);

	void initialize();
	void cleanup();
	            
	void printiter();
	
	void func(int iter);
	
        void resclvar();
	void scalfunc();
	void scalgrad(int ivar);

        double getGnorm() { return gnorm; }

};

//------------------------------------------------------------------------------

class Optalg {
    
     protected:
     
        int num;                         //Number of current Strategy
        int typ;                         //Type of current Srategy

        FILE* optunitout;                //Output file

        Optsol *optsol;                  //Current Optimization Strategy  

	Optgrad optgrad;                 //Gradient-Method

     public:
        virtual ~Optalg() {}
        static Optalg * create( int, int);

        virtual void setOutput(FILE*);
        virtual char* getSuffix(char*);

        virtual void ordercon(); 
        virtual void ordergradcon();         
        virtual void reorderActive(int*,int*);

        virtual int getGradtyp();

        virtual Optgrad* getOptgrad() { return &optgrad; }
	virtual Optsol* getOptsol() { return optsol; }
        
        virtual void buildalg( nlpdata & param)=0;
        virtual void setDefault()=0;

	virtual void solve(Optsol *_optsol)=0;
 	virtual void func(int iter)=0;
	virtual void grad(int* active=0)=0;
        virtual void print()=0;
        virtual void printres()=0;
};

//------------------------------------------------------------------------------

class OptalgNlpql : virtual public Optalg {
	
     public:
   
       int    subtype;         // 0: standard
                               // 1: for reliability 

       double acc;
       double scbou;

       int    maxfun;
       int    maxit;
       int    iprint;
       int    mode;
       int    ifail;
       int    merit;
       int    lql;
      	 
     public:

       OptalgNlpql(int t)  {subtype=t;}

       void buildalg( nlpdata & param);
       void setDefault();

       void solve(Optsol *_optsol);
       void func(int iter);
       void grad(int* active=0);
       void print();
       void printres();
      
};

//------------------------------------------------------------------------------

class OptalgNlpocm : public Optalg {
	
     public:
   
       double acc;
       double beta;
       double delta;	    
       double etha; 	    
       double xdgo;

       int    maxit;
       int    lines;
       int    ifail;	  

     public:
	
       void buildalg( nlpdata & param);
       void setDefault();

       void solve(Optsol *_optsol);
       void func(int iter);
       void grad(int* active=0);
       void print();
       void printres();
};

//------------------------------------------------------------------------------

class OptalgNlpslp : public Optalg {
	
     public:
       
       double qmove;
       double fmove;
       double fback;
       double tol;
       double objtol;
       double amijo;

       int imax;
       int mode;
       int incon;
       int iconv;
       int icycle;
       int iadapt;
       int ilise;
       int iprint;
       int ifail;

     public:
	
       void buildalg( nlpdata & param);
       void setDefault();

       void solve(Optsol *_optsol);
       void func(int iter);
       void grad(int* active=0);
       void print();
       void printres();
};

//------------------------------------------------------------------------------

class OptalgNlpmma : virtual public Optalg {
	
	int subtype;      // 0 - old MMA implementation
	                  // 1 - globally convergent version 
	
     public:

        int maxfun;
        int maxit;
        int isvan;
        int mixu;
        int mixl;
        int iappr;
        int itsub;
        int iprint;
        int ifail;

        double acc;
        double scbou;
        double alm;
        double sa;
        double sb;
        double sc;
        double dstep;
        double sau;
        double sbu;
        double sal;
        double sbl;
        double asscl;
        double assclu;
        double asscll;
        double fixup;
        double fixlow;

     public:
	
       OptalgNlpmma(int t) {subtype=t;}	
	
       void buildalg( nlpdata & param);
       void setDefault();

       void solve(Optsol *_optsol);
       void solveOld(Optsol *_optsol);
       void solveGC(Optsol *_optsol);
       
       void func(int iter);
       virtual void func(int iter, double *xval, double& f0val, double*fval);

       void grad(int* active=0);
       virtual void grad(double* xval, double* df0dx, double**dfdx, int* active=0);

       void print();
       void printres();
};

//------------------------------------------------------------------------------

class OptalgNlpsal : public Optalg {
	
     public:
   
       double acc;
       double accloc;
       double rpen;
       double penfac;
       double uinit;
       double bscal;
       double relbo;

       int    nbfgs;
       int    maxfun;
       int    maxit;
       int    iprint;
       int    ifail;

      	 
     public:
     
       void buildalg( nlpdata & param);
       void setDefault();

       void solve(Optsol *_optsol);
       void func(int iter);
       void grad(int* active=0);
       void print();
       void printres();
      
};

#endif

#endif
