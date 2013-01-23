#ifndef _RELSOL_H_
#define _RELSOL_H_
#ifdef STRUCTOPT

#include <Utils.d/resize_array.h>
#include <Structopt.d/Optinp.h>

class Optpro;
class Structopt;
class Structopt_sd;
class Relsol;
class Absvar;
class Optsol;

struct nlpdata;
struct graddata;

struct RelsolData {

        int numlsf;
        int numdsgvar;

        double *  lsfunction;            //array for limit state function 
        double *  relIndex;              //array for reliability index
        double *  failureprob;           //array for failure probabilities
        double *  pmaValue;              //array for performance measures

        double **  gradlsfunction;       //array for gradients limit state function 
        double **  gradrelIndex;         //array for gradients of reliability
        double **  gradfailureprob;      //array for gradients of failure probabilities
        double **  gradpmaValue;         //array for gradients of performance measures

        RelsolData();

        void initialize(int);
        void initializeFORM();

        void initializeGrad(int);
        void initializeGradFORM();

        void cleanup();

        double getValue(int,int);
        double getGradient(int,int,int);

};


//------------------------------------------------------------------------------

class Relalg {

     public:
        virtual ~Relalg() {}

        enum {MonteCarlo,Form};

        int typ;                      //Type of Algorithm
 
        FILE* relunitout;             //Outputfile

        Relsol* relsol;               

        RelsolData* solData;

     public:

        int num;

        static Relalg * create( int, int, nlpdata&, graddata&);

        virtual int getGradtyp() { return 0; }

	virtual void solve(Relsol*, RelsolData*)=0;
 	virtual void func(int iter)=0;
        virtual void grad()=0;
        virtual void print()=0;
        virtual void printres()=0;
        virtual void setOutput(FILE*)=0; 
};

//------------------------------------------------------------------------------

class RelalgMonteCarlo: public Relalg {

   public:

     double samples;

     int * failure;

   public:

     RelalgMonteCarlo(nlpdata&);

     void setDefault();

     void solve(Relsol*,  RelsolData*);
     void func(int iter);
     void grad();
     void print();
     void printres();
     void setOutput(FILE* op) { relunitout=op; }
};

//------------------------------------------------------------------------------

class RelalgFORM: public Relalg {

   public:

     double** uvals;

     Optpro* optpro;          

   public:
     void initializeOptpro();  
     void cleanOptpro();

     void setLSFactive(int);   
     void setRIAproblem(int);   
     void setPMAproblem(int);   

     void resetLSFactive(int);   
     void resetRIAproblem(int);   
     void resetPMAproblem(int);   

     void saveResults(int);
     void saveRIAresults(int);
     void savePMAresults(int);

   public:

     RelalgFORM(int,nlpdata&,graddata&);

     int getGradtyp();

     void solve(Relsol*,  RelsolData*);
     void func(int iter);
     void grad();
     void print();
     void printres();
     void setOutput(FILE*);

     void analyticDSA(int);

     double normCdf(double);
};

//------------------------------------------------------------------------------

class Relsol {

     public:

        FILE*  relunitout;               //Output file

        char * relprotfile;              //Filename for solution protocol
	int    fsize;

        int numalg;                      //Number of Algorithms

        int giter;                       //Total Number of Iterations
	int liter;                       //Last Iterationnumber for output

        ResizeArray<Relalg*> relalg;     //Optimization Algorithms

        Optpro *relpro;                  //Current Reliability Analysis Problem
        Structopt_sd *structrel;            //Current Structural Prob.

	int numvar;                      //Problem Characteristics
	int numlsf;

        Absvar ** rndvar;                //random variable
        
        int dsgopt;                      //flag for design optimization
        int numdsgvar;                   //number of variables in design optimization

        RelsolData* solData;             //solution data
		
     public:
     
        Relsol();

	void solve( Optpro*, Structopt*);

        void addSolver( int &, int &,nlpdata &, graddata &);

	void initialize();
	void cleanup();
	            
	void printiter();
	
	void func(int iter);
        void scalfunc();
        void generateRandom();

        int getanagrad();

        double getFailprob(int);
        double getRelindex(int);
        double getPMAvalue(int);

        double getgradFailprob(int,int);
        double getgradRelindex(int,int);
        double getgradPMAvalue(int,int);
};

#endif
#endif
