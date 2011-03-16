#ifdef STRUCTOPT

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <Structopt.d/Optsol.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optcon.h>
#include <Structopt.d/Optobj.h>

#include <Structopt.d/Structopt_sd.h>

#include <Structopt.d/MMAgc.d/mma.h>

#include <Math.d/mathUtility.h>

#include <Utils.d/linkfc.h>
#include <Utils.d/DistHelper.h>


#ifdef WITH_IPOPT
#include <Structopt.d/Optsol_ipopt.h>
#endif

#ifdef DISTRIBUTED
#include <Structopt.d/Optsol_dist.h>
#include <Comm.d/Communicator.h>
#endif


//------------------------------------------------------------------------------
//                          Optimization Solution Procedure
//                                         and 
//                       Interfaces to Optimization Algorithms
//------------------------------------------------------------------------------

extern "C"      {

   void  _FORTRAN(nlpql) ( Optalg*, 
                           double*, 
			   double*, double*,
                           double&, double*, 
			   double*, double*, 
			   int&   , int&   ,int&,
                           double*, double*,double*,
                           double*, int*   ,int*, 
			   int&   , int&   ,int&,
			   int&   , int&   ,int&,
			   double&, double&,int&,
			   int&   , int&   ,int&,
			   int&   , int&   ,int&,
			   char*  , int&   ,int&);
}

//------------------------------------------------------------------------------

extern "C"      {

   void  _FORTRAN(relnlpql) ( Optalg*, 
                           double*, 
			   double*, double*,
                           double&, double*, 
			   double*, double*, 
			   int&   , int&   ,int&,
                           double*, double*,double*,
                           double*, int*   ,int*, 
			   int&   , int&   ,int&,
			   int&   , int&   ,int&,
			   double&, double&,int&,
			   int&   , int&   ,int&,
			   int&   , int&   ,int&,
			   char*  , int&   ,int&);
}

//------------------------------------------------------------------------------

extern "C"      {

   void  _FORTRAN(nlpocm) ( Optalg*, 
                            double*, 
			    double*, double*,
                            double*, double*, 
			    double*, double*, 
			    int&,    int&   ,int&,
			    double&, double&,double&,
			    double&, double&,int&,
			    int&,    char*,  int&,
                            int&);
}

//------------------------------------------------------------------------------

extern "C"      {

   void  _FORTRAN(nlpslp) ( Optalg*, 
                            double*, 
			    double*, double*,
                            double&, double*, 
			    double*, double*, 
			    double*, double*, double*,
			    double*, double*, 
                            int*,    int*,
			    int&,    int&   ,int&,
			    int&,    int&   ,int&,
			    int&,    int&   ,int&,
			    int&,    int&   ,char*,
                            int&,    double&,double&,
			    double&, double&,double&,
                            double&, int&);
}
//------------------------------------------------------------------------------

extern "C"      {

   void  _FORTRAN(nlpmma) ( Optalg*, 
                            double*, 
			    double*, double*,
                            double&, double*, 
			    double*, double*,double*,
			    double*, double*,double*,
                            int*,
			    int&,    int&   ,int&,
                            int&,    int&   ,int&,
			    int&,    int&   ,int&,
			    int&,    int&   ,int&,
			    int&,    int&   ,int&,
                            char*,   int&,
			    double&, double&,double&,
			    double&, double&,double&,
			    double&, double&,double&,
			    double&, double&,double&,
			    double&, double&,double&,
			    double&);
}

//------------------------------------------------------------------------------

extern "C"      {

   void  _FORTRAN(nlpsal) ( Optalg*, 
                            double*, 
			    double*, double*,
                            double&, double*, 
			    double*, double*,
			    double*, double*, double*,
                            double*, double*, double*, double*,
                            int*   , int*,    int*   , 
                            int*   ,
			    int&   , int&   , int&   ,
			    double&, double&, double&,
                            double&, double&, double&,
			    double&, int&   , int&   , 
                            int&   , int&   , int&   , 
                            int&   , char*  , int&);
}

//------------------------------------------------------------------------------

extern "C" {

   void _FORTRAN(func) (Optalg *optalg, int& iter) { 
       
       optalg->func(iter);
   }
}

//------------------------------------------------------------------------------

extern "C" {

   void _FORTRAN(grad) (Optalg *optalg, int *active) { 
       
       optalg->grad(active);
   }
}
   
//------------------------------------------------------------------------------

Optsol::Optsol():optalg(0) {

   numalg=0;

   giter=0;
   
   numvar=0;
   numcon=0;
   numeqc=0;
}


//------------------------------------------------------------------------------

void Optsol::solve (Optpro *_optpro, Structopt *_structopt) {


   optpro=_optpro;
   structopt=_structopt;

   optunitout=optpro->optunitout;

   initialize();

   for (iactalg=0;iactalg<numalg;iactalg++) {

      giter=0;
      liter=-1;

      optalg[iactalg]->getOptgrad()->setactpro(optpro,this,structopt);
      optalg[iactalg]->setOutput(optunitout);
      optalg[iactalg]->print();
      optalg[iactalg]->solve(this);
   }
  
   cleanup();
}  

//------------------------------------------------------------------------------

void Optsol::initialize() {

   int i,j;

   numvar=optpro->optvar->nAbs();
   
   var    = new double[numvar+1];
   varup  = new double[numvar+1];
   varlow = new double[numvar+1];
   
   for (i=0;i<numvar;i++) {
     var[i]   = optpro->optvar->getAbsVar(i)->getScaledValue();
     varup[i] = optpro->optvar->getAbsVar(i)->getScaledUpperBound();
     varlow[i]= optpro->optvar->getAbsVar(i)->getScaledLowerBound();
   }
   
   numcon=optpro->optcon->numcon;
   numeqc=optpro->optcon->numeqc;
      
   // allocate memory
      
   con      = new double[numcon];
   gradobj  = new double[numvar+1];
   typcon   = new int[numcon];
   actcon   = new int[numcon];
   
   pgradcon = new double [(numvar+1)*numcon];
   gradcon  = new double*[numvar+1];
   for (i=0;i<(numvar+1);i++) gradcon[i]=pgradcon+i*numcon; 
      
   // initialize arrays
   
   for (i=0;i<(numvar+1);i++) {
     gradobj[i] = 0.0;
     for (j=0;j<numcon;j++) gradcon[i][j] = 0.0; 
   }
   
   // Watch: by default all constraints are set active
   
   for (i=0;i<numcon;i++) {
     con[i]    = 0.0;
     actcon[i] = 1;
     typcon[i] = static_cast<int>(optpro->optcon->typcon[i]);
   }  
}       

//------------------------------------------------------------------------------

void Optsol::cleanup() {

   delete [] var;
   delete [] varup;
   delete [] varlow;

   delete [] con;
   delete [] gradobj;
   delete [] typcon;
   delete [] actcon;

   delete [] pgradcon;
   delete [] gradcon;

   var     = 0;
   varup   = 0;
   varlow  = 0;

   con     = 0;
   gradobj = 0;
   typcon  = 0;
   actcon  = 0;

   pgradcon = 0;
   gradcon  = 0;
   
   structopt->cleanup();   
}

//------------------------------------------------------------------------------

void Optsol::resclvar() {

   int i;

   for (i=0;i<numvar;i++) 
     optpro->optvar->getAbsVar(i)->putScaledValue(var[i]);
}

//------------------------------------------------------------------------------

void Optsol::scalfunc() {
   
   obj = optpro->optobj->getScaledValue();

   int i;
   for (i=0;i<optpro->optcon->numcon;i++)
      con[i] = optpro->optcon->getScaledValue(i);
}

//------------------------------------------------------------------------------

void Optsol::scalgrad(int ivar) {
   
   gradobj[ivar] = optpro->optobj->getScaledGradient();

   // update only active criteria

   int i;
   for (i=0;i<optpro->optcon->numcon;i++) 
     if (actcon[i]) 
       gradcon[ivar][i] = optpro->optcon->getScaledGradient(i);
}

//------------------------------------------------------------------------------

void Optsol::func(int iter) {

   resclvar();                                  //Rescaling of Opt. Variables

   optpro->optvar->updvar();                    //Update Structural Variables
   
   structopt->func();                           //Evaluation of Opt.Crit.

   optpro->optobj->func();                      //Evaluation of Objective
  
   optpro->optcon->func();                      //Evaluation of Constraints

   scalfunc();                                  //Scaling of Objective & Const.

   optpro->optvar->resetvar();                  //Reset Structural Variables
   
   if ( liter < iter ) {
     giter+=1;
     liter=iter;
     filePrint(stderr,"\n ... Optimization Step Nr. %d  (%d)   ...\n",giter,iactalg); 
     structopt->postProcessing(giter);
     printiter();
   }
}

//------------------------------------------------------------------------------

void Optsol::printiter() {

   filePrint(optunitout,"\n\n\t%d. Iteration in Algorithm %d:\n\n",giter,iactalg);
   
   optpro->optobj->printres(optunitout);
  
   optpro->optcon->printres(optunitout);

   optpro->optvar->printres(optunitout);

   fflush(optunitout);
}

//------------------------------------------------------------------------------

int Optsol::getanagrad () 
{
   int i;
   int anagrad=0;

   // Watch: we get "0" form optalg if analytical SA

   for (i=0;i<numalg;i++) { if ( ! optalg[i]->getGradtyp() ) anagrad=1; }
   
   // return "1" if SA is needed

   return anagrad;
}  

//------------------------------------------------------------------------------

Optgrad * Optsol::getActOptgrad()
{
   return optalg[iactalg]->getOptgrad();
}
 
#ifdef DISTRIBUTED
//------------------------------------------------------------------------------
Optalg * Optalg::create( int _num, int _typ) 
{
  Optalg * optalg;
  if(structCom->myID() == Optalg_dist::root)
    {
      switch (_typ) 
	{	    
	case 0:
	  optalg = new OptalgNlpql_dist(0, *(structCom->getCommunicator()));
	  break;
	case 3:
	  optalg = new OptalgNlpmma_dist(0, *(structCom->getCommunicator()));
	  break;
	case 5:
	  optalg = new OptalgNlpmma_dist(1, *(structCom->getCommunicator()));
	  break;
	default:
	  assert(0);
	}
    }
  else
    {
      optalg = new Optalg_dist_slaveFacade(*(structCom->getCommunicator()));
    }
  optalg->num=_num;
  optalg->typ=_typ;
  return optalg;
}   

#else
//------------------------------------------------------------------------------
Optalg * Optalg::create( int _num, int _typ) 
{
  Optalg * optalg;    
  switch (_typ) 
    {	    
    case 0:
      optalg = new OptalgNlpql(0);          
      break;
    case 1:
      optalg = new OptalgNlpocm();
      break;
    case 2:
      optalg = new OptalgNlpslp();
      break;
    case 3:
      optalg = new OptalgNlpmma(0);
      break;
    case 4:
      optalg = new OptalgNlpsal();
      break;
    case 5:
      optalg = new OptalgNlpmma(1);
      break;
#ifdef WITH_IPOPT
    case 6:
      optalg = new OptalgIpopt();
      break;
#endif
    case 10:
      optalg = new OptalgNlpql(1);          
      break;
    default:
      fprintf(stderr," *** ERROR: Optimization algorithm %d unkown\n", _typ);
      exit(-1);
    }
  optalg->num=_num;
  optalg->typ=_typ;
  return optalg;
}   
#endif

//------------------------------------------------------------------------------

void Optalg::setOutput(FILE* op) 
{
   optunitout=op; 
}

//------------------------------------------------------------------------------

int Optalg::getGradtyp()
{
   return optgrad.typ;
}

//------------------------------------------------------------------------------

void Optalg::ordercon() {

   // order constraints such that  1. equality constraints   (typ=0)
   //                              2. inequality constraints (typ=1)   

   int i;
   int ieq=0;
   int eqc=0;
   int numeqc=optsol->numeqc;
   int numcon=optsol->numcon;

   double *xcon = new double[numcon];
    
   for (i=0;i<numcon;i++) {

     if (optsol->typcon[i] == 0) {
       xcon[eqc] = optsol->con[i] ; eqc++;  }
     else { 
       xcon[ieq+numeqc] = optsol->con[i] ; ieq++; }
   }
   
   for (i=0;i<numcon;i++) { optsol->con[i]=xcon[i] ; }
   
   delete [] xcon;
}

//------------------------------------------------------------------------------

void Optalg::ordergradcon() {

   // order constraints such that  1. equality constraints   (typ=0)
   //                              2. inequality constraints (typ=1)   

   int numcon=optsol->numcon;
   int numvar=optsol->numvar;
   int numeqc=optsol->numeqc;

   double *xcon = new double[numcon];
   
   int i,j,ieq,eqc;
   
   for (j=0;j<numvar;j++) {
      
      ieq=0;
      eqc=0;

      for (i=0;i<numcon;i++) {

         if (optsol->typcon[i] == 0) {
           xcon[eqc] = optsol->gradcon[j][i] ; eqc++; }
         else { 
           xcon[ieq+numeqc] = optsol->gradcon[j][i] ;ieq++; }
      }   

      for (i=0;i<numcon;i++) { optsol->gradcon[j][i]=xcon[i] ; }
   }
   
   delete [] xcon;
}

//------------------------------------------------------------------------------

void Optalg::reorderActive(int* active, int* actcon) {

   // set of active constraints is ordered such that  1. equality constraints   (typ=0)
   //                                                 2. inequality constraints (typ=1)   

   int i;
   int ieq=0;
   int eqc=0;

   if(!active) return;

   int numcon=optsol->numcon;
   int numeqc=optsol->numeqc;
  
   for (i=0;i<numcon;i++) {

     if (optsol->typcon[i] == 0) {
       actcon[i] = active[eqc];        eqc++; }
     else { 
       actcon[i] = active[ieq+numeqc]; ieq++; }
   }
}

//------------------------------------------------------------------------------

char * Optalg::getSuffix(char * fname)
{
   char * suffix;

   int    ssize  = strlen(fname);
   char * fpoint = strstr(fname,".");
   int    psize  = fpoint - fname;
   
   if ( psize < ssize ) 
   { 
     suffix  = new char[ssize-psize];
     int i;
     for (i=psize+1;i<ssize;i++) suffix[i-psize-1]=fname[i];
     suffix[ssize-psize-1] = 0;
   } 
   else
   {
     suffix = "";
   }
   
   return suffix;
}

//------------------------------------------------------------------------------

void OptalgNlpql::setDefault()
{
    acc    = 1.0e-4;
    scbou  = 1.0e+6;

    maxit  = 10;
    maxfun = 10;
    iprint =  2;
    mode   =  0;
    merit  =  1;
    lql    =  1;
}

//------------------------------------------------------------------------------

void OptalgNlpql::buildalg( nlpdata & param) {

    setDefault();

    if (param.rflag[0]) acc    = param.rval[0];
    if (param.rflag[1]) scbou  = param.rval[1];

    if (param.iflag[0]) maxit  = param.ival[0];
    if (param.iflag[1]) maxfun = param.ival[1];
    if (param.iflag[2]) iprint = param.ival[2];
    if (param.iflag[3]) merit  = param.ival[3];
    if (param.iflag[4]) lql    = param.ival[4];

    ifail  = 0;
}

//------------------------------------------------------------------------------

void OptalgNlpql::solve(Optsol *_optsol) {

    optsol =_optsol;    

    char * outfile = optsol->optprotfile;
    int    fsize   = optsol->fsize;
    int    iout=45;
    
    // check if optimizer used within optimization

    char * suffix = getSuffix(outfile);

    if ( !strcmp(suffix,"rpo") ) iout=47;

    // allocate memory for algorithm

    int n =  optsol->numvar;
    int m =  optsol->numcon;

    int mmax = max(1,m);
    int nmax = max(2,n+1);

    int mnn2 = m + n + n + 2;
    int n2c  = nmax*nmax;
    int lwa    = mmax*n + 4*mmax + 4*m + 18*n + 55
               + 3 *(n+1)*(n+1) + 20*n + 4*m + 20;
    int lkwa   = 19 + m+n+n;
    int lactiv = 2*mmax + 15;

    double rmem = sizeof(double)*(mnn2+n2c+nmax+lwa) + sizeof(int)*(lkwa+lactiv);

    filePrint(stderr," ... Allocating memory for NLPQL: %8.2f Mb.\n",rmem/1024.0/1024.0);    

    double* up = new double[mnn2];
    double* cp = new double[n2c];
    double* dp = new double[nmax];
    double* wp = new double[lwa];
    int*    kp = new int[lkwa];
    int*    ap = new int[lactiv];

    // call sqp algorithm
    
    switch (subtype)
    { 
      case 0:
        _FORTRAN(nlpql) ( this,
   	    		  optsol->var	 ,
   	    		  optsol->varlow , optsol->varup   ,
        		  optsol->obj	 , optsol->con     ,
        		  optsol->gradobj, optsol->pgradcon,
        		  optsol->numvar , optsol->numcon  , optsol->numeqc,
        		  up		 , cp		   , dp 	   ,
        		  wp		 , kp		   , ap 	   ,
        		  mmax  	 , nmax 	   , mnn2	   ,
        		  lwa		 , lkwa 	   , lactiv	   ,
    	    		  acc		 , scbou	   , maxfun	   ,
   	    		  maxit 	 , iprint	   , mode	   ,
   	    		  ifail 	 , merit	   , lql	   ,
   	    		  outfile	 , fsize	   , iout);
        break;
      case 1:
     	_FORTRAN(relnlpql) ( this,
 	    		  optsol->var	 ,
 	    		  optsol->varlow , optsol->varup   ,
     			  optsol->obj	 , optsol->con     ,
     			  optsol->gradobj, optsol->pgradcon,
     			  optsol->numvar , optsol->numcon  , optsol->numeqc,
     			  up		 , cp		   , dp 	   ,
     			  wp		 , kp		   , ap 	   ,
     			  mmax  	 , nmax 	   , mnn2	   ,
     			  lwa		 , lkwa 	   , lactiv	   ,
     			  acc		 , scbou	   , maxfun	   ,
 	    		  maxit 	 , iprint	   , mode	   ,
 	    		  ifail 	 , merit	   , lql	   ,
 	    		  outfile	 , fsize	   , iout);

       // compute norm of gradient of LSF w.r.t. u

       double gnorm = 0.0;
       int i;
       for (i=0;i<optsol->numvar;i++) {
         double gi= optsol->gradcon[0][i];
         gnorm += gi*gi;
       }
       optsol->gnorm = sqrt(gnorm);

       // check wether mean value design in feasible
     
       optsol->relSwitch = 1.0;

       if (up[0] >  0.0 ) optsol->relSwitch = -1.0;

       break;
    }
		      
    printres();

    // free memory used in algorithm

    delete [] up;
    delete [] cp;
    delete [] dp;
    delete [] wp;
    delete [] kp;
    delete [] ap;
}

//------------------------------------------------------------------------------

void OptalgNlpql::func(int iter) {

   optsol->func(iter);
   
   ordercon();
}

//------------------------------------------------------------------------------

void OptalgNlpql::grad(int* active) {

   // skip update and reset of variables in case of numerical SA
   // as it leads to errors for reliability based analysis

   int gtyp = optgrad.typ; // 0=analytical, 1=numerical

   reorderActive(active,optsol->actcon);   

   if (!gtyp) optsol->optpro->optvar->updvar();      

   optgrad.grad(optsol->actcon);

   if (!gtyp) optsol->optpro->optvar->resetvar();    
   
   ordergradcon();
}

//------------------------------------------------------------------------------

void OptalgNlpql::print() {


   fprintf(optunitout,"\n\tSolution Strateg #%d. : nlpql\n",num);
   fprintf(optunitout,  "\t=============================\n\n");

   fprintf(optunitout,"\tAccuracy ....................... : %12.5f\n",acc);
   fprintf(optunitout,"\tScaling Bound .................. : %12.5f\n\n",scbou);

   fprintf(optunitout,"\tMaximum Number of Function Calls : %d\n",maxfun);
   fprintf(optunitout,"\tMaximum Number of Iterations ... : %d\n\n",maxit);
   
   fprintf(optunitout,"\tPrint Mode ..................... : %d\n",iprint);
   fprintf(optunitout,"\tExecution Mode ................. : %d\n",mode);
   fprintf(optunitout,"\tType of Merit Funtion .......... : %d\n",merit);
   fprintf(optunitout,"\tLQL Type ....................... : %d\n",lql);

   optgrad.print(optunitout);
}

//------------------------------------------------------------------------------

void OptalgNlpql::printres() {

   fprintf(optunitout,"\n\tResult of %d. Strategy: nlpql\n\n",num);
   
   switch (ifail) {
   
     case 0:
        fprintf(optunitout,
	"\tTHE OPTIMALITY CONDITIONS ARE SATISFIED.\n");
	break;
     case 1:
	fprintf(optunitout,
	"\tTHE ALGORITHM HAS BEEN STOPPED AFTER MAXIT ITERATIONS.\n");
	break;
     case 2:
	fprintf(optunitout,
        "\tTHE ALGORITHM COMPUTED AN UPHILL SEARCH DIRECTION.\n");
	break;
     case 3:
	fprintf(optunitout,
        "\tUNDERFLOW OCCURED WHEN DETERMINING A NEW APPROXIMATION MATRIX \n");
	fprintf(optunitout,
	"\tFOR THE HESSIAN OF THE LAGRANGIAN.\n");
	break;
     case 4:
	fprintf(optunitout,
        "\tMORE THAN MAXFUN FUNCTION EVALUATIONS ARE REQUIRED\n");
	fprintf(optunitout,
	"\tDURING THE LINE SEARCH ALGORITHM.\n");
	break;
     case 5:
	fprintf(optunitout,
        "\tLENGTH OF A WORKING ARRAY IS TOO SHORT. MORE DETAILED ERROR\n");
	fprintf(optunitout,
	"\tINFORMATION IS OBTAINED WITH IPRINT .GT. 0 .\n");
	break;
     case 6:
	fprintf(optunitout,
        "\tTHERE ARE FALSE DIMENSIONS, I.E.  M .GT. MMAX ,\n");
	fprintf(optunitout,
	"\tN .GE. NMAX , OR  MNN2 .NE. M+N+N+2 .\n");
	break;
     case 7:
	fprintf(optunitout,
        "\tTHE SEARCH DIRECTION IS CLOSE TO ZERO, BUT THE\n");
	fprintf(optunitout,
	"\tCURRENT ITERATE IS STILL INFEASIBLE.\n");
	break;
     default:
	fprintf(optunitout,
        "\tError Message not specified.\n");
   }
}     

//------------------------------------------------------------------------------

void OptalgNlpocm::setDefault()
{
    acc    = 1.0e-4;
    beta   = 0.7;
    delta  = 0.1;
    etha   = -1;
    xdgo   = 0;

    maxit  = 10;
    lines  = 1;
}

//------------------------------------------------------------------------------

void OptalgNlpocm::buildalg( nlpdata & param) {

    setDefault();

    if (param.rflag[0]) acc    = param.rval[0];
    if (param.rflag[1]) beta   = param.rval[1];
    if (param.rflag[2]) delta  = param.rval[2];
    if (param.rflag[3]) etha   = param.rval[3];
    if (param.rflag[4]) xdgo   = param.rval[4];

    if (param.iflag[0]) maxit  = param.ival[0];
    if (param.iflag[1]) lines  = param.ival[1];
}

//------------------------------------------------------------------------------

void OptalgNlpocm::solve(Optsol *_optsol) {

    optsol =_optsol;    

    char * outfile = optsol->optprotfile;
    int    fsize   = optsol->fsize;
    
    _FORTRAN(nlpocm) ( this,
		      optsol->var,
		      optsol->varlow, optsol->varup,
                      &(optsol->obj), optsol->con,
                      optsol->gradobj,optsol->pgradcon,
                      optsol->numvar, optsol->numcon,    optsol->numeqc,
		      acc,            beta,              delta,
		      etha,           xdgo,  	         maxit, 
		      lines,          outfile,           fsize,
                      ifail);

    printres();
}

//------------------------------------------------------------------------------

void OptalgNlpocm::func(int iter) {

   optsol->func(iter);
}

//------------------------------------------------------------------------------

void OptalgNlpocm::grad(int* active) {

   optsol->optpro->optvar->updvar();      

   optgrad.grad();

   optsol->optpro->optvar->resetvar();    
}

//------------------------------------------------------------------------------

void OptalgNlpocm::print() {


   fprintf(optunitout,"\n\tSolution Strateg #%d. : nlpocm\n",num);
   fprintf(optunitout,  "\t==============================\n\n");

   fprintf(optunitout,"\tAccuracy ....................... : %12.5f\n",acc);
   fprintf(optunitout,"\tBeta exponent................... : %12.5f\n",beta);
   fprintf(optunitout,"\tStep size....................... : %12.5f\n",delta);

   fprintf(optunitout,"\tGradient shift.................. : %12.5f\n",etha);
   fprintf(optunitout,"\tEquality constraint correction.. : %12.5f\n",xdgo);
   							     
   fprintf(optunitout,"\tMaximum number of iterations.... : %12d\n",maxit);
   fprintf(optunitout,"\tLine search..................... : %12d\n",lines);

   optgrad.print(optunitout);
}

//------------------------------------------------------------------------------

void OptalgNlpocm::printres() {


   fprintf(optunitout,"\n\tResult of %d. Strategy: nlpocm\n\n",num);
   
   switch (ifail) {
   
     case 0:
        fprintf(optunitout,
	"\tTHE ALGORITHM HAS CONVERGED.\n");
	break;
     case 1:
	fprintf(optunitout,
	"\tTHE ALGORITHM HAS BEEN STOPPED AFTER MAXIT ITERATIONS.\n");
	break;
     default:
	fprintf(optunitout,
        "\tError Message not specified.\n");
   }
}     

//------------------------------------------------------------------------------

void OptalgNlpslp::setDefault()
{
    qmove  = 0; 
    fmove  = 0;
    fback  = 0;
    tol    = 0;
    objtol = 0;
    amijo  = 0;

    imax   = 0;
    incon  = 0;
    iconv  = 0;
    icycle = 0;
    iadapt = 0;
    ilise  = 0;
    iprint = 0;
}
//------------------------------------------------------------------------------

void OptalgNlpslp::buildalg( nlpdata & param) 
{
    setDefault();

    if (param.rflag[0]) qmove  = param.rval[0];
    if (param.rflag[1]) fmove  = param.rval[1];
    if (param.rflag[2]) fback  = param.rval[2];
    if (param.rflag[3]) tol    = param.rval[3];
    if (param.rflag[4]) objtol = param.rval[4];
    if (param.rflag[5]) amijo  = param.rval[5];

    if (param.iflag[0]) imax   = param.ival[0];
    if (param.iflag[1]) incon  = param.ival[1];
    if (param.iflag[2]) iconv  = param.ival[2];
    if (param.iflag[3]) icycle = param.ival[3];
    if (param.iflag[4]) iadapt = param.ival[4];
    if (param.iflag[5]) ilise  = param.ival[5];
    if (param.iflag[6]) iprint = param.ival[6];
    
    mode   = 0;
}

//------------------------------------------------------------------------------

void OptalgNlpslp::solve(Optsol *_optsol) {

    optsol =_optsol;    

    char * outfile = optsol->optprotfile;
    int    fsize   = optsol->fsize;

    // allocate memory for algorithm

    int n =  optsol->numvar;
    int m =  optsol->numcon;

    int mmax = max(1,m);

    int mnn    = m +2*n;
    int ny     = m+2*n+2;
    int lwa    = ((mnn+3)*(n+2) + 5*n + mnn + 1)*2 + (mnn + n + 5)*1;
 
    double rmem = sizeof(double)*(3*n+mmax+ny) + sizeof(int)*(lwa+mmax);

    filePrint(stderr," ... Allocating memory for SLP: %8.2f Mb.\n",rmem/1024.0/1024.0);    

    double* dmin   = new double[n];
    double* dmax   = new double[n];
    double* scx    = new double[n];
    double* scg    = new double[mmax];
    double* ypp    = new double[ny];

    int* iwa       = new int[lwa];
    int* active    = new int[mmax];
    
     _FORTRAN(nlpslp) ( this,
     			optsol->var,
     			optsol->varlow, optsol->varup,
     			optsol->obj,	optsol->con,
     			optsol->gradobj,optsol->pgradcon,
     			dmin,		dmax,		  scx,
     			scg,		ypp,
     			iwa,		active,
     			optsol->numvar, optsol->numcon,   optsol->numeqc,
     			imax,		mode,		  incon,
     			iconv,  	icycle, 	  iadapt,
     			ilise,  	iprint, 	  outfile,
     			fsize,  	qmove,  	  fmove,
     			fback,  	tol,		  objtol,
     			amijo,  	ifail);


    delete [] dmin   ;
    delete [] dmax   ;
    delete [] scx    ;
    delete [] scg    ;
    delete [] active ;
    delete [] ypp    ;
    delete [] iwa    ; 

    printres();
}

//------------------------------------------------------------------------------

void OptalgNlpslp::func(int iter) {

   optsol->func(iter);

   ordercon();
}

//------------------------------------------------------------------------------

void OptalgNlpslp::grad(int* active) {

   reorderActive(active,optsol->actcon);   

   optsol->optpro->optvar->updvar();      

   optgrad.grad(optsol->actcon);

   optsol->optpro->optvar->resetvar();    

   ordergradcon();
}

//------------------------------------------------------------------------------

void OptalgNlpslp::print() {


   fprintf(optunitout,"\n\tSolution Strateg #%d. : nlpslp\n",num);
   fprintf(optunitout,  "\t==============================\n\n");

   optgrad.print(optunitout);
}

//------------------------------------------------------------------------------

void OptalgNlpslp::printres() {


   fprintf(optunitout,"\n\tResult of %d. Strategy: nlpslp\n\n",num);
   
   switch (ifail) {
   
     case 0:
        fprintf(optunitout,
	"\tTHE ALGORITHM HAS CONVERGED.\n");
	break;
     case 1:
	fprintf(optunitout,
	"\tTHE ALGORITHM HAS BEEN STOPPED AFTER MAXIT ITERATIONS.\n");
	break;
     default:
	fprintf(optunitout,
        "\tError Message not specified.\n");
   }
}     
//------------------------------------------------------------------------------

void OptalgNlpmma::setDefault()
{
    maxfun = 20;
    maxit  = 500;
    isvan  = 1;
    mixu   = 1;
    mixl   = 1;
    iappr  = 3;
    itsub  = 20;
    iprint = 1;

    acc    = 1.0e-1;
    scbou  = 0.0e+0;
    alm    = 1.0e+0;
    sa     = 0.5;
    sb     = 0.7;
    sc     =-1.0;
    dstep  = 1.0;
    sau    = 0.5;
    sbu    = 0.7;
    sal    = 0.5;
    sbl    = 0.7;
    asscl  = 0.5;
    assclu = 0.5;
    asscll = 0.5;
    fixup  = 1.0;
    fixlow = 0.01;
}

//------------------------------------------------------------------------------

void OptalgNlpmma::buildalg( nlpdata & param) 
{
    setDefault();

    if (param.iflag[0]) maxfun = param.ival[0];
    if (param.iflag[1]) maxit  = param.ival[1];
    if (param.iflag[2]) isvan  = param.ival[2];
    if (param.iflag[3]) mixu   = param.ival[3];
    if (param.iflag[4]) mixl   = param.ival[4];
    if (param.iflag[5]) iappr  = param.ival[5];
    if (param.iflag[6]) itsub  = param.ival[6];
    if (param.iflag[7]) iprint = param.ival[7];

    if (param.rflag[0])  acc	= param.rval[0];
    if (param.rflag[1])  scbou  = param.rval[1];
    if (param.rflag[2])  alm	= param.rval[2];
    if (param.rflag[3])  sa	= param.rval[3];
    if (param.rflag[4])  sb	= param.rval[4];
    if (param.rflag[5])  sc	= param.rval[5];
    if (param.rflag[6])  dstep	= param.rval[6];
    if (param.rflag[7])  sau	= param.rval[7];
    if (param.rflag[8])  sbu	= param.rval[8];
    if (param.rflag[9])  sal	= param.rval[9];
    if (param.rflag[10]) sbl	= param.rval[10];
    if (param.rflag[11]) asscl  = param.rval[11];
    if (param.rflag[12]) assclu = param.rval[12];
    if (param.rflag[13]) asscll = param.rval[13];
    if (param.rflag[14]) fixup  = param.rval[14];
    if (param.rflag[15]) fixlow = param.rval[15];

    if ( sc < 0 )  sc = 1.0/sb;
}

//------------------------------------------------------------------------------

void OptalgNlpmma::solve(Optsol *_optsol) {

   if (subtype) 
     solveGC(_optsol);
   else
     solveOld(_optsol);

}

//------------------------------------------------------------------------------

void OptalgNlpmma::solveOld(Optsol *_optsol) {

    optsol =_optsol;    

    char * outfile = optsol->optprotfile;
    int    fsize   = optsol->fsize;

    // allocate memory for algorithm

    int n =  optsol->numvar;
    int m =  optsol->numcon;

    int mmax = max(1,m);
    int nmax = max(1,n);

    int lwa  = 6*mmax+8*nmax+nmax*mmax;

    double rmem = sizeof(double)*(lwa+mmax+nmax+mmax) + sizeof(int)*mmax ;

    filePrint(stderr," ... Allocating memory for MMA: %8.2f Mb.\n",rmem/1024.0/1024.0);    

    double* up   = new double[mmax];
    double* scxp = new double[nmax];
    double* scgp = new double[mmax];
    double* wap  = new double[lwa];
    int*    actp = new int[mmax];
    
     _FORTRAN(nlpmma) ( this,
	   	        optsol->var,
 		        optsol->varlow, optsol->varup,
                        optsol->obj,    optsol->con,
                        optsol->gradobj,optsol->pgradcon, up,
                        scxp,           scgp,             wap,
                        actp,
                        optsol->numvar, optsol->numcon,   optsol->numeqc,
                        nmax,           mmax,             lwa,            
     			maxfun,         maxit,            isvan,
                        mixu,           mixl,             iappr,
                        itsub,          iprint,           ifail,
                        outfile,        fsize,
     			acc,            scbou,            alm,
                        sau,            sbu,              sal,
                        sbl,            asscl,            assclu,
                        asscll,         fixup,            fixlow,
                        sa,             sb,               sc,
			dstep);

    printres();

    delete [] up;
    delete [] scxp;
    delete [] scgp;
    delete [] wap;
    delete [] actp;

}
//------------------------------------------------------------------------------

void OptalgNlpmma::solveGC(Optsol *_optsol) {

    optsol =_optsol;    
    
    FILE* optprot = fopen(optsol->optprotfile,"a");

    MMAgc mmaAlg(this,
                 optsol->var, optsol->varup, optsol->varlow,
		 optsol->numvar, optsol->numcon, 
		 maxit,itsub,acc,sa,sb,sc,dstep,optprot);
    
    ifail = mmaAlg.solve();

    mmaAlg.cleanup();

    fclose(optprot);

    printres();

}

//------------------------------------------------------------------------------

void OptalgNlpmma::func(int iter) 
{
  optsol->func(iter);
  ordercon();
  return;
}

//------------------------------------------------------------------------------

void OptalgNlpmma::func(int iter, double *xval, double& f0val, double*fval)
{
   // update variables
  for(int i=0;i<optsol->numvar;i++) { optsol->var[i] = xval[i]; }
  func(iter);
  // save results 
  // contraints are multiplied by "-1" due to contraint definition
  f0val = optsol->obj;  
  for(int j=0;j<optsol->numcon;j++) { fval[j] = -1.0 * optsol->con[j]; }
  return;
}

//------------------------------------------------------------------------------

void OptalgNlpmma::grad(int* active) 
{
  reorderActive(active,optsol->actcon);
  optsol->optpro->optvar->updvar();      
  optgrad.grad(optsol->actcon);
  optsol->optpro->optvar->resetvar();
  ordergradcon();
  return;
}

//------------------------------------------------------------------------------

void OptalgNlpmma::grad(double* xval, double* df0dx, double** dfdx, int* active) 
{
  // update variables (should not be necessary but we do it anyhow)
  for(int i=0;i<optsol->numvar;i++) { optsol->var[i] = xval[i]; }  
  grad(active);
  // save results 
  // contraints are multiplied by "-1" due to contraint definition  
  for (int i=0;i<optsol->numvar;i++) { df0dx[i] = optsol->gradobj[i]; }
  for(int j=0;j<optsol->numcon;j++) 
    {
      for(int i=0;i<optsol->numvar;i++) { dfdx[j][i] = -1.0 * optsol->gradcon[i][j]; }
    }
  return;
}

//------------------------------------------------------------------------------

void OptalgNlpmma::print() {


   fprintf(optunitout,"\n\tSolution Strateg #%d. : nlpmma\n",num);
   fprintf(optunitout,  "\t==============================\n\n");

   optgrad.print(optunitout);
}

//------------------------------------------------------------------------------

void OptalgNlpmma::printres() {


   fprintf(optunitout,"\n\tResult of %d. Strategy: nlpmma\n\n",num);
   
   switch (ifail) {
   
     case 0:
        fprintf(optunitout,
	"\tTHE ALGORITHM HAS CONVERGED.\n");
	break;
     case 1:
	fprintf(optunitout,
	"\tTHE ALGORITHM HAS BEEN STOPPED AFTER MAXIT ITERATIONS.\n");
	break;
     default:
	fprintf(optunitout,
        "\tError Message not specified.\n");
   }
}     
//------------------------------------------------------------------------------

void OptalgNlpsal::setDefault()
{
    acc    = 1.0e-6;
    rpen   = 1,0;
    penfac = 0.1;
    uinit  = 0.0e0;
    accloc = 1.0e-6;
    bscal  = 1.0e2;
    relbo  = 1.0;

    nbfgs  = 100;
    maxit  = 10;
    maxfun = 10;
    iprint = 1;
}

//------------------------------------------------------------------------------

void OptalgNlpsal::buildalg( nlpdata & param) 
{
    setDefault();

    if (param.rflag[0]) acc    = param.rval[0];
    if (param.rflag[1]) rpen   = param.rval[1];
    if (param.rflag[2]) penfac = param.rval[2];
    if (param.rflag[3]) uinit  = param.rval[3];
    if (param.rflag[4]) accloc = param.rval[4];
    if (param.rflag[5]) bscal  = param.rval[5];
    if (param.rflag[6]) relbo  = param.rval[6];

    if (param.iflag[0]) nbfgs  = param.ival[0];
    if (param.iflag[1]) maxit  = param.ival[1];
    if (param.iflag[2]) maxfun = param.ival[2];
    if (param.iflag[3]) iprint = param.ival[3];
} 

//------------------------------------------------------------------------------

void OptalgNlpsal::solve(Optsol *_optsol) {

    optsol =_optsol;    

    char * outfile = optsol->optprotfile;
    int    fsize   = optsol->fsize;

    int n = optsol->numvar;
    int m = optsol->numcon; 

    int lenwa = 2*nbfgs*n + 4*n + 11*nbfgs*nbfgs + 8*nbfgs;

    double rmem = sizeof(double)*(m+n+lenwa+29) + sizeof(int)*(4*n+m+44);

    filePrint(stderr," ... Allocating memory for SAL: %8.2f Mb.\n",rmem/1024.0/1024.0);    

    double *u     = new double[m];
    double *wa    = new double[lenwa];

    double *dfloc = new double[n];
    double *xupp  = new double[n];
    double *xlow  = new double[n];
    double *xsave = new double[n];

    double *dsave = new double[29];

    int *nbd      = new int[n];
    int *iwa      = new int[3*n];
    int *iactive  = new int[m];
    int *isave    = new int[44];

    _FORTRAN(nlpsal) (this,
		      optsol->var,
		      optsol->varlow, optsol->varup,
                      optsol->obj,    optsol->con,
                      optsol->gradobj,optsol->pgradcon,
                      u,              wa,                dfloc,
                      xupp,           xlow,              dsave,   xsave,         
                      nbd,            iwa,               iactive,        
                      isave,
                      optsol->numvar, optsol->numcon,    optsol->numeqc,
		      acc,            rpen,              penfac,
                      uinit,          accloc,            bscal,
                      relbo,          nbfgs,             lenwa,
                      maxit,          maxfun,            ifail,
                      iprint,         outfile,           fsize);
		      
    printres();

    // free memory used in algorithm

    delete [] u    ;
    delete [] wa   ;
    delete [] dfloc;
    delete [] dsave;
    delete [] xupp ;
    delete [] xlow ;
    delete [] xsave;
    delete [] nbd    ;
    delete [] iwa    ;
    delete [] iactive;
    delete [] isave  ;
}

//------------------------------------------------------------------------------

void OptalgNlpsal::func(int iter) {

   optsol->func(iter);
   
   ordercon();
}

//------------------------------------------------------------------------------

void OptalgNlpsal::grad(int *active) {

   reorderActive(active,optsol->actcon);   

   optsol->optpro->optvar->updvar();      

   optgrad.grad(optsol->actcon);

   optsol->optpro->optvar->resetvar();    
   
   ordergradcon();
}

//------------------------------------------------------------------------------

void OptalgNlpsal::print() {

   fprintf(optunitout,"\n\tSolution Strateg #%d. : nlpsal\n",num);
   fprintf(optunitout,  "\t==============================\n\n");

   fprintf(optunitout,"\tAccuracy .......................  : %12.5e\n",acc);
   fprintf(optunitout,"\tPenalty factor..................  : %12.5e\n",rpen);
   fprintf(optunitout,"\tPenalty increasing factor.......  : %12.5e\n",penfac);
   fprintf(optunitout,"\tInitial lagrange mulitplier.....  : %12.5e\n",uinit);
   fprintf(optunitout,"\tAccuracy in local optimizer.....  : %12.5e\n",accloc);
   fprintf(optunitout,"\tScaling in local optimizer......  : %12.5e\n",bscal);
   fprintf(optunitout,"\tStep size adaption factor.......  : %12.5e\n",relbo);

   fprintf(optunitout,"\tNumber of BFGS Correction Vectors : %12d\n"  ,maxfun);
   fprintf(optunitout,"\tMaximum Number of Function Calls  : %12d\n"  ,maxfun);
   fprintf(optunitout,"\tMaximum Number of Iterations ...  : %12d\n\n",maxit);
   
   fprintf(optunitout,"\tPrint Mode .....................  : %d\n",iprint);

   optgrad.print(optunitout);
}

//------------------------------------------------------------------------------

void OptalgNlpsal::printres() {

   fprintf(optunitout,"\n\tResult of %d. Strategy: nlpsal\n\n",num);
   
}     

#endif
