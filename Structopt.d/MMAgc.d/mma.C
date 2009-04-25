#ifdef STRUCTOPT

#include <stdio.h>

#include "mma.h"


MMAgc::MMAgc(OptalgNlpmma* _optalg,
             double* _xval, double* _xmax, double* _xmin,
             int _numvar, int _numcon,
	     int _maxiter, int _maxinnerit,
	     double& _acc, double& _saFac, double& _sbFac, 
	     double& _scFac, double& _dstep, FILE* _ofile)
{	     
   int i,j;

   // parameter 

   saFac = _saFac;
   sbFac = _sbFac;
   scFac = _scFac;
   
   dstep = _dstep;

   optalg = _optalg;

   numvar = _numvar;
   numcon = _numcon;

   maxiter    = _maxiter;
   maxinnerit = _maxinnerit;
   
   acc = _acc;
   
   // output file
   
   ofile =_ofile;

   // allocate memory and initialize 

   initialize(_xval,_xmax,_xmin);
   
   // print parameters
   
   fprintf(ofile,"\n\n");
   
   fprintf(ofile,"Globally Convergent MMA\n");
   fprintf(ofile,"=======================\n\n");

   fprintf(ofile,"maximum number of iterations  : %d\n",maxiter);   
   fprintf(ofile,"maximum number of sub cycles  : %d\n",maxinnerit);
   fprintf(ofile,"accuracy                      : %e\n",acc);
   fprintf(ofile,"asymptotes adaptation factor  : %e\n",saFac);
   fprintf(ofile,"asymptotes adaptation factor  : %e\n",sbFac);
   fprintf(ofile,"asymptotes adaptation factor  : %e\n",scFac);
   fprintf(ofile,"step size in local problem    : %e\n",dstep);

   fprintf(ofile,"\n\n");

   fflush(ofile);
}

#endif
