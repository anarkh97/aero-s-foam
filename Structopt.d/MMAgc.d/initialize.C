#ifdef STRUCTOPT

#include <cstdio>
#include <cmath>

#include "mma.h"


void MMAgc::initialize(double* _xval,double* _xmax,double* _xmin)
{
   int i,j;

   epsimin = sqrt(static_cast<double>(numcon+numvar))*1.0e-9;
    
   // MATLAB  vers 1: raa0eps  = 1.0e-5;
   // MATLAB  vers 2: raa0eps  = 1.0e-6;
     
   raa0eps  = 1.0e-7;
   a0       = 1.0;
   f0app    = 0.0;

   // initialize initial solution
   
   xval  = new double[numvar];
   xmax  = new double[numvar];
   xmin  = new double[numvar];
   xold1 = new double[numvar];
   xold2 = new double[numvar];
   low   = new double[numvar];
   upp   = new double[numvar];
   alfa  = new double[numvar];
   beta  = new double[numvar];

   
   for (i=0;i<numvar;i++) {
     xval[i]  = _xval[i];
     xmax[i]  = _xmax[i];
     xmin[i]  = _xmin[i];
     xold1[i] = xval[i];
     xold2[i] = xval[i];
     low[i]   = xmin[i];
     upp[i]   = xmax[i];
     alfa[i]  = xmin[i];
     beta[i]  = xmax[i];
   }
   
   fval    = new double[numcon];
   fvalnew = new double[numcon];
   df0dx   = new double[numvar];
   dfdat   = new double[numvar*numcon];
   dfdx    = new double*[numcon];

   for(j=0;j<numcon;j++) dfdx[j] = dfdat + j*numvar;

   a = new double[numcon];
   b = new double[numcon];
   c = new double[numcon];
   d = new double[numcon];
   r = new double[numcon];
  
   raa    = new double[numcon];

   for (j=0;j<numcon;j++) {
     a[j]   = 0.0;
     d[j]   = 0.0;
   }

   p0  = new double[numvar];
   q0  = new double[numvar];

   Pd = new double[numvar*numcon];
   P  = new double*[numcon];

   for(i=0;i<numcon;i++) P[i] = Pd + i*numvar;

   Qd = new double[numvar*numcon];
   Q  = new double*[numcon];
   
   for(i=0;i<numcon;i++) Q[i] = Qd + i*numvar;

   GGd = new double[numvar*numcon];
   GG  = new double*[numcon];
   
   for(i=0;i<numcon;i++) GG[i] = GGd + i*numvar;
   
   if ( numcon < numvar ) {
     int numcon1 = numcon + 1;
     AAd = new double[numcon1*numcon1];
     AA  = new double*[numcon1];
     for(i=0;i<numcon1;i++) AA[i] = AAd + i*numcon1;
     bb  = new double[numcon1];
   }
   else {
     int numvar1 = numvar + 1;
     AAd = new double[numvar1*numvar1];
     AA  = new double*[numvar1];
     for(i=0;i<numvar1;i++) AA[i] = AAd + i*numvar1;
     bb  = new double[numvar1];
   }
   // arrays for subsolve
   
   x        = new double[numvar];
   xold     = new double[numvar];
   dx       = new double[numvar+1];
   delx     = new double[numvar];
   xmma     = new double[numvar];
   diagx    = new double[numvar];
   diagxinv = new double[numvar];
   xsi      = new double[numvar];  
   xsiold   = new double[numvar];  
   dxsi     = new double[numvar];  
   xsimma   = new double[numvar]; 
   eta      = new double[numvar];
   etaold   = new double[numvar];
   deta     = new double[numvar];
   etamma   = new double[numvar];

   y            = new double[numcon];
   yold         = new double[numcon];
   dy           = new double[numcon];
   dely         = new double[numcon];
   diagy        = new double[numcon];
   ymma         = new double[numcon];
   diagyinv     = new double[numcon];
   lam          = new double[numcon]; 
   lamold       = new double[numcon]; 
   dlam         = new double[numcon+1]; 
   dellam       = new double[numcon];
   dellamyi     = new double[numcon];
   lamma        = new double[numcon];
   diaglam      = new double[numcon];
   diaglamyi    = new double[numcon];
   diaglamyiinv = new double[numcon];
   mu           = new double[numcon];
   muold        = new double[numcon];
   dmu          = new double[numcon];
   mumma        = new double[numcon];
   s            = new double[numcon]; 
   sold         = new double[numcon]; 
   ds           = new double[numcon]; 
   smma         = new double[numcon]; 
   gvec         = new double[numcon];
 
   fapp    	= new double[numcon]; 
  
   // pseudo active strategy set

   active = new int[numcon];

   for (i=0;i<numcon;i++) active[i]=1;
   
   // compute total memory requirement
   
   double memsize = 27 * sizeof(double) * numvar
                  + 35 * sizeof(double) * numcon
                  +  1 * sizeof(int)    * numcon
	          +  4 * sizeof(double) * numcon * numvar
	          +  3 * sizeof(double*) * numcon;
	       
   if ( numcon < numvar )
      memsize +=  sizeof(double)  * (numcon+1) * (numcon+2)
               +  sizeof(double*) * (numcon+1);
   else
      memsize +=  sizeof(double)  * (numvar+1) * (numvar+2)
               +  sizeof(double*) * (numvar+1);
	       
   fprintf(stderr," ... Allocating memory for GCM: %f MB\n",memsize/1024.0/1024.0);
         
}

#endif
