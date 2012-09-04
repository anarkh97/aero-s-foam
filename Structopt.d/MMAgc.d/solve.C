#ifdef STRUCTOPT

#include <cmath>
#include <cstdio>

#include "mma.h"
#include <Structopt.d/Optsol.h>

int MMAgc::solve()
{

   int i,j;
   
   // set counter

   int iter  = 0;
   int istop = 0;

   int nfunc = 0;
   int nsens = 0;

   // evaluate functions and gradients

   optalg->func(iter,xval,f0val,fval);
   optalg->grad(xval,df0dx,dfdx,active);

   nfunc++;
   nsens++;

   // initialize raa0 and raa
   
   // MATLAB version
   // raa0 = 1.0;
   // for (j=0;j<numcon;j++)  raa[j] = 1.0;

   // paper version

   raa0 = 1.0e-4;
   for (j=0;j<numcon;j++) raa[j] = 1.0e-4;

   double sclraa = 0.1/numvar;

   for (i=0;i<numvar;i++) {
     raa0 += sclraa*abs(df0dx[i])/(xmax[i]-xmin[i]); 
     for (j=0;j<numcon;j++) raa[j] += sclraa*abs(dfdx[j][i])/(xmax[i]-xmin[i]); 
   }

   // initialize c_j based on initial solution
   
   // MATLAB:  c[j] = 1.Oe3;
   
   for (j=0;j<numcon;j++) c[j] = abs(f0val);

   // print initial setup

   printsol(iter);

   // Main loop

   while ( iter < maxiter && istop == 0) {

     iter++;

     // Calculate low, upp, raa0 and raa:

     asymp(iter);

     // Generate and solve the subproblem:
  
     gcmmasub(iter);
  
     // Check if conservative approximations:

     optalg->func(iter,xmma,f0valnew,fvalnew);

     nfunc = nfunc + 1;

     int conserv = 1;

     if ( (f0app+epsimin) < f0valnew)  conserv = 0;

     for(j=0;j<numcon;j++) 
       if ( (fapp[j] + epsimin) < fvalnew[j])  conserv = 0;
    
     int innerit=0;

     // While non-conservative approximations, solve again:

     if ( conserv == 0 ) {

       while ( conserv == 0 && innerit < maxinnerit ) {

         innerit = innerit+1;
  
         // Calculate new values on raa0 and raa:
  
         raaupdate();
        
         // Solve the subproblem with these new raa0 and raa:
   
         gcmmasub(iter);
  
         // Evaluate criteria
        
         optalg->func(iter,xmma,f0valnew,fvalnew);
  
         nfunc = nfunc + 1;
       
         // check if conservative approximations:

         conserv = 1;
  
         if ( (f0app+epsimin) < f0valnew)  conserv = 0;
  
         for(j=0;j<numcon;j++) {
           if ( (fapp[j] + epsimin) < fvalnew[j] )  conserv = 0;
         }
       }
     }
    
    fprintf(ofile,"\nNUMBER OF CYCLES: %d\n\n",innerit);

    // swap variables

    for (i=0;i<numvar;i++) {
      xold2[i] = xold1[i];
      xold1[i] = xval[i];
      xval[i]  = xmma[i];
    }
    
    // swap solution
    
    f0val = f0valnew;
    
    for (j=0;j<numcon;j++) fval[j] = fvalnew[j];
    
    // evaluate gradients

    optalg->grad(xval,df0dx,dfdx,active);
 
    nsens = nsens + 1;
 
    printsol(iter);

    // Calculate the residual vector of the KKT conditions:
 
    istop = kktcheck();

  }

  // print final solution and solver statistics

  fprintf(ofile,"\n\n");

  
  fprintf(ofile," FINAL SOLUTION OF GCMMA\n");
  fprintf(ofile," =======================\n\n");

  printsol(iter);

  fprintf(ofile,"NUMBER OF OUTER ITERATIONS     : %4d\n",iter);
  fprintf(ofile,"NUMBER OF FUNCTION EVALUATIONS : %4d\n",nfunc);
  fprintf(ofile,"NUMBER OF GRADIENT EVALUATIONS : %4d\n",nsens);
  
  return (iter < maxiter) ?  0 : 1;
}

#endif
