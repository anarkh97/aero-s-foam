#ifdef STRUCTOPT

#include "mma.h"

void MMAgc::raaupdate()
{
   int i,j;

   // MATLAB raacofmin = 1.0e-13

   double raacofmin = 1.0e-7;

   double raacof = 0.0;
    
   for (i=0;i<numvar;i++) {
     double yyux = (xmma[i]-xval[i])/(upp[i]-xmma[i]);
     double yyxl = (xmma[i]-xval[i])/(xmma[i]-low[i]);
     raacof += 0.5*yyux*yyxl;
   }

   raacof = max(raacof,raacofmin);

   if (f0valnew > (f0app+0.5*epsimin) ) {
     double deltaraa0 = (1.0/raacof)*(f0valnew-f0app);
     double zz0 = 1.1*(raa0 + deltaraa0);
     zz0 = min(zz0,10.0*raa0);
     raa0 = zz0;
   }

   for (j=0;j<numcon;j++) {
     double fappe = fapp[j] +0.5*epsimin;
     double fdelta = fvalnew[j]-fappe;
     if (fdelta > 0.0 ) {
       double deltaraa = (1.0/raacof)*(fvalnew[j]-fapp[j]);
       double zzz = 1.1*(raa[j] + deltaraa);
       raa[j] = min(zzz,10.0*raa[j]);
     }
   }
}

#endif
