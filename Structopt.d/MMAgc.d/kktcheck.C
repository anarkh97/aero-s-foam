#ifdef STRUCTOPT

#include <cmath>
#include <cstdio>

#include "mma.h"

int MMAgc::kktcheck()
{
  int i,j;
  
  double residunorm = 0.0;
  double residumax  =-1.0;
  
  for (i=0;i<numvar;i++) {

    double rex   = df0dx[i] - xsi[i] + eta[i];
    double rexsi = xsi[i]*(xmma[i]-xmin[i]);
    double reeta = eta[i]*(xmax[i]-xmma[i]);

    for (j=0;j<numcon;j++) rex += dfdx[j][i]*lam[j];

    residunorm += rex*rex;
    residumax   = max(residumax,abs(rex));
   
    residunorm += rexsi*rexsi;
    residumax   = max(residumax,abs(rexsi));

    residunorm += reeta*reeta;
    residumax   = max(residumax,abs(reeta));
  }
  
  double rez   = a0 - zet;
  double rezet = zet*z;
  
  for (j=0;j<numcon;j++) {
    double rey   = c[j] + d[j]*y[j]  - mu[j]  - lam[j] ;
    double relam = fval[j] - z*a[j] - y[j] + s[j];
    double remu  = mu[j]*y[j];
    double res   = lam[j]*s[j];
           rez  -= a[j]*lam[j];

    residunorm += rey*rey;
    residumax   = max(residumax,abs(rey));

    residunorm += relam*relam;
    residumax   = max(residumax,abs(relam));

    residunorm += remu*remu;
    residumax   = max(residumax,abs(remu));

    residunorm += res*res;
    residumax   = max(residumax,abs(res));
  }

  residunorm += rez*rez;
  residumax   = max(residumax,abs(rez));

  residunorm  = sqrt(residunorm);

  fprintf(ofile,"RESIDUAL NORM OF K.K.T.: %20.10e\n",residunorm);
  fprintf(ofile,"MAXIMUM RESIDUAL COMP. : %20.10e\n",residumax);

  return ( max(residunorm,residumax) > acc ) ?  0 : 1;
}

#endif
