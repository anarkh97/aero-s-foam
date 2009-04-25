#ifdef STRUCTOPT

#include "mma.h"

#include <math.h>

void MMAgc::subsolve()
{
  //
  // This function subsolv solves the MMA subproblem:
  //         
  // minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
  //          + SUM[ ci*yi + 0.5*di*(yi)^2 ],
  //
  // subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
  //            alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
  //        
  // Input:  m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
  // Output: xmma,ymma,zmma, slack variables and Lagrange multiplers.


  int i,j,ii,jj;

  // initialize local arrays


  for (i=0;i<numvar;i++) {
      x[i] = 0.5*(alfa[i]+beta[i]);
    xsi[i] = max(1.0, 1.0/(x[i]-alfa[i]));
    eta[i] = max(1.0, 1.0/(beta[i]-x[i]));
  }
  
  for (j=0;j<numcon;j++) {
      y[j] = 1.0;
    lam[j] = 1.0;
     mu[j] = max(1.0,0.5*c[j]);
     s[j]  = 1.0;
  }

  z    = 1.0;
  zet  = 1.0;

  double epsi = 1.0;

  int itera   = 0;

  // -------------------------------------------------------------
      
  while (epsi > epsimin) {

    double residunorm = 0.0;
    double residumax  =-1.0;
 
    for (j=0;j<numcon;j++) gvec[j] = 0.0;
  
    // -------------------------------------------------------------
     
    for (i=0;i<numvar;i++) {
  
      double ux1    = upp[i]-x[i];
      double xl1    = x[i]-low[i];
      double ux2    = ux1*ux1;
      double xl2    = xl1*xl1;
      double uxinv1 = 1.0/ux1;
      double xlinv1 = 1.0/xl1;

      double plam = p0[i]; 
      for (j=0;j<numcon;j++) plam += P[j][i]*lam[j];

      double qlam = q0[i];
      for (j=0;j<numcon;j++) qlam += Q[j][i]*lam[j];

      double dpsidx = plam/ux2 - qlam/xl2 ;

      double rex    = dpsidx - xsi[i] + eta[i];
      double rexsi  = xsi[i]*(x[i]-alfa[i]) - epsi;
      double reeta  = eta[i]*(beta[i]-x[i]) - epsi;

      for (j=0;j<numcon;j++) gvec[j] += P[j][i]*uxinv1 + Q[j][i]*xlinv1;
   
      residunorm += rex*rex; 
      residumax   = max(residumax,abs(rex));

      residunorm += rexsi*rexsi; 
      residumax   = max(residumax,abs(rexsi));
   
      residunorm += reeta*reeta; 
      residumax   = max(residumax,abs(reeta));
   
    }

    // -------------------------------------------------------------

    double rez   = a0 - zet; 
    double rezet = zet*z - epsi;

    for (j=0;j<numcon;j++) {
      double rey   = c[j] + d[j]*y[j] - mu[j] - lam[j];
      double relam = gvec[j] - z*a[j] - y[j] + s[j] - b[j];
      double remu  = mu[j]*y[j] - epsi;
      double res   = lam[j]*s[j] - epsi;
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
    residumax	= max(residumax,abs(rez));

    residunorm += rezet*rezet; 
    residumax	= max(residumax,abs(rezet));

    residunorm  = sqrt(residunorm);

    // -------------------------------------------------------------
     
    int ittt = 0;

    while (residumax > 0.9*epsi && ittt < 200) {

      ittt++;
      itera++;

      for (j=0;j<numcon;j++) gvec[j] = 0.0;

      for (i=0;i<numvar;i++) {
  
  	double ux1 = upp[i]-x[i];
  	double xl1 = x[i]-low[i];
  	double ux2 = ux1*ux1;
  	double xl2 = xl1*xl1;
  	double ux3 = ux1*ux2;
  	double xl3 = xl1*xl2;
     
  	double uxinv1 = 1.0/ux1;
  	double xlinv1 = 1.0/xl1;
  	double uxinv2 = 1.0/ux2;
  	double xlinv2 = 1.0/xl2;
  
  	double plam = p0[i]; 
  	for (j=0;j<numcon;j++) plam += P[j][i]*lam[j];

  	double qlam = q0[i];
  	for (j=0;j<numcon;j++) qlam += Q[j][i]*lam[j];

  	double dpsidx = plam/ux2 - qlam/xl2 ;
  	delx[i] = dpsidx - epsi/(x[i]-alfa[i]) + epsi/(beta[i]-x[i]);
  
  	diagx[i]    = plam/ux3 + qlam/xl3;
  	diagx[i]    = 2.0*diagx[i] + xsi[i]/(x[i]-alfa[i]) + eta[i]/(beta[i]-x[i]);
  	diagxinv[i] = 1.0/diagx[i];

  	for (j=0;j<numcon;j++) {
  	  gvec[j]  += P[j][i]*uxinv1 + Q[j][i]*xlinv1;
  	  GG[j][i]  = P[j][i]*uxinv2 - Q[j][i]*xlinv2;
  	}
      }
 
      delz = a0 - epsi/z;
 
      for (j=0;j<numcon;j++) {
  	     dely[j] = c[j] + d[j]*y[j] - lam[j] - epsi/y[j];
  	   dellam[j] = gvec[j] - z*a[j] - y[j] - b[j] + epsi/lam[j];
  	    diagy[j] = d[j] + mu[j]/y[j];
  	 diagyinv[j] = 1.0/diagy[j];
  	  diaglam[j] = s[j]/lam[j];
  	diaglamyi[j] = diaglam[j]+diagyinv[j];
  	     delz   -= a[j]*lam[j]; 
      }

      // -------------------------------------------------------------
   
      if ( numcon < numvar ) {

  	AA[numcon][numcon] = -zet/z;

  	for (j=0;j<numcon;j++) {
  	  bb[j] =  dellam[j] + dely[j]/diagy[j];
  	  for (i=0;i<numvar;i++) bb[j] -= GG[j][i]*(delx[i]/diagx[i]);
  	}
    
  	bb[numcon] = delz;
     
  	for (j=0;j<numcon;j++) {
  	  for (jj=0;jj<numcon;jj++) {
  	    AA[j][jj] = 0.0;
  	    for (i=0;i<numvar;i++) AA[j][jj]+= GG[j][i]*GG[jj][i]*diagxinv[i];
  	  }
  	  AA[j][j] += diaglamyi[j];
  	  AA[numcon][j] = a[j];
          AA[j][numcon] = a[j];
  	}

  	linsolve(AA,dlam,bb,numcon+1);
  	dz = dlam[numcon];
       
  	for (i=0;i<numvar;i++) {
  	  dx[i] = -delx[i]/diagx[i];
          for (j=0;j<numcon;j++) dx[i] -= (GG[j][i]*dlam[j])/diagx[i];
  	}
      }
      else {

  	AA[numvar][numvar] = zet/z;

  	for (j=0;j<numcon;j++) {
  	  diaglamyiinv[j]     = 1.0/diaglamyi[j];
  	  dellamyi[j]	      = dellam[j] + dely[j]/diagy[j];
          AA[numvar][numvar] += a[j]*a[j]/diaglamyi[j];
  	}
     
  	for (i=0;i<numvar;i++) {       
          for (ii=0;ii<numvar;ii++) {
            AA[i][ii] = 0.0;
            for (j=0;j<numcon;j++) AA[i][ii] += GG[j][i]*GG[j][ii]*diaglamyiinv[j];
          }   
  	  AA[i][i] += diagx[i];
          AA[numvar][i] = 0.0;
          AA[i][numvar] = 0.0;
  	  for (j=0;j<numcon;j++) AA[numvar][i] -= GG[j][i]*a[j]/diaglamyi[j];
  	  for (j=0;j<numcon;j++) AA[i][numvar] -= GG[j][i]*a[j]/diaglamyi[j];
          bb[i]  = -delx[i]; 
  	  for (j=0;j<numcon;j++) bb[i] -= GG[j][i]*dellamyi[j]/diaglamyi[j];
  	}
  	bb[numvar]  = -delz;
  	for (j=0;j<numcon;j++) bb[numvar] += a[j]*dellamyi[j]/diaglamyi[j];

  	linsolve(AA,dx,bb,numvar+1);
  	dz = dx[numvar];

  	for (j=0;j<numcon;j++) {
  	  dlam[j] =  - dz*a[j]/diaglamyi[j] + dellamyi[j]/diaglamyi[j];
          for (i=0;i<numvar;i++) dlam[j] += GG[j][i]*dx[i]/diaglamyi[j];
  	}
      }
    
      // -------------------------------------------------------------

      double dzet    = -zet + epsi/z - zet*dz/z;
   
      double stmxx   = -1.0e99;
      double stmalfa = -1.0e99;
      double stmbeta = -1.0e99; 

      stmxx = max(stmxx,-1.01*dz/z);
      stmxx = max(stmxx,-1.01*dzet/zet);
  
      for (i=0;i<numvar;i++) {
  	dxsi[i] = -xsi[i] + epsi/(x[i]-alfa[i]) - (xsi[i]*dx[i])/(x[i]-alfa[i]);
  	deta[i] = -eta[i] + epsi/(beta[i]-x[i]) + (eta[i]*dx[i])/(beta[i]-x[i]);

  	stmxx	= max(stmxx,-1.01*dxsi[i]/xsi[i]);
  	stmxx	= max(stmxx,-1.01*deta[i]/eta[i]);

  	stmalfa = max(stmalfa,-1.01*dx[i]/(x[i]-alfa[i]));
  	stmbeta = max(stmbeta, 1.01*dx[i]/(beta[i]-x[i]));
      }
   
      for (j=0;j<numcon;j++) {
  	 dy[j] = -dely[j]/diagy[j] + dlam[j]/diagy[j];
  	dmu[j] = -mu[j] + epsi/y[j] - (mu[j]*dy[j])/y[j];
  	 ds[j] = -s[j] + epsi/lam[j] - (s[j]*dlam[j])/lam[j] ;

  	stmxx  = max(stmxx,-1.01*  dy[j]/  y[j]);
  	stmxx  = max(stmxx,-1.01*dlam[j]/lam[j]);
  	stmxx  = max(stmxx,-1.01* dmu[j]/ mu[j]);
  	stmxx  = max(stmxx,-1.01*  ds[j]/  s[j]);
      }
   
      double stmalbe   = max(stmalfa,stmbeta);
      double stmalbexx = max(stmalbe,stmxx);
      double stminv    = max(stmalbexx,1.0);
      double steg      = 1.0/stminv;

      // -------------------------------------------------------------

      zold   =   z;
      zetold =  zet;

      for (i=0;i<numvar;i++) {
  	  xold[i] =   x[i];
  	xsiold[i] = xsi[i];
  	etaold[i] = eta[i];
      }
     
      for (j=0;j<numcon;j++) {
  	  yold[j] =   y[j];
  	lamold[j] = lam[j];
  	 muold[j] =  mu[j];
  	  sold[j] =   s[j];
      }

      // -------------------------------------------------------------

      int itto = 0;
   
      double resinew = 2.0*residunorm;
      double resimax = residumax;

      while (resinew > residunorm && itto < 50) {

  	itto++;
     
  	resinew =  0.0;
  	resimax = -1.0;

  	z   =	zold + steg*dz;
  	zet = zetold + steg*dzet;

  	for (i=0;i<numvar;i++) {
  	    x[i] =   xold[i] +   steg*dx[i];
  	  xsi[i] = xsiold[i] + steg*dxsi[i];
  	  eta[i] = etaold[i] + steg*deta[i];
  	}
 
  	for (j=0;j<numcon;j++) {
  	    y[j] =   yold[j] +   steg*dy[j];
  	  lam[j] = lamold[j] + steg*dlam[j];
  	   mu[j] =  muold[j] +  steg*dmu[j];
  	    s[j] =   sold[j] +   steg*ds[j];
      
  	 gvec[j] = 0.0;
  	}

  	for (i=0;i<numvar;i++) {

  	  double ux1 = upp[i]-x[i];
  	  double xl1 = x[i]-low[i];
  	  double ux2 = ux1*ux1;
  	  double xl2 = xl1*xl1;
  	  double uxinv1 = 1.0/ux1;
  	  double xlinv1 = 1.0/xl1;

  	  double plam = p0[i]; 
  	  for (j=0;j<numcon;j++) plam += P[j][i]*lam[j];

  	  double qlam = q0[i];
  	  for (j=0;j<numcon;j++) qlam += Q[j][i]*lam[j];

  	  double dpsidx = plam/ux2 - qlam/xl2 ;

  	  double rex   = dpsidx - xsi[i] + eta[i];
  	  double rexsi = xsi[i]*(x[i]-alfa[i]) - epsi;
  	  double reeta = eta[i]*(beta[i]-x[i]) - epsi;

  	  for (j=0;j<numcon;j++) gvec[j] += P[j][i]*uxinv1 + Q[j][i]*xlinv1;

          resinew += rex*rex; 
          resimax  = max(resimax,abs(rex));

  	  resinew += rexsi*rexsi; 
          resimax  = max(resimax,abs(rexsi));

  	  resinew += reeta*reeta; 
          resimax  = max(resimax,abs(reeta));
  	}
     
  	rez   = a0 - zet; 
  	rezet = zet*z - epsi;
    
  	for (j=0;j<numcon;j++) {

  	  double rey   = c[j] + d[j]*y[j] - mu[j] - lam[j];
  	  double relam = gvec[j] - z*a[j] - y[j] + s[j] - b[j];
  	  double remu  = mu[j]*y[j] - epsi;
  	  double res   = lam[j]*s[j] - epsi;
  		 rez  -= a[j]*lam[j];

  	  resinew += rey*rey;
          resimax  = max(resimax,abs(rey));

  	  resinew += relam*relam; 
  	  resimax  = max(resimax,abs(relam));

  	  resinew += remu*remu; 
  	  resimax  = max(resimax,abs(remu));

  	  resinew += res*res;
          resimax  = max(resimax,abs(res));
  	} 
     
  	resinew += rez*rez; 
  	resimax  = max(resimax,abs(rez));

  	resinew += rezet*rezet; 
  	resimax  = max(resimax,abs(rezet));

  	resinew  = sqrt(resinew);
  	steg     = steg/2.0;
      }
      
      residunorm = resinew;
      residumax  = resimax;
      steg      *= 2.0;
    }
    epsi *= 0.1;
  }

  zmma   =   z;
  zetmma =  zet;


  for (i=0;i<numvar;i++){
      xmma[i] =   x[i];
    xsimma[i] = xsi[i];
    etamma[i] = eta[i];
  }

  for (j=0;j<numcon;j++){
     ymma[j] =   y[j];
    lamma[j] = lam[j];
    mumma[j] =  mu[j];
     smma[j] =   s[j];
  }
}    

#endif
