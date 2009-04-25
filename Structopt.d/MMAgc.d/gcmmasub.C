#ifdef STRUCTOPT

#include "mma.h"

void MMAgc::gcmmasub(int iter)
{

  // Calculation of the bounds alfa and beta;

  int i,j;
  
  for(i=0;i<numvar;i++) {
   alfa[i] = max(xmin[i], low[i] + 0.1*(xval[i]-low[i]));
   alfa[i] = max(alfa[i], xval[i] - dstep*(xmax[i]-xmin[i])); 
   
   beta[i] = min(xmax[i], upp[i] - 0.1*(upp[i]-xval[i]));
   beta[i] = min(beta[i], xval[i] + dstep*(xmax[i]-xmin[i]));
  }
  
  // Calculations of p0, q0, r0, P, Q, r and b.

  r0  = f0val;

  for(j=0;j<numcon;j++) r[j] = fval[j];
  
  for(i=0;i<numvar;i++) {

    double ux1 = upp[i]-xval[i];
    double ux2 = ux1*ux1;
    double xl1 = xval[i]-low[i];
    double xl2 = xl1*xl1;
    double ul1 = upp[i]-low[i];

    double uxinv = 1.0/ux1;
    double xlinv = 1.0/xl1;
    double ulinv = 1.0/ul1;

    p0[i] = 0.0;
    if (df0dx[i] > 0) p0[i] = df0dx[i];
    p0[i] += 0.5*raa0*ulinv;
    p0[i] *= ux2;

    q0[i] = 0.0;
    if (df0dx[i] < 0) q0[i] = -df0dx[i];
    q0[i] += 0.5*raa0*ulinv;
    q0[i] *= xl2;

    r0    -=  p0[i]*uxinv + q0[i]*xlinv;
  
    for (j=0;j<numcon;j++) {

      P[j][i] = 0.0;
      if (dfdx[j][i] > 0) P[j][i] = dfdx[j][i];
      P[j][i] += 0.5*ulinv*raa[j];
      P[j][i] *= ux2;

      Q[j][i] = 0.0;
      if (dfdx[j][i] < 0) Q[j][i] = -dfdx[j][i];
      Q[j][i] += 0.5*ulinv*raa[j];
      Q[j][i] *= xl2;

      r[j]    -= uxinv*P[j][i] + xlinv*Q[j][i];
    }
  }

  for(j=0;j<numcon;j++) b[j] = -r[j];
 
  // Solving the subproblem by a primal-dual Newton method

  subsolve();

  // Calculations of f0app and fapp.

  f0app = r0;

  for(j=0;j<numcon;j++) fapp[j] = r[j];

  for(i=0;i<numvar;i++) {

    double ux1 = upp[i]-xmma[i];
    double xl1 = xmma[i]-low[i];

    double uxinv = 1.0/ux1;
    double xlinv = 1.0/xl1;

    f0app += uxinv*p0[i] + xlinv*q0[i];

    for(j=0;j<numcon;j++) fapp[j]  += uxinv*P[j][i] + xlinv*Q[j][i];
  }
}  

#endif
