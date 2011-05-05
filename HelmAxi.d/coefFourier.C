#include <cstdio>
#include <cmath>
#include <HelmAxi.d/coefFourier.h>


// UH - 9/19/00
//
// We are using the long double versions of the bessel function
// Indeed for large n, at x fixed, jn(x) goes to 0 but the double
// function gives not a number.
// According to IRIX manual, long double operations are not portable.
// Long double operations on this system are only supported in round to
// nearest rounding mode (the default).


DComplex coefExpDir(int mode, double k, double rho, double phi) {

  DComplex Ex(0.0,0.0);
  double jbessel;

  if (fabs(rho)<1e-15) {
     if (mode==0) {
        Ex=DComplex(1.0,0.0); 
     }
  }
  else {
     if (mode==0) {
        Ex = DComplex(besselj(0, k*rho),0.0); 
     }
     else {
        if (mode%2==0) {
          int n = mode/2;
          phi *= n;
          jbessel = besselj(n, k*rho);
          switch (n%4) {
            case 0:
               Ex = 2*sin(phi)*jbessel*DComplex(1.0,0.0);
               break;
            case 1:
               Ex = 2*sin(phi)*jbessel*DComplex(0.0,1.0);
               break;
            case 2:
               Ex = 2*sin(phi)*jbessel*DComplex(-1.0,0.0);
               break;
            case 3:
               Ex = 2*sin(phi)*jbessel*DComplex(0.0,-1.0);
               break; 
          } 
        }
        else { 
          int n = (mode+1)/2;
          phi *= n;
          jbessel = besselj(n, k*rho);
          switch (n%4) {
            case 0:
               Ex = 2*cos(phi)*jbessel*DComplex(1.0,0.0);
               break;
            case 1:
               Ex = 2*cos(phi)*jbessel*DComplex(0.0,1.0);
               break;
            case 2:
               Ex = 2*cos(phi)*jbessel*DComplex(-1.0,0.0);
               break;
            case 3: 
               Ex = 2*cos(phi)*jbessel*DComplex(0.0,-1.0);
               break; 
          }
        }
     }
  }

  return Ex; 

} 


DComplex coefExpNeu(int mode, double k, double rho, double phi, double nr,
                    double nz, double dx, double dy, double dz) {

  DComplex result(0.0,0.0);
  DComplex E1(0.0,0.0);
  DComplex E2(0.0,0.0);
  DComplex E3(0.0,0.0);

  if (fabs(rho)<1e-15) {
     if (mode==0) {
        E1=DComplex(1.0,0.0);
     }
  }
  else {
     if (mode==0) {
        double kRho = k*rho;
        double jOne = besselj(1,kRho);
        E1 = DComplex(besselj(0, kRho),0.0);
        E2 = cos(phi)*jOne*DComplex(0.0,1.0); 
        E3 = sin(phi)*jOne*DComplex(0.0,1.0);
     }
     else {
        if (mode%2==0) {
          int n = mode/2;
          double nPhi = n*phi;
          double kRho = k*rho;
          double jPlus = besselj(n+1,kRho);
          double jZero = besselj(n,kRho);
          double jMinus = besselj(n-1,kRho);
          switch (n%4) {
            case 0:
               E1 = 2*sin(nPhi)*jZero*DComplex(1.0,0.0);
               E2 = jPlus*sin(nPhi+phi)*DComplex(0.0,1.0);
               E2 -= jMinus*sin(nPhi-phi)*DComplex(0.0,1.0);
               E3 = jPlus*cos(nPhi+phi)*DComplex(0.0,-1.0);
               E3 += jMinus*cos(nPhi-phi)*DComplex(0.0,-1.0);
               break;
            case 1:
               E1 = 2*sin(nPhi)*jZero*DComplex(0.0,1.0);
               E2 = jPlus*sin(nPhi+phi)*DComplex(-1.0,0.0);
               E2 -= jMinus*sin(nPhi-phi)*DComplex(-1.0,0.0);
               E3 = jPlus*cos(nPhi+phi)*DComplex(1.0,0.0);
               E3 += jMinus*cos(nPhi-phi)*DComplex(1.0,0.0);
               break;
            case 2:
               E1 = 2*sin(nPhi)*jZero*DComplex(-1.0,0.0);
               E2 = jPlus*sin(nPhi+phi)*DComplex(0.0,-1.0);
               E2 -= jMinus*sin(nPhi-phi)*DComplex(0.0,-1.0);
               E3 = jPlus*cos(nPhi+phi)*DComplex(0.0,1.0);
               E3 += jMinus*cos(nPhi-phi)*DComplex(0.0,1.0);
               break;
            case 3:
               E1 = 2*sin(nPhi)*jZero*DComplex(0.0,-1.0);
               E2 = jPlus*sin(nPhi+phi)*DComplex(1.0,0.0);
               E2 -= jMinus*sin(nPhi-phi)*DComplex(1.0,0.0);
               E3 = jPlus*cos(nPhi+phi)*DComplex(-1.0,0.0);
               E3 += jMinus*cos(nPhi-phi)*DComplex(-1.0,0.0);
               break;
          }
        }
        else {
          int n = (mode+1)/2;
          double nPhi = n*phi;
          double kRho = k*rho;
          double jPlus = besselj(n+1,kRho);
          double jZero = besselj(n,kRho);
          double jMinus = besselj(n-1,kRho);
          switch (n%4) {
            case 0:
               E1 = 2*cos(nPhi)*jZero*DComplex(1.0,0.0);
               E2 = jPlus*cos(nPhi+phi)*DComplex(0.0,1.0);
               E2 -= jMinus*cos(nPhi-phi)*DComplex(0.0,1.0);
               E3 = jPlus*sin(nPhi+phi)*DComplex(0.0,1.0);
               E3 += jMinus*sin(nPhi-phi)*DComplex(0.0,1.0);
               break;
            case 1:
               E1 = 2*cos(nPhi)*jZero*DComplex(0.0,1.0);
               E2 = jPlus*cos(nPhi+phi)*DComplex(-1.0,0.0);
               E2 -= jMinus*cos(nPhi-phi)*DComplex(-1.0,0.0);
               E3 = jPlus*sin(nPhi+phi)*DComplex(-1.0,0.0);
               E3 += jMinus*sin(nPhi-phi)*DComplex(-1.0,0.0);
               break;
            case 2:
               E1 = 2*cos(nPhi)*jZero*DComplex(-1.0,0.0);
               E2 = jPlus*cos(nPhi+phi)*DComplex(0.0,-1.0);
               E2 -= jMinus*cos(nPhi-phi)*DComplex(0.0,-1.0);
               E3 = jPlus*sin(nPhi+phi)*DComplex(0.0,-1.0);
               E3 += jMinus*sin(nPhi-phi)*DComplex(0.0,-1.0);
               break;
            case 3:
               E1 = 2*cos(nPhi)*jZero*DComplex(0.0,-1.0);
               E2 = jPlus*cos(nPhi+phi)*DComplex(1.0,0.0);
               E2 -= jMinus*cos(nPhi-phi)*DComplex(1.0,0.0);
               E3 = jPlus*sin(nPhi+phi)*DComplex(1.0,0.0);
               E3 += jMinus*sin(nPhi-phi)*DComplex(1.0,0.0);
               break;
          }
        }
     }
  }

  result = nz*dz*E1 + nr*dx*E2 + nr*dy*E3;
  result *= DComplex(0.0,k);

  return result;

}


double besselj(int n, double x) {

 double result = jn(n, x);

#if defined(sgi) && ! defined(_OPENMP)
 if (isnan(result)) {
   result = (double) jnl(n, (long double) x);
 }
#endif

 return result;

}
