#include <stdlib.h>
#include <Utils.d/dbg_alloca.h>
#include <math.h>
#include <HelmAxi.d/coefFourier.h>
#include <HelmAxi.d/Line2AxiSommer.h>
#include <Utils.d/MyComplex.h>
#include <Math.d/FullSquareMatrix.h>

// Sommerfeld b.c. contribution for 2-D elements
// Note that the boundary condition will be always 
// on an element edge. 

// UH - 9 august 2000 
//
// The node 2 is the middle point of the edge [0-1]


Line2AxiSommer::Line2AxiSommer(int n1, int n2, int n3) {

 nn[0] = n1;
 nn[1] = n2;
 nn[2] = n3;

}


void
Line2AxiSommer::setType(int t) {

 type = t;

}


void
Line2AxiSommer::setSurf(double aR, double aZ) {

 surfR0 = aR;
 surfZ0 = aZ;

}


void
Line2AxiSommer::renum(int *table) {

  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];

}


FullSquareMatrix
Line2AxiSommer::sommerMatrix(CoordSet &cs) {

 return sommerMatrix(cs,new double [3*3]);

}


FullSquareMatrix
Line2AxiSommer::sommerMatrix(CoordSet &cs, double *d) {

 FullSquareMatrix sommerM(3,d);
 return sommerM;

}


FullSquareMatrixC
Line2AxiSommer::turkelMatrix(CoordSet &cs, double kappa, int mode) {

 return turkelMatrix(cs,kappa,mode,new DComplex [3*3]);

}


FullSquareMatrixC
Line2AxiSommer::turkelMatrix(CoordSet &cs, double kappa, int mode, 
                DComplex *d) {

 Node nd1 = cs.getNode(nn[0]);
 Node nd2 = cs.getNode(nn[1]);
 Node nd3 = cs.getNode(nn[2]);

 double x[3], y[3];
 double t0, t1;
 double length;

 x[0] = nd1.x; y[0] = nd1.y;
 x[1] = nd2.x; y[1] = nd2.y;
 x[2] = nd3.x; y[2] = nd3.y;

 FullSquareMatrixC turkelM(3,d);
 turkelM.zero();

 int i;

 switch (type) {
   default:
     fprintf(stderr,"BC parameter not available -> Sommerfeld\n");
   case 0:

     // Sommerfeld condition

     length = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));

     for (i=0; i<3;++i) {
       double xi;
       double weight;

       switch (i) {
         case 0:
           xi = -sqrt(0.6);
           weight = 5.0/9.0;
           break;
         case 1:
           xi = 0.0;
           weight = 8.0/9.0;
           break;
         case 2: 
           xi = sqrt(0.6);
           weight = 5.0/9.0;
           break;
       }

       weight *= 0.5*length;

       double r = 0.5*(x[0]+x[1])+0.5*xi*(x[1]-x[0]);

       double n0 = 0.5*(xi-1)*xi;
       double n1 = 0.5*(1+xi)*xi;
       double n2 = 1-xi*xi;

       turkelM[0][0] += weight*n0*n0*r*DComplex(0.0,kappa);
       turkelM[0][1] += weight*n0*n1*r*DComplex(0.0,kappa);
       turkelM[0][2] += weight*n0*n2*r*DComplex(0.0,kappa);
       turkelM[1][0] = turkelM[0][1];
       turkelM[1][1] += weight*n1*n1*r*DComplex(0.0,kappa);
       turkelM[1][2] += weight*n1*n2*r*DComplex(0.0,kappa);
       turkelM[2][0] = turkelM[0][2];
       turkelM[2][1] = turkelM[1][2];
       turkelM[2][2] += weight*n2*n2*r*DComplex(0.0,kappa);

     }
     break;
   case 1:

      // BTL 1 for an ellipsoid of revolution parameterized by
      //       r(t) = surfR0*cos(t)               (surfR0>=0.0)
      //       z(t) = - surfZ0*sin(t)             (surfZ0>=0.0)

      if (fabs(surfR0*surfZ0)<1e-12) {
        fprintf(stderr,"Inconsistent parameters for the ellipse. Aborting\n");
        exit(1);
      }

      if (fabs(x[0])<1e-12) {
        if (y[0]>=0.0)
          t0 = -0.5*M_PI;
        else
          t0 = +0.5*M_PI;
      }
      else {
        t0 = -atan(y[0]*surfR0/(x[0]*surfZ0));
      }

      if (fabs(x[1])<1e-12) {
        if (y[1]>=0.0)
          t1 = -0.5*M_PI;
        else
          t1 = +0.5*M_PI;
      }
      else {
        t1 = -atan(y[1]*surfR0/(x[1]*surfZ0));
      }

      length = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
  
      for (i=0; i<3;++i) {

         double xi;
         double weight;

         switch (i) {
           case 0:
             xi = -sqrt(0.6);
             weight = 5.0/9.0;
             break;
           case 1:
             xi = 0.0;
             weight = 8.0/9.0;
             break;
           case 2:
             xi = sqrt(0.6);
             weight = 5.0/9.0;
             break;
         }

         weight *= 0.5*length;

         double r = 0.5*(x[0]+x[1])+0.5*xi*(x[1]-x[0]);

         double t = 0.5*(t0+t1)+0.5*xi*(t1-t0);
         double H;
         double tmp1 = pow(surfR0*sin(t),2.0) + pow(surfZ0*cos(t),2.0);
         double tmp2 = sqrt(tmp1);
         H = 0.5*surfZ0/(surfR0*tmp2) + 0.5*surfR0*surfZ0/(tmp1*tmp2);

         double n0 = 0.5*(xi-1)*xi;
         double n1 = 0.5*(1+xi)*xi;
         double n2 = 1-xi*xi;

         turkelM[0][0] += weight*n0*n0*r*DComplex(-H,kappa);
         turkelM[0][1] += weight*n0*n1*r*DComplex(-H,kappa);
         turkelM[0][2] += weight*n0*n2*r*DComplex(-H,kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*n1*n1*r*DComplex(-H,kappa);
         turkelM[1][2] += weight*n1*n2*r*DComplex(-H,kappa);
         turkelM[2][0] = turkelM[0][2];
         turkelM[2][1] = turkelM[1][2];
         turkelM[2][2] += weight*n2*n2*r*DComplex(-H,kappa);

      }
      break;
    case -1:

      // BTL 1 for a cylinder of revolution parameterized by surfR0, surfZ0
      // The corner is replaced by a circle. We assume that the 2 edges
      // on the corner have same length.

      if (fabs(surfR0*surfZ0)<1e-12) {
        fprintf(stderr,"Inconsistent parameters for the cylinder. Aborting\n");
        exit(1);
      }

      length = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));

      for (i=0; i<3;++i) {

         double xi;
         double weight;

         switch (i) {
           case 0:
             xi = -sqrt(0.6);
             weight = 5.0/9.0;
             break;
           case 1:
             xi = 0.0;
             weight = 8.0/9.0;
             break;
           case 2:
             xi = sqrt(0.6);
             weight = 5.0/9.0;
             break;
         }

         double H;

         if (fabs(x[0]-x[1])<1.0e-5) {
           H = 0.5/surfR0;
           if (fabs(fabs(y[0])-surfZ0)<1.0e-5) {
             double st = sin(M_PI*0.125*(3+xi));
             H = 0.5/length + 0.5*st/(surfR0-length*(1-st));
           }
           if (fabs(fabs(y[1])-surfZ0)<1.0e-5) {
             double st = sin(M_PI*0.125*(3-xi));
             H = 0.5/length + 0.5*st/(surfR0-length*(1-st));
           }
         }
         else {
           H = 0.0;
           if (fabs(x[0]-surfR0)<1.0e-5) {
             double st = sin(M_PI*0.125*(1.0-xi));
             H = 0.5/length + 0.5*st/(surfR0-length*(1-st));
           }
           if (fabs(x[1]-surfR0)<1.0e-5) {
             double st = sin(M_PI*0.125*(1.0+xi));
             H = 0.5/length + 0.5*st/(surfR0-length*(1-st));
           }
         }

         weight *= 0.5*length;

         double r = 0.5*(x[0]+x[1])+0.5*xi*(x[1]-x[0]);

         double n0 = 0.5*(xi-1)*xi;
         double n1 = 0.5*(1+xi)*xi;
         double n2 = 1-xi*xi;

         turkelM[0][0] += weight*n0*n0*r*DComplex(-H,kappa);
         turkelM[0][1] += weight*n0*n1*r*DComplex(-H,kappa);
         turkelM[0][2] += weight*n0*n2*r*DComplex(-H,kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*n1*n1*r*DComplex(-H,kappa);
         turkelM[1][2] += weight*n1*n2*r*DComplex(-H,kappa);
         turkelM[2][0] = turkelM[0][2];
         turkelM[2][1] = turkelM[1][2];
         turkelM[2][2] += weight*n2*n2*r*DComplex(-H,kappa);

      }
      break;
    case 2:

      // BTL 2 for an ellipsoid of revolution parameterized by
      //       r(t) = surfR0*cos(t)               (surfR0>=0.0)
      //       z(t) = - surfZ0*sin(t)             (surfZ0>=0.0)

      if (fabs(surfR0*surfZ0)<1e-12) {
        fprintf(stderr,"Inconsistent parameters for the ellipse. Aborting\n");
        exit(1);
      }

      if (fabs(x[0])<1e-12) {
        if (y[0]>=0.0)
          t0 = -0.5*M_PI;
        else
          t0 = +0.5*M_PI;
      }
      else {
        t0 = -atan(y[0]*surfR0/(x[0]*surfZ0));
      }

      if (fabs(x[1])<1e-12) {
        if (y[1]>=0.0)
          t1 = -0.5*M_PI;
        else
          t1 = +0.5*M_PI;
      }
      else {
        t1 = -atan(y[1]*surfR0/(x[1]*surfZ0));
      }

      length = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));

      for (i=0; i<3;++i) {

         double xi;
         double weight;

         switch (i) {
           case 0:
             xi = -sqrt(0.6);
             weight = 5.0/9.0;
             break;
           case 1:
             xi = 0.0;
             weight = 8.0/9.0;
             break;
           case 2:
             xi = sqrt(0.6);
             weight = 5.0/9.0;
             break;
         }

         weight *= 0.5*length;

         double r = 0.5*(x[0]+x[1])+0.5*xi*(x[1]-x[0]);

         double t = 0.5*(t0+t1)+0.5*xi*(t1-t0);
         double H, K;
         double R1, R2;
         double Hdot, Hdotdot;
         double tmp1 = pow(surfR0*sin(t),2.0) + pow(surfZ0*cos(t),2.0);
         double tmp2 = sqrt(tmp1);

         R1 = surfZ0/(surfR0*tmp2);
         R2 = surfR0*surfZ0/(tmp1*tmp2);
         H = 0.5*R1 + 0.5*R2;
         K = pow(surfZ0/tmp1,2.0);

         Hdot = (surfZ0-surfR0)*(surfZ0+surfR0)*sin(2*t)*
                (0.25*surfZ0/(surfR0*tmp1*tmp2) +
                 0.75*surfR0*surfZ0/(tmp1*tmp1*tmp2));

         Hdotdot = 2*(surfZ0-surfR0)*(surfZ0+surfR0)*cos(2*t)*
                   (0.25*surfZ0/(surfR0*tmp1*tmp2) +
                   0.75*surfR0*surfZ0/(tmp1*tmp1*tmp2));
         Hdotdot += (surfZ0-surfR0)*(surfZ0+surfR0)*sin(2*t)*
                   (0.25*surfZ0/(surfR0) + 0.75*surfR0*surfZ0/tmp1)*
                   1.5*(surfZ0-surfR0)*(surfZ0+surfR0)*sin(2*t)*
                   pow(tmp1,-2.5);
         Hdotdot += (surfZ0-surfR0)*(surfZ0+surfR0)*sin(2*t)*
                   (0.0 + 0.75*surfR0*surfZ0*sin(2*t)*(surfZ0-surfR0)*
                   (surfZ0+surfR0)*pow(tmp1,-2.0))*pow(tmp1,-1.5);

         double n0 = 0.5*(xi-1)*xi;
         double n1 = 0.5*(1+xi)*xi;
         double n2 = 1-xi*xi;

         DComplex coeff;
         coeff = DComplex(-H,kappa) - (K-H*H)*DComplex(0.0, 0.5/kappa)/
                 DComplex(1.0,2.0*H/kappa);

         turkelM[0][0] += weight*n0*n0*r*coeff;
         turkelM[0][1] += weight*n0*n1*r*coeff;
         turkelM[0][2] += weight*n0*n2*r*coeff;
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*n1*n1*r*coeff;
         turkelM[1][2] += weight*n1*n2*r*coeff;
         turkelM[2][0] = turkelM[0][2];
         turkelM[2][1] = turkelM[1][2];
         turkelM[2][2] += weight*n2*n2*r*coeff;

         double Hlaplace;
         Hlaplace = Hdot*(surfZ0-surfR0)*(surfZ0+surfR0)*0.5*sin(2*t)/
                    (tmp1*tmp1);
         Hlaplace += Hdotdot/tmp1;

         turkelM[0][0] += weight*n0*n0*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[0][1] += weight*n0*n1*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[0][2] += weight*n0*n2*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*n1*n1*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[1][2] += weight*n1*n2*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[2][0] = turkelM[0][2];
         turkelM[2][1] = turkelM[1][2];
         turkelM[2][2] += weight*n2*n2*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);

         turkelM[0][0] += weight*n0*n0*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);
         turkelM[0][1] += weight*n0*n1*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);
         turkelM[0][2] += weight*n0*n2*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*n1*n1*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);
         turkelM[1][2] += weight*n1*n2*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);
         turkelM[2][0] = turkelM[0][2];
         turkelM[2][1] = turkelM[1][2];
         turkelM[2][2] += weight*n2*n2*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);

         double dn0 = 0.5*(2.0*xi-1.0);
         double dn1 = 0.5*(2.0*xi+1.0);
         double dn2 = -2.0*xi;

         turkelM[0][0] += weight*r*dn0*dn0*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);
         turkelM[0][1] += weight*r*dn0*dn1*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);
         turkelM[0][2] += weight*r*dn0*dn2*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*r*dn1*dn1*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);
         turkelM[1][2] += weight*r*dn1*dn2*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);
         turkelM[2][0] = turkelM[0][2];
         turkelM[2][1] = turkelM[1][2];
         turkelM[2][2] += weight*r*dn2*dn2*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);

      }
      break;
    case -2:

      // BTL 2 for a cylinder of revolution parameterized by surfR0, surfZ0
      // The corner is replaced by a circle. We assume that the 2 edges
      // on the corner have same length.

      if (fabs(surfR0*surfZ0)<1e-12) {
        fprintf(stderr,"Inconsistent parameters for the ellipse. Aborting\n");
        exit(1);
      }

      length = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));

      for (i=0; i<3;++i) {

         double xi;
         double weight;

         switch (i) {
           case 0:
             xi = -sqrt(0.6);
             weight = 5.0/9.0;
             break;
           case 1:
             xi = 0.0;
             weight = 8.0/9.0;
             break;
           case 2:
             xi = sqrt(0.6);
             weight = 5.0/9.0;
             break;
         }

         weight *= 0.5*length;

         double r = 0.5*(x[0]+x[1])+0.5*xi*(x[1]-x[0]);

         double H, K;
         double R1, R2;
         double Hdot = 0.0;
         double Hdotdot = 0.0;
         double Hlaplace;

         if (fabs(x[0]-x[1])<1.0e-4) {
           R1 = 1.0/surfR0;
           R2 = 0.0;
           if (fabs(fabs(y[0])-surfZ0)<1.0e-5) {
             double st = sin(M_PI*0.125*(3+xi));
             double ct = cos(M_PI*0.125*(3+xi));
             R1 = st/(surfR0-length*(1-st));
             R2 = 1.0/length;
             Hdotdot = -0.5*(2*length*ct*ct + st*(surfR0-length*(1-st)))/
                       pow(surfR0-length*(1-st), 3.0);
             Hdotdot *= (surfR0-length)*pow(1.0/length, 2.0);
           }
           if (fabs(fabs(y[1])-surfZ0)<1.0e-5) {
             double st = sin(M_PI*0.125*(3-xi));
             double ct = cos(M_PI*0.125*(3-xi));
             R1 = st/(surfR0-length*(1-st));
             R2 = 1.0/length;
             Hdotdot = -0.5*(2*length*ct*ct + st*(surfR0-length*(1-st)))/
                       pow(surfR0-length*(1-st), 3.0);
             Hdotdot *= (surfR0-length)*pow(1.0/length, 2.0);
           }
           H = 0.5*R1 + 0.5*R2;
           K = R1*R2;
         }
         else {
           R1 = 0.0;
           R2 = 0.0;
           if (fabs(x[0]-surfR0)<1.0e-5) {
             double st = sin(M_PI*0.125*(1.0-xi));
             double ct = cos(M_PI*0.125*(1.0-xi));
             R1 = st/(surfR0-length*(1-st));
             R2 = 1.0/length;
             Hdotdot = -0.5*(2*length*ct*ct + st*(surfR0-length*(1-st)))/
                       pow(surfR0-length*(1-st), 3.0);
             Hdotdot *= (surfR0-length)*pow(1.0/length, 2.0);
           }
           if (fabs(x[1]-surfR0)<1.0e-5) {
             double st = sin(M_PI*0.125*(1+xi));
             double ct = cos(M_PI*0.125*(1+xi));
             R1 = st/(surfR0-length*(1-st));
             R2 = 1.0/length;
             Hdotdot = -0.5*(2*length*ct*ct + st*(surfR0-length*(1-st)))/
                       pow(surfR0-length*(1-st), 3.0);
             Hdotdot *= (surfR0-length)*pow(1.0/length, 2.0);
           }
           H = 0.5*R1 + 0.5*R2;
           K = R1*R2;
         }

         double n0 = 0.5*(xi-1)*xi;
         double n1 = 0.5*(1+xi)*xi;
         double n2 = 1-xi*xi;

         DComplex coeff;
         coeff = DComplex(-H,kappa) - (K-H*H)*DComplex(0.0, 0.5/kappa)/
                 DComplex(1.0,2.0*H/kappa);

         turkelM[0][0] += weight*n0*n0*r*coeff;
         turkelM[0][1] += weight*n0*n1*r*coeff;
         turkelM[0][2] += weight*n0*n2*r*coeff;
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*n1*n1*r*coeff;
         turkelM[1][2] += weight*n1*n2*r*coeff;
         turkelM[2][0] = turkelM[0][2];
         turkelM[2][1] = turkelM[1][2];
         turkelM[2][2] += weight*n2*n2*r*coeff;

         Hlaplace = Hdotdot - Hdot;

         turkelM[0][0] += weight*n0*n0*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[0][1] += weight*n0*n1*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[0][2] += weight*n0*n2*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*n1*n1*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[1][2] += weight*n1*n2*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[2][0] = turkelM[0][2];
         turkelM[2][1] = turkelM[1][2];
         turkelM[2][2] += weight*n2*n2*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);

         turkelM[0][0] += weight*n0*n0*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);
         turkelM[0][1] += weight*n0*n1*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);
         turkelM[0][2] += weight*n0*n2*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*n1*n1*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);
         turkelM[1][2] += weight*n1*n2*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);
         turkelM[2][0] = turkelM[0][2];
         turkelM[2][1] = turkelM[1][2];
         turkelM[2][2] += weight*n2*n2*r*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);

         double dn0 = 0.5*(2.0*xi-1.0);
         double dn1 = 0.5*(2.0*xi+1.0);
         double dn2 = -2.0*xi;

         turkelM[0][0] += weight*r*dn0*dn0*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);
         turkelM[0][1] += weight*r*dn0*dn1*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);
         turkelM[0][2] += weight*r*dn0*dn2*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*r*dn1*dn1*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);
         turkelM[1][2] += weight*r*dn1*dn2*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);
         turkelM[2][0] = turkelM[0][2];
         turkelM[2][1] = turkelM[1][2];
         turkelM[2][2] += weight*r*dn2*dn2*
                          (4.0/(pow(x[1]-x[0],2.0) + pow(y[1]-y[0],2.0)))/
                          DComplex(-2*R2,2*kappa);

      }
      break;
 }

 return turkelM;

}


FullSquareMatrix 
Line2AxiSommer::interfMatrixConsistent(CoordSet &cs) {

 return interfMatrixConsistent(cs, new double[3*3]);

}


FullSquareMatrix 
Line2AxiSommer::interfMatrixConsistent(CoordSet &cs, double *d) {

 Node nd1 = cs.getNode(nn[0]);
 Node nd2 = cs.getNode(nn[1]);
 Node nd3 = cs.getNode(nn[2]);

 double x[3], y[3];

 x[0] = nd1.x; y[0] = nd1.y;
 x[1] = nd2.x; y[1] = nd2.y;
 x[2] = nd3.x; y[2] = nd3.y;

 FullSquareMatrix interfM(3,d);
 interfM.zero();

 int i;

 for (i=0; i<3;++i) {
   double xi;
   double weight;

   switch (i) {
     case 0:
       xi = -sqrt(0.6);
       weight = 5.0/9.0;
       break;
     case 1:
       xi = 0.0;
       weight = 8.0/9.0;
       break;
     case 2: 
       xi = sqrt(0.6);
       weight = 5.0/9.0;
       break;
   }

   weight *= 0.5*sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));

   double r = 0.5*(x[0]+x[1])+0.5*xi*(x[1]-x[0]);

   double n0 = 0.5*(xi-1)*xi;
   double n1 = 0.5*(1+xi)*xi;
   double n2 = 1-xi*xi;

   interfM[0][0] += weight*n0*n0*r;
   interfM[0][1] += weight*n0*n1*r;
   interfM[0][2] += weight*n0*n2*r;
   interfM[1][0] = interfM[0][1];
   interfM[1][1] += weight*n1*n1*r;
   interfM[1][2] += weight*n1*n2*r;
   interfM[2][0] = interfM[0][2];
   interfM[2][1] = interfM[1][2];
   interfM[2][2] += weight*n2*n2*r;

 }

 return interfM;

}


FullSquareMatrix 
Line2AxiSommer::interfMatrixLumped(CoordSet &cs) {

 return interfMatrixLumped(cs, new double[3*3]);

}


FullSquareMatrix 
Line2AxiSommer::interfMatrixLumped(CoordSet &cs, double *d) {

 FullSquareMatrix interfM = interfMatrixConsistent(cs,d);

 interfM[0][0] += interfM[0][1] + interfM[0][2];
 interfM[0][1] = 0.0;
 interfM[0][2] = 0.0;
 interfM[1][1] += interfM[1][0] + interfM[1][2];
 interfM[1][0] = 0.0;
 interfM[1][2] = 0.0;
 interfM[2][2] += interfM[2][0] + interfM[2][1];
 interfM[2][0] = 0.0;
 interfM[2][1] = 0.0;

 return interfM;

}


int Line2AxiSommer::numDofs() {

 return 3;

}


int* Line2AxiSommer::dofs(DofSetArray &dsa, int *p) {

 if (p == 0) p = new int[3];

 dsa.number(nn[0],DofSet::Helm , p);
 dsa.number(nn[1],DofSet::Helm , p+1);
 dsa.number(nn[2],DofSet::Helm , p+2);

 return p;

}


int* Line2AxiSommer::nodes(int *p) {

 if (p == 0) 
   p = new int[3];

 p[0] = nn[0];
 p[1] = nn[1];
 p[2] = nn[2];

 return p;

}


int Line2AxiSommer::numNodes() {

 return 3;

}


void 
Line2AxiSommer::ffpAxiNeum(int ndir, DComplex *ffp, CoordSet& cs, DComplex **u,
                      double k, double (*dir)[3], double* idir, int numMode) {

 int i, j, mode;

 Node nd1 = cs.getNode(nn[0]);
 Node nd2 = cs.getNode(nn[1]);
 Node nd3 = cs.getNode(nn[2]);

 double x[3], y[3];
 x[0] = nd1.x; y[0] = nd1.y;
 x[1] = nd2.x; y[1] = nd2.y;
 x[2] = nd3.x; y[2] = nd3.y;

 double nx, ny, l;
 nx = -y[1] + y[0];
 ny = x[1] - x[0];
 l = sqrt(nx*nx+ny*ny);
 nx = nx/l;
 ny = ny/l;

 double &dx = idir[0],  &dy = idir[1], &dz = idir[2];

 DComplex expDir;
 DComplex expNeu;

 DComplex *dudn = (DComplex*) dbg_alloca(sizeof(DComplex)*3*numMode);

 for (i=0; i<3; ++i) {

   DComplex C = DComplex(-cos(k*y[i]*dz), -sin(k*y[i]*dz));

   double rho = fabs(x[i])*sqrt(dx*dx + dy*dy);
   rho *= (rho>1e-15);
   double ratio = (rho<1e-15) ? 1.0 : dx*x[i]/rho;
   ratio = (ratio>1.0)  ?  1.0 : ratio;
   ratio = (ratio<-1.0) ? -1.0 : ratio;
   double phi;
   phi = acos(ratio);
   if (dy<0.0)
     phi *= -1.0;
     
   for (mode=0; mode<numMode; ++mode) 
     dudn[3*mode+i] = coefExpNeu(mode, k, rho, phi, nx, ny, dx, dy, dz) * C;

 }

 double d[9];
 FullSquareMatrix mass=interfMatrixConsistent(cs,d);

 for (int iDir = 0; iDir < ndir; iDir++) {

   double &dirx = dir[iDir][0],  &diry = dir[iDir][1], &dirz = dir[iDir][2];

   for (i=0; i<3; ++i) {

     DComplex C = DComplex(cos(k*y[i]*dirz),-sin(k*y[i]*dirz));

     double rho = fabs(x[i])*sqrt(dirx*dirx + diry*diry);
     rho *= (rho>1e-15);
     double ratio = (rho<1e-15) ? 1.0 : dirx*x[i]/rho;
     ratio = (ratio>1.0)  ?  1.0 : ratio;
     ratio = (ratio<-1.0) ? -1.0 : ratio;
     double phi;
     phi = acos(ratio);
     if (diry<0.0)
       phi *= -1.0;
     
     for (mode=0; mode<numMode; ++mode)  {

       expNeu = coefExpNeu(mode, k, -rho, phi, nx, ny, dirx, diry, dirz) * C;

       if ((real(expNeu)!=0.0) || (imag(expNeu)!=0.0)) {
         if (mode==0) {
           expNeu *= DComplex(0.5, 0.0);
         }
         else {
           expNeu *= DComplex(0.25, 0.0);
         }
         for (j=0; j<3; ++j) {
           ffp[2*iDir] += u[mode][j] * mass[j][i] * expNeu;
           ffp[2*iDir+1] += u[mode][j]*mass[j][i]*
                            DComplex(real(expNeu),-imag(expNeu));
         }
       }

       expDir = coefExpDir(mode, k, -rho, phi) * C;
  
       if ((real(expDir)!=0.0) || (imag(expDir)!=0.0)) {
         if (mode==0) {
           expDir *= DComplex(0.5, 0.0);
         }
         else {
           expDir *= DComplex(0.25, 0.0);
         }
         for (j=0; j<3; ++j) {
           ffp[2*iDir] += dudn[3*mode+j] * mass[j][i] * expDir;
           ffp[2*iDir+1] += dudn[3*mode+j] * mass[j][i] *
                            DComplex(real(expDir),-imag(expDir));
         }
       }

     }

   }

 }

}


void 
Line2AxiSommer::ffpAxiDir(int ndir, DComplex *ffp, CoordSet& cs, DComplex **u,
                      double k, double (*dir)[3], double* idir, int numMode) {

 int i, j, mode;

 Node nd1 = cs.getNode(nn[0]);
 Node nd2 = cs.getNode(nn[1]);
 Node nd3 = cs.getNode(nn[2]);

 double x[3], y[3];
 x[0] = nd1.x; y[0] = nd1.y;
 x[1] = nd2.x; y[1] = nd2.y;
 x[2] = nd3.x; y[2] = nd3.y;

 double nx, ny, l;
 nx = -y[1] + y[0];
 ny = x[1] - x[0];
 l = sqrt(nx*nx+ny*ny);
 nx = nx/l;
 ny = ny/l;

 DComplex expNeuEl;

 double d[9];
 FullSquareMatrix mass=interfMatrixConsistent(cs,d);

 DComplex tmp = DComplex(0.25, 0.0);

 for (int iDir = 0; iDir < ndir; iDir++) {

   double &dirx = dir[iDir][0],  &diry = dir[iDir][1], &dirz = dir[iDir][2];

   for (i=0; i<3; ++i) {

     DComplex C = DComplex(cos(k*y[i]*dirz),-sin(k*y[i]*dirz));

     double rho = fabs(x[i])*sqrt(dirx*dirx + diry*diry);
     rho *= (rho>1e-15);
     double ratio = (rho<1e-15) ? 1.0 : dirx*x[i]/rho;
     ratio = (ratio>1.0)  ?  1.0 : ratio;
     ratio = (ratio<-1.0) ? -1.0 : ratio;
     double phi;
     phi = acos(ratio);
     if (diry<0.0)
       phi *= -1.0;
     
     for (mode=0; mode<numMode; ++mode)  {

       expNeuEl = coefExpNeu(mode, k, -rho, phi, nx, ny, dirx, diry, 
                         dirz) * C;
       
       if ((real(expNeuEl)==0.0) && (imag(expNeuEl)==0.0))
         continue;

       if (mode==0) 
         expNeuEl *= DComplex(0.5, 0.0);
       else
         expNeuEl *= DComplex(0.25, 0.0);

       for (j=0; j<3; ++j) {
         ffp[2*iDir] += u[mode][j] * mass[j][i] * expNeuEl;
         ffp[2*iDir+1] += u[mode][j] * mass[j][i] * 
                          DComplex(real(expNeuEl), -imag(expNeuEl));
       }

     }

   }

 }

}
