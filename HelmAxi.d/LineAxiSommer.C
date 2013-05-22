#include <cstdlib>
#include <Utils.d/dbg_alloca.h>
#include <cmath>
#include <HelmAxi.d/coefFourier.h>
#include <HelmAxi.d/LineAxiSommer.h>
#include <Utils.d/MyComplex.h>
#include <Math.d/FullSquareMatrix.h>

#include <Timers.d/GetTime.h>

// Sommerfeld b.c. contribution for 2-D elements
// Note that the boundary condition will be always 
// on an element edge. 


LineAxiSommer::LineAxiSommer(int n1, int n2) {

 nn[0] = n1;
 nn[1] = n2;

}


void
LineAxiSommer::setType(int t) {

 type = t;

}


void
LineAxiSommer::setSurf(double aR, double aZ) {

 surfR0 = aR;
 surfZ0 = aZ;

}


void
LineAxiSommer::renum(int *table) {

 nn[0] = table[nn[0]];
 nn[1] = table[nn[1]];

}

void
LineAxiSommer::renum(EleRenumMap &table) {

 nn[0] = table[nn[0]];
 nn[1] = table[nn[1]];

}


FullSquareMatrix
LineAxiSommer::sommerMatrix(CoordSet &cs) {

 return sommerMatrix(cs,new double [2*2]);

}


FullSquareMatrix
LineAxiSommer::sommerMatrix(CoordSet &cs, double *d) {

 FullSquareMatrix sommerM(2,d);
 return sommerM;

}


FullSquareMatrixC
LineAxiSommer::turkelMatrix(CoordSet &cs, double kappa, int mode) {

 return turkelMatrix(cs,kappa,mode,new DComplex [2*2]);

}


FullSquareMatrixC
LineAxiSommer::turkelMatrix(CoordSet &cs, double kappa, int mode, 
               DComplex *d) {

 Node nd1 = cs.getNode(nn[0]);
 Node nd2 = cs.getNode(nn[1]);

 double x[2], y[2];
 double t0, t1;
 double length;

 x[0] = nd1.x; y[0] = nd1.y;
 x[1] = nd2.x; y[1] = nd2.y;

 FullSquareMatrixC turkelM(2,d);
 turkelM.zero();

 int i;

 switch (type) {
   default:
      fprintf(stderr,"BC parameter not available -> Sommerfeld\n");
   case 0:
      length = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
      for (i=0; i<2;++i) {
         double xi;
         double weight;

         switch (i) {
           case 0:
             xi = -1/sqrt(3.0);
             weight = 1;
             break;
           case 1:
             xi = 1/sqrt(3.0);
             weight = 1;
             break;
        }

        weight *= 0.5*length;

        double r = 0.5*(x[0]+x[1])+0.5*xi*(x[1]-x[0]);

        double n = 0.5*(1-xi);

        turkelM[0][0] += weight*n*n*r*DComplex(0.0,kappa);
        turkelM[0][1] += weight*n*(1-n)*r*DComplex(0.0,kappa);
        turkelM[1][0] = turkelM[0][1];
        turkelM[1][1] += weight*(1-n)*(1-n)*r*DComplex(0.0,kappa);
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

         double n = 0.5*(1-xi);

         turkelM[0][0] += weight*n*n*r*DComplex(-H,kappa);
         turkelM[0][1] += weight*n*(1-n)*r*DComplex(-H,kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*(1-n)*(1-n)*r*DComplex(-H,kappa);
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

         double n = 0.5*(1-xi);
   
         turkelM[0][0] += weight*n*n*r*DComplex(-H,kappa);
         turkelM[0][1] += weight*n*(1-n)*r*DComplex(-H,kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*n*n*r*DComplex(-H,kappa);

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

         double n = 0.5*(1-xi);

         DComplex coeff;
         coeff = DComplex(-H,kappa) - (K-H*H)*DComplex(0.0, 0.5/kappa)/
                 DComplex(1.0,2.0*H/kappa);

         turkelM[0][0] += weight*n*n*r*coeff;
         turkelM[0][1] += weight*n*(1-n)*r*coeff;
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*(1-n)*(1-n)*r*coeff;

         double Hlaplace;
         Hlaplace = Hdot*(surfZ0-surfR0)*(surfZ0+surfR0)*0.5*sin(2*t)/
                    (tmp1*tmp1);
         Hlaplace += Hdotdot/tmp1;

         turkelM[0][0] += weight*n*n*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[0][1] += weight*n*(1-n)*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*(1-n)*(1-n)*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);

         turkelM[0][0] += weight*n*n*mode*mode/DComplex(-2*R1*r,2*r*kappa);
         turkelM[0][1] += weight*n*(1-n)*mode*mode/DComplex(-2*R1*r,2*r*kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*(1-n)*(1-n)*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);

         turkelM[0][0] += weight*r*(1.0/(length*length))/
                          DComplex(-2*R2,2*kappa);
         turkelM[0][1] += weight*r*(-1.0/(length*length))/
                          DComplex(-2*R2,2*kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*r*(1.0/(length*length))/
                          DComplex(-2*R2,2*kappa);

      }

      break;
   case -2:

      // BTL 2 for a cylinder of revolution parameterized by surfR0, surfZ0
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

         double n = 0.5*(1-xi);

         DComplex coeff;
         coeff = DComplex(-H,kappa) - (K-H*H)*DComplex(0.0, 0.5/kappa)/
                 DComplex(1.0,2.0*H/kappa);

         turkelM[0][0] += weight*n*n*r*coeff;
         turkelM[0][1] += weight*n*(1-n)*r*coeff;
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*(1-n)*(1-n)*r*coeff;

         Hlaplace = Hdotdot - Hdot;

         turkelM[0][0] += weight*n*n*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[0][1] += weight*n*(1-n)*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*(1-n)*(1-n)*r*Hlaplace*
                          DComplex(0.25/(kappa*kappa), 0.0);

         turkelM[0][0] += weight*n*n*mode*mode/DComplex(-2*R1*r,2*r*kappa);
         turkelM[0][1] += weight*n*(1-n)*mode*mode/DComplex(-2*R1*r,2*r*kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*(1-n)*(1-n)*mode*mode/
                          DComplex(-2*R1*r,2*r*kappa);

         turkelM[0][0] += weight*r*(1.0/(length*length))/
                          DComplex(-2*R2,2*kappa);
         turkelM[0][1] += weight*r*(-1.0/(length*length))/
                          DComplex(-2*R2,2*kappa);
         turkelM[1][0] = turkelM[0][1];
         turkelM[1][1] += weight*r*(1.0/(length*length))/
                          DComplex(-2*R2,2*kappa);

      }

      break;
 }

 return turkelM;

}


FullSquareMatrix 
LineAxiSommer::interfMatrixConsistent(CoordSet &cs) {

  return interfMatrixConsistent(cs, new double[2*2]);

}


FullSquareMatrix 
LineAxiSommer::interfMatrixConsistent(CoordSet &cs, double *d) {

   Node nd1 = cs.getNode(nn[0]);
   Node nd2 = cs.getNode(nn[1]);

   double x[2], y[2];

   x[0] = nd1.x; y[0] = nd1.y;
   x[1] = nd2.x; y[1] = nd2.y;

   FullSquareMatrix interfM(2,d);
   interfM.zero();

   int i;

   for (i=0; i<2;++i) {
      double xi;
      double weight;

      switch (i) {
        case 0:
          xi = -1/sqrt(3.0);
          weight = 1;
          break;
        case 1:
          xi = 1/sqrt(3.0);
          weight = 1;
          break;
      }

      weight *= 0.5*sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));

      double r = 0.5*(x[0]+x[1])+0.5*xi*(x[1]-x[0]);

      double n = 0.5*(1-xi);

      interfM[0][0] += weight*n*n*r;
      interfM[0][1] += weight*n*(1-n)*r;
      interfM[1][0] = interfM[0][1];
      interfM[1][1] += weight*(1-n)*(1-n)*r;
   }

   return interfM;

}


FullSquareMatrix 
LineAxiSommer::interfMatrixLumped(CoordSet &cs) {

  return interfMatrixLumped(cs, new double[2*2]);

}


FullSquareMatrix 
LineAxiSommer::interfMatrixLumped(CoordSet &cs, double *d) {

 FullSquareMatrix interfM = interfMatrixConsistent(cs,d);

 interfM[0][0] += interfM[0][1];
 interfM[0][1] = 0.0;
 interfM[1][1] += interfM[1][0];
 interfM[1][0] = 0.0;

 return interfM;

}


int LineAxiSommer::numDofs() {

 return 2;

}


int* LineAxiSommer::dofs(DofSetArray &dsa, int *p) {

 if(p == 0) p = new int[2];

 dsa.number(nn[0],DofSet::Helm , p);
 dsa.number(nn[1],DofSet::Helm , p+1);

 return p;

}


int* LineAxiSommer::nodes(int *p) {

 if (p == 0) 
   p = new int[2];

 p[0] = nn[0];
 p[1] = nn[1];

 return p;

}


int LineAxiSommer::numNodes() {

 return 2;

}


void 
LineAxiSommer::ffpAxiNeum(int ndir, DComplex *ffp, CoordSet& cs, DComplex **u,
                      double k, double (*dir)[3], double* idir, int numMode) {

 int i, j, mode;

 Node nd1 = cs.getNode(nn[0]);
 Node nd2 = cs.getNode(nn[1]);

 double x[2], y[2];
 x[0] = nd1.x; y[0] = nd1.y;
 x[1] = nd2.x; y[1] = nd2.y;

 double nx, ny, l;
 nx = -y[1] + y[0];
 ny = x[1] - x[0];
 l = sqrt(nx*nx+ny*ny);
 nx = nx/l;
 ny = ny/l;

 double &dx = idir[0],  &dy = idir[1], &dz = idir[2];

 DComplex expDir;
 DComplex expNeu;

 DComplex *dudn = (DComplex*) dbg_alloca(sizeof(DComplex)*2*numMode);

 for (i=0; i<2; ++i) {

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
     dudn[2*mode+i] = coefExpNeu(mode, k, rho, phi, nx, ny, dx, dy, dz) * C;

 }

 double d[4];
 FullSquareMatrix mass=interfMatrixConsistent(cs,d);

 for (int iDir = 0; iDir < ndir; iDir++) {

   double &dirx = dir[iDir][0],  &diry = dir[iDir][1], &dirz = dir[iDir][2];

   for (i=0; i<2; ++i) {

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
         for (j=0; j<2; ++j) {
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
         for (j=0; j<2; ++j) {
           ffp[2*iDir] += dudn[2*mode+j] * mass[j][i] * expDir;
           ffp[2*iDir+1] += dudn[2*mode+j] * mass[j][i] *
                            DComplex(real(expDir),-imag(expDir));
         }
       }

     }

   }

 }

}


void 
LineAxiSommer::ffpAxiDir(int ndir, DComplex *ffp, CoordSet& cs, DComplex **u,
                      double k, double (*dir)[3], double* idir, int numMode) {

 int i, j, mode;

 Node nd1 = cs.getNode(nn[0]);
 Node nd2 = cs.getNode(nn[1]);

 double x[2], y[2];
 x[0] = nd1.x; y[0] = nd1.y;
 x[1] = nd2.x; y[1] = nd2.y;

 double nx, ny, l;
 nx = -y[1] + y[0];
 ny = x[1] - x[0];
 l = sqrt(nx*nx+ny*ny);
 nx = nx/l;
 ny = ny/l;

 DComplex expNeuEl;

 double d[4];
 FullSquareMatrix mass=interfMatrixConsistent(cs,d);

 DComplex tmp = DComplex(0.25, 0.0);

 for (int iDir = 0; iDir < ndir; iDir++) {

   double &dirx = dir[iDir][0],  &diry = dir[iDir][1], &dirz = dir[iDir][2];

   for (i=0; i<2; ++i) {

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

       for (j=0; j<2; ++j) {
         ffp[2*iDir] += u[mode][j] * mass[j][i] * expNeuEl;
         ffp[2*iDir+1] += u[mode][j] * mass[j][i] * 
                          DComplex(real(expNeuEl), -imag(expNeuEl));
       }

     }

   }

 }

}
