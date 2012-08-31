#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Utils.d/MyComplex.h>
#include <HelmAxi.d/coefFourier.h>
#include <HelmAxi.d/FourierHelmBCs.h>
#include <HelmAxi.d/LineAxiSommer.h>
#include <HelmAxi.d/Line2AxiSommer.h>


FourierHelmBCs * fourHelmBC = 0;


ComplexFDBC zeroFDBC = FDBC(-1);
ComplexFNBC zeroFNBC = FNBC(-1, -1, -1);


FourierHelmBCs::FourierHelmBCs() : 
          dirBC(zeroFDBC), neuBC(zeroFNBC), somBC(NULL) {

  numDir = 0;
  numNeu = 0;
  numSomm = 0;
  numModes = 0;
  numSlices = 0;
  surfR = 0.0;
  surfZ = 0.0;

}


void FourierHelmBCs::addDirichlet(ComplexFDBC &boundary) {
  dirBC[numDir++] = boundary;
}


void FourierHelmBCs::addNeuman(ComplexFNBC &boundary) {
  neuBC[numNeu++] = boundary;
}


void FourierHelmBCs::addSommer(SommerElement *ele) {

  somBC[numSomm++] = ele;

  LineAxiSommer *ele1=dynamic_cast<LineAxiSommer *>(ele);
  Line2AxiSommer *ele2 = dynamic_cast<Line2AxiSommer *>(ele);

  if (ele1 != 0) {
    ele1->setSurf(surfR, surfZ);
    ele1->setType(somType);
  }

  if (ele2 != 0) {
    ele2->setSurf(surfR, surfZ);
    ele2->setType(somType);
  }

}


void FourierHelmBCs::setModes(int mode) {
  numModes = mode;
}


void FourierHelmBCs::setSlices(int slice) {
  numSlices = slice;
}


void FourierHelmBCs::setSurf(double aR, double aZ) {
  surfR = aR;
  surfZ = aZ;
}


void FourierHelmBCs::setConst(DComplex C) {
  cstant = C;
}


void FourierHelmBCs::setSomType(int type) {
  somType = type;
}


void FourierHelmBCs::setDir(double dx, double dy, double dz) {
  dirX = dx;
  dirY = dy;
  dirZ = dz;
}


DComplex 
FourierHelmBCs::Dir2Four(int mode, int pos, double k, double x, double y) {

  DComplex result;
  DComplex Ez(cos(k*y*dirZ),sin(k*y*dirZ));

  double rho = sqrt(x*x*(dirX*dirX + dirY*dirY));

  double ratio = (rho<1e-15) ? 1.0 : dirX*x/rho;
  ratio = (ratio>1.0)  ?  1.0 : ratio;
  ratio = (ratio<-1.0) ? -1.0 : ratio;

  double phi;
  phi = acos(ratio);
  if (dirY<0.0)
    phi *= -1.0;

  result = coefExpDir(mode, k, rho, phi);
  result *= Ez;
  result *= cstant;

  return result; 

} 


DComplex
FourierHelmBCs::Neu2Four(int mode, int pos, double k, double x, double y, 
                         double nr, double nz) {

  DComplex result;
  DComplex Ez(cos(k*y*dirZ),sin(k*y*dirZ));

  double rho=sqrt(x*x*(dirX*dirX + dirY*dirY));

  double ratio = (rho<1e-15) ? 1.0 : dirX*x/rho;
  ratio = (ratio>1.0)  ?  1.0 : ratio;
  ratio = (ratio<-1.0) ? -1.0 : ratio;

  double phi;
  phi = acos(ratio);
  if (dirY<0.0)
    phi *= -1.0;

  result = coefExpNeu(mode, k, rho, phi, nr, nz, dirX, dirY, dirZ);
  result *= Ez;
  result *= cstant;

  return result;

}


void 
FourierHelmBCs::IntNeu2(int mode, int pos, double k, double *x, double *y, 
                       double nx, double ny, DComplex *T) {

  DComplex V1(0.0, 0.0);
  DComplex V2(0.0, 0.0);

  double rho = sqrt(nx*nx + ny*ny);

  if (rho==0.0) {
    int one=1;
    fprintf(stderr,"Error in AXIHNEU - node appearing simultaneously."
                   " Aborting\n");
    exit(one);
  }

  rho = 1/rho;
 
  nx *= rho;
  ny *= rho;

  for (int i=0; i<4;++i) {
     double xi;
     double weight;

     switch (i) {
       case 0:
         xi =     -0.861136311594053;
         weight =  0.347854845137454;
         break; 
       case 1:
         xi =     -0.339981043584856;
         weight =  0.652145154862546;
         break;
       case 2:
         xi =      0.339981043584856;
         weight =  0.652145154862546;
         break;
       case 3:
         xi =      0.861136311594053;
         weight =  0.347854845137454;
         break;
     }

     weight *= 0.5*sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));

     double r = 0.5*(x[1]+x[0])+0.5*xi*(x[1]-x[0]);
     double z = 0.5*(y[1]+y[0])+0.5*xi*(y[1]-y[0]);

     weight *= r;

     DComplex coeff = Neu2Four(mode,pos,k,r,z,nx,ny);
     
     double n = 0.5*(1-xi);
     V1 += weight*n*coeff;
     V2 += weight*(1-n)*coeff;
  }

  V1 *= M_PI;
  V2 *= M_PI;

  if (mode==0) {
    V1 *= 2;
    V2 *= 2;
  }  

  T[0] += V1;
  T[1] += V2;

}


void 
FourierHelmBCs::IntNeu3(int mode, int pos, double k, double *x, double *y, 
                       double nx, double ny, DComplex *T) {

  DComplex V1(0.0, 0.0);
  DComplex V2(0.0, 0.0);
  DComplex V3(0.0, 0.0);

  double rho = sqrt(nx*nx + ny*ny);

  if (rho==0.0) {
    int one=1;
    fprintf(stderr,"Error in AXIHNEU - node appearing simultaneously."
                   " Aborting\n");
    exit(one);
  }

  rho = 1/rho;
 
  nx *= rho;
  ny *= rho;

  for (int i=0; i<4;++i) {
     double xi;
     double weight;

     switch (i) {
       case 0:
         xi =     -0.861136311594053;
         weight =  0.347854845137454;
         break; 
       case 1:
         xi =     -0.339981043584856;
         weight =  0.652145154862546;
         break;
       case 2:
         xi =      0.339981043584856;
         weight =  0.652145154862546;
         break;
       case 3:
         xi =      0.861136311594053;
         weight =  0.347854845137454;
         break;
     }

     weight *= 0.5*sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));

     double r = 0.5*(x[1]+x[0])+0.5*xi*(x[1]-x[0]);
     double z = 0.5*(y[1]+y[0])+0.5*xi*(y[1]-y[0]);

     weight *= r;

     DComplex coeff = Neu2Four(mode,pos,k,r,z,nx,ny);
     
     double n1 = 0.5*(xi-1)*xi;
     double n2 = 0.5*(1+xi)*xi;
     double n3 = (1-xi*xi);
     V1 += weight*n1*coeff;
     V2 += weight*n2*coeff;
     V3 += weight*n3*coeff;
  }

  V1 *= M_PI;
  V2 *= M_PI;
  V3 *= M_PI;

  if (mode==0) {
    V1 *= 2;
    V2 *= 2;
    V3 *= 2;
  }  

  T[0] += V1;
  T[1] += V2;
  T[2] += V3;

}


