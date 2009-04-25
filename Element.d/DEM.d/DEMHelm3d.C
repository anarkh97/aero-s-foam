#include <stdio.h>
#include <stdlib.h>

#include <Element.d/DEM.d/DEMHelm3d.h>
#include <Element.d/Helm.d/IsoParamUtils.h>
#include <math.h>

class HelmDGMEMatricesFunction3d : public IntegFunctionA3d {
 int o;
 double kappa;
 int ndir;
 complex<double> *dirs;
 int nldir;
 int *nldirs;
 complex<double> *ldirs;
 complex<double> *K;
 complex<double> *L;
 complex<double> *PL;
 complex<double> *PP;
 complex<double> *PE;
 int arbflag;
 int fi;
 double *xsc;
 double *xc;
public:
 HelmDGMEMatricesFunction3d(int _o, double _kappa,
                      int _ndir, complex<double> *_dirs,
                      int _nldir, complex<double> *_ldirs,
                      double *_xsc, double *_xc, int _arbflag, int _fi,
                      complex<double> *_K, complex<double> *_L,
                      complex<double> *_PL,complex<double> *_PP,
                      complex<double> *_PE) {
   o = _o; kappa = _kappa; ndir = _ndir; dirs = _dirs; nldir = _nldir; 
   ldirs = _ldirs; K = _K; L = _L; PL = _PL; PP = _PP; PE =_PE;
   arbflag = _arbflag; fi = _fi;
   xsc = _xsc; xc = _xc;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   int oc = o*o*o;
  
   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1]) +
                                      dirs[i*3+2]*(x[2]-xc[2]) );

   double wc = w*sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);

   for(int j=0;j<nldir;j++) {
     complex<double> l = exp(ldirs[j*3+0]*(x[0]-xsc[0])+
                             ldirs[j*3+1]*(x[1]-xsc[1])+
                             ldirs[j*3+2]*(x[2]-xsc[2]) );
     l *= wc;
     for(int i=0;i<ndir;i++)
       L[j*ndir+i] += e[i]*l;
     if (PL!=0) for(int i=0;i<oc;i++)
       PL[j*oc+i] += N[i]*l;
   }

   complex<double> ikc =
     (arbflag==0)? 0.0: complex<double>(0.0, w*kappa*
                 sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]));
   for(int j=0;j<ndir;j++) {
     complex<double> ic = w*nsign*
              (cross[0]*dirs[j*3+0]+cross[1]*dirs[j*3+1]+cross[2]*dirs[j*3+2]);
     ic -= ikc;
     for(int i=j;i<ndir;i++)
       K[j*ndir+i] += ic*e[i]*e[j];
   }
// Sommerfeld part of bc
   if (PP!=0)  {
     for(int j=0;j<oc;j++) {
       for(int i=j;i<oc;i++) PP[j*oc+i] -= ikc*N[i]*N[j];
       for(int i=0;i<ndir;i++) PE[i*oc+j] -= ikc*N[j]*e[i];
     }
   }
   
   delete[] e;
 }
};


class HelmDGMELMatrixFunction3d : public IntegFunctionA3d {
 double kappa;
 int ndir;
 complex<double> *dirs;
 int nldir;
 int *nldirs;
 complex<double> *ldirs;
 complex<double> *L;
 int fi;
 double *xsc;
 double *xc;
public:
 HelmDGMELMatrixFunction3d(double _kappa,
                      int _ndir, complex<double> *_dirs,
                      int _nldir, complex<double> *_ldirs, double *_xsc, double *_xc,
                      complex<double> *_L, int _fi) {
   kappa = _kappa; ndir = _ndir; dirs = _dirs; nldir = _nldir; ldirs = _ldirs;
   L = _L; fi = _fi; xsc = _xsc; xc = _xc;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
  
   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1]) +
                                      dirs[i*3+2]*(x[2]-xc[2]) );

   double wc = w*sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
   for(int j=0;j<nldir;j++) {
     complex<double> l = exp(ldirs[j*3+0]*(x[0]-xsc[0])+
                             ldirs[j*3+1]*(x[1]-xsc[1])+
                             ldirs[j*3+2]*(x[2]-xsc[2]) );
     l *= wc;
     for(int i=0;i<ndir;i++)
       L[j*ndir+i] += e[i]*l;
   }
   delete[] e;
 }
};


class HelmDGMSomEEMatrixFunction3d : public IntegFunctionA3d {
 double kappa;
 int ndir;
 complex<double> *dirs;
 double *xc;
 complex<double> *K;
public:
 HelmDGMSomEEMatrixFunction3d(double _kappa,
                      int _ndir, complex<double> *_dirs, double *_xc,
                      complex<double> *_K) {
   kappa = _kappa; ndir = _ndir; dirs = _dirs; xc = _xc;
   K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1]) +
                                      dirs[i*3+2]*(x[2]-xc[2]) );

   complex<double> ikc = complex<double>(0.0, w*kappa*
           sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]));
   for(int j=0;j<ndir;j++)
     for(int i=0;i<ndir;i++)
       K[j*ndir+i] -= ikc*e[i]*e[j];
   delete[] e; 
 }
};


class HelmDGMEeMatrixFunction3d : IntegFunctionA3d {
 double kappa;
 int ndir;
 complex<double> *dirs;
 double *xc;
 complex<double> *K;
public:
 HelmDGMEeMatrixFunction3d(double _kappa,
                      int _ndir, complex<double> *_dirs, double *_xc,
                      complex<double> *_K) {
   kappa = _kappa; ndir = _ndir; dirs = _dirs; xc = _xc;
   K = _K;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   complex<double> *e = new complex<double>[2*ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1]) +
                                      dirs[i*3+2]*(x[2]-xc[2]) );

   complex<double> *ikc = e+ndir;
   for(int i=0;i<ndir;i++) 
     ikc[i] = w*nsign*
              (cross[0]*dirs[i*3+0]+cross[1]*dirs[i*3+1]+cross[2]*dirs[i*3+2]);

   for(int j=0;j<ndir;j++)
     for(int i=0;i<ndir;i++)
       K[j*ndir+i] += ikc[i]*e[i]*e[j];
   delete[] e;
 }
};


void HelmDGMEMatrices3d(int order, double *xyz,
                    int ndir, complex<double> *dirs,
                    int *nldirs, complex<double> *ldirs,
                    double kappa, int *sflags, double *xsc,
                    double *xc,
                    complex<double> *K, complex<double> *L) {


 IsoParamUtils *ipu = (order>0)? new IsoParamUtils(order):new IsoParamUtilsTetra(-order);
 int gorder = (order>0)?7:13;

 int nldir = nldirs[0] + nldirs[1] + nldirs[2] + nldirs[3];
 int nf = 4;
 if (order>0) {
   nldir += nldirs[4] + nldirs[5];
   nf = 6;
 }
 for(int i=0;i<nldir*ndir;i++) L[i] = 0.0;
 for(int i=0;i<ndir*ndir;i++) K[i] = 0.0;

// int os = ipu->getordersq();
// double *fxyz = new double[os*3];
// int *fi = new int[os];

 int c = 0;
 for(int faceindex=1;faceindex<=nf;faceindex++) {
//     ipu->faceindeces(faceindex,fi);
//     for(int i=0;i<os;i++) {
//       fxyz[i*3+0] = xyz[fi[i]*3+0];
//       fxyz[i*3+1] = xyz[fi[i]*3+1];
//       fxyz[i*3+2] = xyz[fi[i]*3+2];
//     }
     HelmDGMEMatricesFunction3d f(order,kappa,ndir,dirs,
                  nldirs[faceindex-1],ldirs+c*3,
                  xsc + (faceindex-1)*3, xc, sflags[faceindex-1], faceindex,
                  K,L+c*ndir,0,0,0);
     ipu->surfInt3d(xyz, faceindex, f,gorder);
     c += nldirs[faceindex-1];
     
/*
     HelmDGMEeMatrixFunction3d f(kappa,ndir,dirs,xc,K);
//     ipu.surfSurfInt3d(fxyz, f);
     ipu.surfInt3d(xyz, faceindex, f);
     if (nldirs[faceindex-1]!=0) {
       HelmDGMELMatrixFunction3d f(kappa,ndir,dirs,
                  nldirs[faceindex-1],ldirs+c*3, xsc + (faceindex-1)*3, xc,
                  L+c*ndir,faceindex);
//       ipu.surfSurfInt3d(fxyz, f);
       ipu.surfInt3d(xyz, faceindex, f);
       c += nldirs[faceindex-1];
     }
     if (sflags[faceindex-1]!=0) {
       HelmDGMSomEEMatrixFunction3d f(kappa,ndir,dirs,xc,K);
//       ipu.surfSurfInt3d(fxyz, f);
       ipu.surfInt3d(xyz, faceindex, f);
     }
*/
 }
 for(int j=0;j<ndir;j++) {
   for(int i=0;i<j;i++)
     K[j*ndir+i] = K[i*ndir+j];
 }
// delete[] fxyz;
// delete[] fi;
 delete ipu;
}


int *DGMHelm3d::faceCorners(int fi) {
 if (o>0) {
   int *fc = new int[4];
   int osq = o*o;
   int oc = osq*o;
   if (fi==1) {
     fc[0] = nn[1 - 1];
     fc[1] = nn[osq-o+1 - 1];
     fc[2] = nn[osq - 1];
     fc[3] = nn[o - 1];
   } else if (fi==2) {
     fc[0] = nn[1 - 1];
     fc[1] = nn[o - 1];
     fc[2] = nn[osq*(o-1)+o - 1];
     fc[3] = nn[osq*(o-1)+1 - 1];
   } else if (fi==3) {
     fc[0] = nn[o - 1];
     fc[1] = nn[osq - 1];
     fc[2] = nn[oc - 1];
     fc[3] = nn[osq*(o-1)+o - 1];
   } else if (fi==4) {
     fc[0] = nn[osq-o+1 - 1];
     fc[1] = nn[oc-o+1 - 1];
     fc[2] = nn[oc - 1];
     fc[3] = nn[osq - 1];
   } else if (fi==5) {
     fc[0] = nn[1 - 1];
     fc[1] = nn[osq*(o-1)+1 - 1];
     fc[2] = nn[oc-o+1 - 1];
     fc[3] = nn[osq-o+1 - 1];
   } else {
     fc[0] = nn[osq*(o-1)+1 - 1];
     fc[1] = nn[osq*(o-1)+o - 1];
     fc[2] = nn[oc - 1];
     fc[3] = nn[oc-o+1 - 1];
   }
   return fc;
 } else {
   int *fc = new int[3];
   int osq = ((-o)*(-o+1))/2;
   int oc = ((-o)*(-o+1)*(-o+2))/6;
   if (fi==1) {
     fc[0] = nn[-o - 1];
     fc[1] = nn[osq - 1];
     fc[2] = nn[oc - 1];
   } else if (fi==2) {
     fc[0] = nn[1 - 1];
     fc[1] = nn[oc - 1];
     fc[2] = nn[osq - 1];
   } else if (fi==3) {
     fc[0] = nn[1 - 1];
     fc[1] = nn[-o - 1];
     fc[2] = nn[oc - 1];
   } else {
     fc[0] = nn[1 - 1];
     fc[1] = nn[osq - 1];
     fc[2] = nn[-o - 1];
   }
   return fc;
 }
}


DGMHelm3d_6::DGMHelm3d_6(int _nnodes, int* nodenums) :
  DGMHelm3d(_nnodes,nodenums) {
 ndir = 6;
}

DGMHelm3d_6t::DGMHelm3d_6t(int _nnodes, int* nodenums) :
  DGMHelm3d(-_nnodes,nodenums) {
 ndir = 6;
}


DGMHelm3d_26::DGMHelm3d_26(int _nnodes, int* nodenums) :
  DGMHelm3d(_nnodes,nodenums) {
 ndir = 26;
}

DGMHelm3d_26t::DGMHelm3d_26t(int _nnodes, int* nodenums) :
  DGMHelm3d(-_nnodes,nodenums) {
 ndir = 26;
}

DGMHelm3d_56::DGMHelm3d_56(int _nnodes, int* nodenums) :
  DGMHelm3d(_nnodes,nodenums) {
 ndir = 56;
}

DGMHelm3d_56t::DGMHelm3d_56t(int _nnodes, int* nodenums) :
  DGMHelm3d(-_nnodes,nodenums) {
 ndir = 56;
}

DGMHelm3d_98::DGMHelm3d_98(int _nnodes, int* nodenums) :
  DGMHelm3d(_nnodes,nodenums) {
 ndir = 98;
}

void DGMHelm3d_6::dir(int n, complex<double> *d) {
 double kappa = prop ->kappaHelm;
 double a[][3] =  {
                    {1,0,0},
                    {-1,0,0},
                    {0,1,0},
                    {0,-1,0},
                    {0,0,1},
                    {0,0,-1}
                  };
 d[0] = complex<double>(0.0,kappa*a[n][0]);
 d[1] = complex<double>(0.0,kappa*a[n][1]);
 d[2] = complex<double>(0.0,kappa*a[n][2]);
}

void DGMHelm3d_6t::dir(int n, complex<double> *d) {
 double kappa = prop ->kappaHelm;
 double a[][3] =  {
                    {1,0,0},
                    {-1,0,0},
                    {0,1,0},
                    {0,-1,0},
                    {0,0,1},
                    {0,0,-1}
                  };
 d[0] = complex<double>(0.0,kappa*a[n][0]);
 d[1] = complex<double>(0.0,kappa*a[n][1]);
 d[2] = complex<double>(0.0,kappa*a[n][2]);
}


double* DGMHelm3d::getCubeDir(int n) {
 static double *a[3] = {0,0,0}; 
 if (a[n-2]==0) {
   a[n-2] = new double[(6*n*n+2)*3];
   int c = 0;
   for(int kk=0;kk<=n;kk++) for (int jj=0;jj<=n;jj++) for(int ii=0;ii<=n;ii++) {
     double cube[3];
     if (kk==0 || kk==n) {
       cube[0] = 0.5*tan(M_PI/4.0*(-1.0+2.0*ii/double(n)));
       cube[1] = 0.5*tan(M_PI/4.0*(-1.0+2.0*jj/double(n)));
       cube[2] = kk/double(n) -0.5;
     } else if (jj==0 || jj==n) {
       cube[0] = 0.5*tan(M_PI/4.0*(-1.0+2.0*ii/double(n)));
       cube[1] = jj/double(n) -0.5;
       cube[2] = 0.5*tan(M_PI/4.0*(-1.0+2.0*kk/double(n)));
     } else if (ii==0 || ii==n) {
       cube[0] = ii/double(n) -0.5;
       cube[1] = 0.5*tan(M_PI/4.0*(-1.0+2.0*jj/double(n)));
       cube[2] = 0.5*tan(M_PI/4.0*(-1.0+2.0*kk/double(n)));
     } else continue;
     double l = sqrt(cube[0]*cube[0]+cube[1]*cube[1]+cube[2]*cube[2]);
     a[n-2][c*3+0] = cube[0]/l;
     a[n-2][c*3+1] = cube[1]/l;
     a[n-2][c*3+2] = cube[2]/l;
     c++;
   }
 } 
 return a[n-2];
}

void DGMHelm3d_26::dir(int n, complex<double>* d) {
 double kappa = prop ->kappaHelm;
 double *a = getCubeDir(2);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_26t::dir(int n, complex<double>* d) {
 double kappa = prop ->kappaHelm;
 double *a = getCubeDir(2);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_56::dir(int n, complex<double>* d) {
 double kappa = prop ->kappaHelm;
 double *a = getCubeDir(3);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_56t::dir(int n, complex<double>* d) {
 double kappa = prop ->kappaHelm;
 double *a = getCubeDir(3);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DGMHelm3d_98::dir(int n, complex<double>* d) {
 double kappa = prop ->kappaHelm;
 double *a = getCubeDir(4);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}


void DGMHelm3d_1_LM::ldir(int n,double *tau1, double *tau2, complex<double>* d) {
 double lkappa = max(e1->getProperty()->kappaHelm,
                     e2->getProperty()->kappaHelm);
 d[0] = 0.0;
 d[1] = 0.0;
 d[2] = 0.0;
}

void DGMHelm3d_4_LM::ldir(int n,double *tau1, double *tau2, complex<double>* d) {
 double lkappa = max(e1->getProperty()->kappaHelm,
                     e2->getProperty()->kappaHelm);
 double a[][2] = { 
    { 1.0/sqrt(2.0), 1.0/sqrt(2.0) },
    { -1.0/sqrt(2.0), 1.0/sqrt(2.0) }, 
    { 1.0/sqrt(2.0), -1.0/sqrt(2.0) },
    { -1.0/sqrt(2.0), -1.0/sqrt(2.0) } };
 d[0] = complex<double>(0.0,0.6*lkappa*(a[n][0]*tau1[0]+a[n][1]*tau2[0]));
 d[1] = complex<double>(0.0,0.6*lkappa*(a[n][0]*tau1[1]+a[n][1]*tau2[1]));
 d[2] = complex<double>(0.0,0.6*lkappa*(a[n][0]*tau1[2]+a[n][1]*tau2[2]));
}

void DGMHelm3d_8_LM::ldir(int n,double *tau1, double *tau2, complex<double>* d) {
 double lkappa = max(e1->getProperty()->kappaHelm,
                     e2->getProperty()->kappaHelm);
 if (n<4) {
   double a[][2] = { 
      { 1.0/sqrt(2.0), 1.0/sqrt(2.0) },
      { -1.0/sqrt(2.0), 1.0/sqrt(2.0) }, 
      { 1.0/sqrt(2.0), -1.0/sqrt(2.0) },
      { -1.0/sqrt(2.0), -1.0/sqrt(2.0) } };
   d[0] = complex<double>(0.0,0.8*lkappa*(a[n][0]*tau1[0]+a[n][1]*tau2[0]));
   d[1] = complex<double>(0.0,0.8*lkappa*(a[n][0]*tau1[1]+a[n][1]*tau2[1]));
   d[2] = complex<double>(0.0,0.8*lkappa*(a[n][0]*tau1[2]+a[n][1]*tau2[2]));
 } else {
   n -= 4;
   double a[][2] = { 
      { 1,0} , { -1,0}, {0,1}, {0,-1} };
   d[0] = complex<double>(0.0,0.5*lkappa*(a[n][0]*tau1[0]+a[n][1]*tau2[0]));
   d[1] = complex<double>(0.0,0.5*lkappa*(a[n][0]*tau1[1]+a[n][1]*tau2[1]));
   d[2] = complex<double>(0.0,0.5*lkappa*(a[n][0]*tau1[2]+a[n][1]*tau2[2]));
 }
}

void DGMHelm3d_12_LM::ldir(int n,double *tau1, double *tau2, complex<double>* d) {
 double lkappa = max(e1->getProperty()->kappaHelm,
                     e2->getProperty()->kappaHelm);
 if (n<6) {
   double c = cos(n*M_PI/6.0+M_PI/12.0);
   double s = sin(n*M_PI/6.0+M_PI/12.0);
   d[0] = complex<double>(0.0,0.8*lkappa*(c*tau1[0]+s*tau2[0]));
   d[1] = complex<double>(0.0,0.8*lkappa*(c*tau1[1]+s*tau2[1]));
   d[2] = complex<double>(0.0,0.8*lkappa*(c*tau1[2]+s*tau2[2]));
 } else {
   double c = cos(n*M_PI/6.0);
   double s = sin(n*M_PI/6.0);
   d[0] = complex<double>(0.0,0.5*lkappa*(c*tau1[0]+s*tau2[0]));
   d[1] = complex<double>(0.0,0.5*lkappa*(c*tau1[1]+s*tau2[1]));
   d[2] = complex<double>(0.0,0.5*lkappa*(c*tau1[2]+s*tau2[2]));
 }
}


DGMHelm3d::DGMHelm3d(int _nnodes, int* nodenums) {
 if (_nnodes<0) {
   _nnodes=-_nnodes;
   if (_nnodes==4) o = -2;
   else if (_nnodes==10) o = -3;
   else if (_nnodes==20) o = -4;
   else if (_nnodes==35) o = -5;
   else if (_nnodes==56) o = -6;
   else {
     fprintf(stderr,"DGMHelm3d::DGMHelm3d: order too high\n"); exit(-1);
   }
 }
 else {
   o = int(pow(double(_nnodes),1.0/3.0));
 }
 nn = new int[_nnodes]; 
 for(int i=0;i<_nnodes;i++) nn[i] = nodenums[i];
 lm = new DEMLM*[nFaces()];
 for(int i=0;i<nFaces();i++) lm[i] = 0;
 bc = new int[nFaces()];
 for(int i=0;i<nFaces();i++) bc[i] = 0;
 ndir = 0;
}


void DGMHelm3d::getRef(double *xyz,double *cxyz) {
 if (o>0) {
   IsoParamUtils ipu(o);
   ipu.elementcenter(xyz,cxyz);  
 } else {
   IsoParamUtilsTetra ipu(-o);
   ipu.elementcenter(xyz,cxyz);  
 }
 cxyz[0] = cxyz[1] = cxyz[2] = 0.0;
}

void DGMHelm3d::createM(CoordSet &cs, complex<double>*M) {

 IsoParamUtils *ipu = (o>0)? new IsoParamUtils(o):new IsoParamUtilsTetra(-o);
 int os = ipu->getordersq();
 int oc = ipu->getorderc();
 double *xyz= new double[3*oc];
 cs.getCoordinates(nn,oc,xyz,xyz+oc,xyz+2*oc);

 double kappa = prop ->kappaHelm;
 double rho = prop ->rho;

 double xref[3];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) dir(i,cdir+i*3);

 int nldir[6] = {0,0,0,0,0,0};
 int nf,nv;
 int tnldir = 0;
 if (o>0) { nf=6; nv = 4; }
 else { nf = 4; nv = 3; }
 for(int fi=0;fi<nf;fi++) if (lm[fi]!=0) {
   nldir[fi] = lm[fi]->nDofs();
   tnldir += nldir[fi];
   DGMHelm3d_LM *l = dynamic_cast<DGMHelm3d_LM*>(lm[fi]);
//   DGMHelm3d_Eva_LM *le = dynamic_cast<DGMHelm3d_Eva_LM*>(lm[fi]);
//   if (l==0 && le==0) {
   if (l==0) {
     fprintf(stderr,"DGMHelm3d::createM: illegal LM type %d. Exiting.\n",
              lm[fi]->type());
   }
 }

// DGMHelm3d_LM type specific
 double xlref[18] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
// for(int fi=0;fi<nf;fi++) 
// ipu->sidecenter(xyz,fi+1,xlref+fi*3);

 complex<double>* kee = new complex<double>[ndir*ndir];
 complex<double>* kel = new complex<double>[ndir*tnldir];

 complex<double>* cldir = new complex<double>[3*tnldir];
 int corner[6] = { 1-1, o-1, os-1, os-o+1-1, 1-1, o-1};
 if (o<0) {
   corner[1] = -o-1;
   corner[3] = 1-1;
   corner[4] = -o-1;
 }
 int *fi = new int[os];
 int cc = 0;
 double sign[6] = { 1,1,1,1,1,1};
 for(int i=0;i<nf;i++) {
   ipu->faceindeces(i+1,fi);
   int cmin = 1;
   for(int j=2;j<=nv;j++) if ( nn[fi[corner[j]]]<nn[fi[corner[cmin]]]) cmin = j;
   double tau1[3], tau2[3];
   if (nn[fi[corner[cmin-1]]] < nn[fi[corner[cmin+1]]]) {
     tau1[0] = xyz[0*oc+fi[corner[cmin-1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau1[1] = xyz[1*oc+fi[corner[cmin-1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau1[2] = xyz[2*oc+fi[corner[cmin-1]]]-xyz[2*oc+fi[corner[cmin]]];
     tau2[0] = xyz[0*oc+fi[corner[cmin+1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau2[1] = xyz[1*oc+fi[corner[cmin+1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau2[2] = xyz[2*oc+fi[corner[cmin+1]]]-xyz[2*oc+fi[corner[cmin]]];
   } else {
     tau2[0] = xyz[0*oc+fi[corner[cmin-1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau2[1] = xyz[1*oc+fi[corner[cmin-1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau2[2] = xyz[2*oc+fi[corner[cmin-1]]]-xyz[2*oc+fi[corner[cmin]]];
     tau1[0] = xyz[0*oc+fi[corner[cmin+1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau1[1] = xyz[1*oc+fi[corner[cmin+1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau1[2] = xyz[2*oc+fi[corner[cmin+1]]]-xyz[2*oc+fi[corner[cmin]]];
   }
   double xsc[3];
   double xc[3];
   ipu->sidecenter(xyz,i+1,xsc);
   ipu->elementcenter(xyz,xc);
   xsc[0] -= xc[0];
   xsc[1] -= xc[1];
   xsc[2] -= xc[2];
   double cross[3] = {
        tau1[1]*tau2[2]-tau1[2]*tau2[1],
        tau1[2]*tau2[0]-tau1[0]*tau2[2],
        tau1[0]*tau2[1]-tau1[1]*tau2[0]};
   if (cross[0]*xsc[0]+cross[1]*xsc[1]+cross[2]*xsc[2]>0.0) sign[i] = -1.0;
   double l = sqrt(tau1[0]*tau1[0]+tau1[1]*tau1[1]+tau1[2]*tau1[2]);
   tau1[0] /= l;
   tau1[1] /= l;
   tau1[2] /= l;
   double p = tau1[0]*tau2[0]+tau1[1]*tau2[1]+tau1[2]*tau2[2];
   tau2[0] -= p*tau1[0];
   tau2[1] -= p*tau1[1];
   tau2[2] -= p*tau1[2];
   l = sqrt(tau2[0]*tau2[0]+tau2[1]*tau2[1]+tau2[2]*tau2[2]);
   tau2[0] /= l;
   tau2[1] /= l;
   tau2[2] /= l;
   if (nldir[i]>0) {
     DGMHelm3d_LM *l = dynamic_cast<DGMHelm3d_LM*>(lm[i]);
     for(int j=0;j<lm[i]->nDofs();j++) {
       l->ldir(j,tau1,tau2,cldir+3*cc);
       cc++;
     }
   }
 }
 delete[] fi;

 int arbFlag[6] = { 0,0,0,0,0,0};
 for(int i=0;i<nf;i++) if (bc[i]==1) arbFlag[i] = 1;

 HelmDGMEMatrices3d(o, xyz, ndir, cdir, nldir, cldir,
                    kappa, arbFlag, xlref,  xref, kee, kel); 

 cc = 0;
 for(int i=0;i<nf;i++) {
  for(int j=cc;j<cc+nldir[i];j++) for(int k=0;k<ndir;k++)
    kel[j*ndir+k] *= sign[i];
  cc += nldir[i];
 }
// for(int i=0;i<ndir*ndir;i++) fprintf(stderr,"K %d %d: %e %e\n",i/ndir,i%ndir,real(kee[i]),imag(kee[i]));
// for(int i=0;i<tnldir*ndir;i++) fprintf(stderr,"L %d %d: %e %e\n",i/ndir,i%ndir,real(kel[i]),imag(kel[i]));

 for(int i=0;i<ndir;i++) for(int j=0;j<ndir;j++)
   M[(tnldir+j)*(ndir+tnldir)+(tnldir+i)] = kee[j*ndir+i]/rho;
 for(int i=0;i<ndir;i++) for(int j=0;j<tnldir;j++) {
   M[(tnldir+i)*(ndir+tnldir)+j] = kel[j*ndir+i]/rho;
   M[j*(ndir+tnldir)+(tnldir+i)] = kel[j*ndir+i]/rho;
 }
 for(int i=0;i<tnldir;i++) for(int j=0;j<tnldir;j++) 
   M[j*(ndir+tnldir)+i] = 0.0;

 delete[] kee;
 delete[] kel;
 delete[] xyz;
 delete[] cldir;
 delete[] cdir;
 delete ipu;
}

class HelmDGMENeumVFunction3d : public IntegFunctionA3d {
 int ndir;
 complex<double> *dirs;
 complex<double> *incdir;
 double *xc;
 complex<double> *v;
public:
 HelmDGMENeumVFunction3d(
                      int _ndir, complex<double> *_dirs,
                      complex<double> *_incdir, double *_xc,
                      complex<double> *_v) {
   ndir = _ndir; dirs = _dirs; incdir = _incdir; xc = _xc;
   v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   int i;
   complex<double> *e = new complex<double>[ndir];
   for(i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                  dirs[i*3+1]*(x[1]-xc[1]) +
                                  dirs[i*3+2]*(x[2]-xc[2]));

   complex<double> ince = w* (nsign*
        (incdir[0]*cross[0]+incdir[1]*cross[1]+incdir[2]*cross[2])) *
        exp(incdir[0]*x[0] + incdir[1]*x[1] + incdir[2]*x[2]);

   for(i=0;i<ndir;i++)
     v[i] += ince*e[i];
   delete[] e; 
 }
};


void HelmDGMENeumV(int order, double *xyz,
                     int ndir, complex<double> *dirs,
                     double kappa, complex<double> *incdir, int faceindex,
                     double *xc,
                     complex<double>* v) {
 IsoParamUtils *ipu = (order>0)? new IsoParamUtils(order):new IsoParamUtilsTetra(-order);
 int gorder = (order>0)?7:13;
// int os = ipu->getordersq();
// double *fxyz = new double[os*3];
// int *fi = new int[os];

// ipu.faceindeces(faceindex,fi);
// for(int i=0;i<os;i++) {
//   fxyz[i*3+0] = xyz[fi[i]*3+0];
//   fxyz[i*3+1] = xyz[fi[i]*3+1];
//   fxyz[i*3+2] = xyz[fi[i]*3+2];
// }
 HelmDGMENeumVFunction3d f(ndir,dirs,incdir,xc,v);
// ipu.surfSurfInt3d(fxyz, f);
 ipu->surfInt3d(xyz,faceindex, f,gorder);

// delete[] fxyz;
// delete[] fi;
 delete ipu;
}

class HelmDGMERobVFunction3d : public IntegFunctionA3d {
 int ndir;
 complex<double> *dirs;
 complex<double> *incdir;
 double *xc;
 complex<double> *v;
 double kappa;
public:
 HelmDGMERobVFunction3d(
                      int _ndir, complex<double> *_dirs,
                      complex<double> *_incdir, double *_xc, double _kappa,
                      complex<double> *_v) {
   ndir = _ndir; dirs = _dirs; incdir = _incdir; xc = _xc; kappa = _kappa;
   v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {

   int i;
   complex<double> *e = new complex<double>[ndir];
   for(i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                  dirs[i*3+1]*(x[1]-xc[1]) +
                                  dirs[i*3+2]*(x[2]-xc[2]));

   complex<double> ince = w* ((nsign*
        (incdir[0]*cross[0]+incdir[1]*cross[1]+incdir[2]*cross[2])) -
complex<double>(0.0,kappa*
 sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2])) ) *
        exp(incdir[0]*x[0] + incdir[1]*x[1] + incdir[2]*x[2]);

   for(i=0;i<ndir;i++)
     v[i] += ince*e[i];
   delete[] e; 
 }
};


void HelmDGMERobV(int order, double *xyz,
                     int ndir, complex<double> *dirs,
                     double kappa, complex<double> *incdir, int faceindex,
                     double *xc,
                     complex<double>* v) {
 IsoParamUtils *ipu = (order>0)? new IsoParamUtils(order):new IsoParamUtilsTetra(-order);
 int gorder = (order>0)?7:13;

 HelmDGMERobVFunction3d f(ndir,dirs,incdir,xc,kappa,v);
 ipu->surfInt3d(xyz,faceindex, f, gorder);
 delete ipu;
}



void DGMHelm3d::createRHS(CoordSet &cs, complex<double>*v) {
 int oc;
 if (o>0) {
   IsoParamUtils ipu(o);
   oc = ipu.getorderc();
 } else {
   IsoParamUtilsTetra ipu(-o);
   oc = ipu.getorderc();
 }
 double *xyz= new double[3*oc];
 cs.getCoordinates(nn,oc,xyz,xyz+oc,xyz+2*oc);

 double kappa = prop ->kappaHelm;
 double rho = prop ->rho;

 double xref[3];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) dir(i,cdir+i*3);

 complex<double> incdir[3] = {complex<double>(0.0,kappa*0.0),
                              complex<double>(0.0,kappa*0.0),
                              complex<double>(0.0,kappa*1.0)};
 incdir[0] = complex<double>(0.0,kappa*2.672612419124244e-01);
 incdir[1] = complex<double>(0.0,kappa*5.345224838248488e-01);
 incdir[2] = complex<double>(0.0,kappa*8.017837257372732e-01);

 complex<double> *vv = new complex<double>[ndir];
 for(int i=0;i<ndir;i++) vv[i] = 0.0;

 for(int i=0;i<nFaces();i++) {
//   if (bc[i]==2)
//     HelmDGMENeumV(o, xyz, ndir, cdir, kappa, incdir, i+1 , xref, vv);
   if (bc[i]==1)
     HelmDGMERobV(o, xyz, ndir, cdir, kappa, incdir, i+1 , xref, vv);
 }

 delete[] cdir;
 delete[] xyz;

 int tnldir = nLagrangeDofs();
 for(int i=0;i<tnldir;i++) v[i] = 0;
 for(int i=0;i<ndir;i++) 
   v[tnldir+i] = vv[i]/rho;
 delete[] vv;
}


void DGMHelm3d::createSol(double *xyz,
                        complex<double>* sol, complex<double> *nodalSol) {

 for(int i=0;i<8;i++) nodalSol[i] = 0.0;
 for(int i=0;i<ndir;i++) {
   complex<double> d[3];
   dir(i,d);
   nodalSol[0] += sol[i]*exp( d[0]*xyz[0]+d[1]*xyz[1]+d[2]*xyz[2] );
 }
}


DEMHelm3d::DEMHelm3d(int _nnodes, int* nodenums) : DGMHelm3d(_nnodes, nodenums) {}


DEMHelm3d_6::DEMHelm3d_6(int _nnodes, int* nodenums) :
  DEMHelm3d(_nnodes,nodenums) {
 ndir = 6;
}

DEMHelm3d_26::DEMHelm3d_26(int _nnodes, int* nodenums) :
  DEMHelm3d(_nnodes,nodenums) {
 ndir = 26;
}

DEMHelm3d_56::DEMHelm3d_56(int _nnodes, int* nodenums) :
  DEMHelm3d(_nnodes,nodenums) {
 ndir = 56;
}

DEMHelm3d_98::DEMHelm3d_98(int _nnodes, int* nodenums) :
  DEMHelm3d(_nnodes,nodenums) {
 ndir = 98;
}

void DEMHelm3d_6::dir(int n, complex<double> *d) {
 double kappa = prop ->kappaHelm;
 double a[][3] =  {
                    {1,0,0},
                    {-1,0,0},
                    {0,1,0},
                    {0,-1,0},
                    {0,0,1},
                    {0,0,-1}
                  };
 d[0] = complex<double>(0.0,kappa*a[n][0]);
 d[1] = complex<double>(0.0,kappa*a[n][1]);
 d[2] = complex<double>(0.0,kappa*a[n][2]);
}

void DEMHelm3d_26::dir(int n, complex<double>* d) {
 double kappa = prop ->kappaHelm;
 double *a = getCubeDir(2);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DEMHelm3d_56::dir(int n, complex<double>* d) {
 double kappa = prop ->kappaHelm;
 double *a = getCubeDir(3);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}

void DEMHelm3d_98::dir(int n, complex<double>* d) {
 double kappa = prop ->kappaHelm;
 double *a = getCubeDir(4);
 d[0] = complex<double>(0.0,kappa*a[n*3+0]);
 d[1] = complex<double>(0.0,kappa*a[n*3+1]);
 d[2] = complex<double>(0.0,kappa*a[n*3+2]);
}


class HelmPMatricesFunction3d : public IntegFunctionV3d {
 int o;
 int ndir;
 complex<double> *dirs;
 double *xc;
 double kappa;
 complex<double> *PP;
 complex<double> *PE;
public:
 HelmPMatricesFunction3d(int _o,  double _kappa,
                         int _ndir, complex<double> *_dirs, double *_xc,
                      complex<double> *_PE, complex<double> *_PP) {
   o = _o; kappa = _kappa; ndir = _ndir; dirs = _dirs; xc = _xc;
   PP = _PP; PE = _PE;
 }
 void evaluate(double *x, double *N, double (*dNdx)[3], double w, double det) {
    
   int oc = o*o*o;

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1])+
                                      dirs[i*3+2]*(x[2]-xc[2]));
   for(int j=0;j<ndir;j++)
     for(int i=0;i<oc;i++)
       PE[j*oc+i] += w*det*e[j]*(-kappa*kappa*N[i] +
     dNdx[i][0]*dirs[j*3+0]+dNdx[i][1]*dirs[j*3+1]+dNdx[i][2]*dirs[j*3+2]);
   delete[] e;
   for(int j=0;j<oc;j++)
     for(int i=j;i<oc;i++)
       PP[j*oc+i] += w*det*(-kappa*kappa*N[i]*N[j]
        +dNdx[i][0]*dNdx[j][0]+dNdx[i][1]*dNdx[j][1]+dNdx[i][2]*dNdx[j][2]);
 }
};


void HelmDEMMatrices3d(int order, double *xyz,
                    int ndir, complex<double> *dirs,
                    int *nldirs, complex<double> *ldirs,
                    double kappa, int *sflags, double *xsc,
                    double *xc,
                    complex<double> *kee, complex<double> *kel,
                    complex<double> *kpp, complex<double> *kpl,
                    complex<double> *kpe) {

 IsoParamUtils ipu(order);
 int oc = ipu.getorderc();

 int nldir =nldirs[0] + nldirs[1] + nldirs[2] +
            nldirs[3] + nldirs[4] + nldirs[5];
 for(int i=0;i<nldir*ndir;i++) kel[i] = 0.0;
 for(int i=0;i<ndir*ndir;i++) kee[i] = 0.0;
 for(int i=0;i<oc*ndir;i++) kpe[i] = 0.0;
 for(int i=0;i<oc*nldir;i++) kpl[i] = 0.0;
 for(int i=0;i<oc*oc;i++) kpp[i] = 0.0;

 int c = 0;
 for(int faceindex=1;faceindex<=6;faceindex++) {
   HelmDGMEMatricesFunction3d f(order,kappa,ndir,dirs,
                nldirs[faceindex-1],ldirs+c*3,
                xsc + (faceindex-1)*3, xc, sflags[faceindex-1], faceindex,
                kee,kel+c*ndir,kpl+c*ndir,kpp,kpe);
   ipu.surfInt3d(xyz, faceindex, f);
   c += nldirs[faceindex-1];
 }
 for(int j=0;j<ndir;j++) {
   for(int i=0;i<j;i++)
     kee[j*ndir+i] = kee[i*ndir+j];
 }

 HelmPMatricesFunction3d f(order,kappa,ndir,dirs, xc, kpe,kpp);
 ipu.volumeInt3d(xyz, f);

 for(int j=0;j<oc;j++) {
   for(int i=0;i<j;i++)
     kpp[j*oc+i] = kpp[i*oc+j];
 }
}


void DEMHelm3d::createM(CoordSet &cs, complex<double>*M) {

 IsoParamUtils ipu(o);
 int os = ipu.getordersq();
 int oc = ipu.getorderc();
 double *xyz= new double[3*oc];
 cs.getCoordinates(nn,oc,xyz,xyz+oc,xyz+2*oc);

 double kappa = prop ->kappaHelm;
 double rho = prop ->rho;

 double xref[3];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) dir(i,cdir+i*3);

 int nldir[6] = {0,0,0,0,0,0};
 for(int fi=0;fi<6;fi++) if (lm[fi]!=0) nldir[fi] = lm[fi]->nDofs();

 int tnldir = 0;
 for(int i=0;i<6;i++) tnldir += nldir[i];

 for(int fi=0;fi<6;fi++) if (lm[fi]!=0) {
   DGMHelm3d_LM *l = dynamic_cast<DGMHelm3d_LM*>(lm[fi]);
   if (l==0) {
     fprintf(stderr,"DGMHelm3d::createM: illegal LM type %d. Exiting.\n",
              lm[fi]->type());
   }
 }

// DGMHelm3d_LM type specific
 double xlref[18] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };
// for(int fi=0;fi<6;fi++) 
// ipu.sidecenter(xyz,fi+1,xlref+fi*3);

 complex<double>* cldir = new complex<double>[3*tnldir];
 int corner[6] = { 1-1, o-1, os-1, os-o+1-1, 1-1, o-1};
 int *fi = new int[os];
 int cc = 0;
 double sign[6] = { 1,1,1,1,1,1};
 for(int i=0;i<6;i++) {
   ipu.faceindeces(i+1,fi);
   int cmin = 1;
   for(int j=2;j<5;j++) if ( nn[fi[corner[j]]]<nn[fi[corner[cmin]]]) cmin = j;
   double tau1[3], tau2[3];
   if (nn[fi[corner[cmin-1]]] < nn[fi[corner[cmin+1]]]) {
     tau1[0] = xyz[0*oc+fi[corner[cmin-1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau1[1] = xyz[1*oc+fi[corner[cmin-1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau1[2] = xyz[2*oc+fi[corner[cmin-1]]]-xyz[2*oc+fi[corner[cmin]]];
     tau2[0] = xyz[0*oc+fi[corner[cmin+1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau2[1] = xyz[1*oc+fi[corner[cmin+1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau2[2] = xyz[2*oc+fi[corner[cmin+1]]]-xyz[2*oc+fi[corner[cmin]]];
   } else {
     tau2[0] = xyz[0*oc+fi[corner[cmin-1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau2[1] = xyz[1*oc+fi[corner[cmin-1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau2[2] = xyz[2*oc+fi[corner[cmin-1]]]-xyz[2*oc+fi[corner[cmin]]];
     tau1[0] = xyz[0*oc+fi[corner[cmin+1]]]-xyz[0*oc+fi[corner[cmin]]];
     tau1[1] = xyz[1*oc+fi[corner[cmin+1]]]-xyz[1*oc+fi[corner[cmin]]];
     tau1[2] = xyz[2*oc+fi[corner[cmin+1]]]-xyz[2*oc+fi[corner[cmin]]];
   }
   double xsc[3];
   double xc[3];
   ipu.sidecenter(xyz,i+1,xsc);
   ipu.elementcenter(xyz,xc);
   xsc[0] -= xc[0];
   xsc[1] -= xc[1];
   xsc[2] -= xc[2];
   double cross[3] = {
        tau1[1]*tau2[2]-tau1[2]*tau2[1],
        tau1[2]*tau2[0]-tau1[0]*tau2[2],
        tau1[0]*tau2[1]-tau1[1]*tau2[0]};
   if (cross[0]*xsc[0]+cross[1]*xsc[1]+cross[2]*xsc[2]>0.0) sign[i] = -1.0;
   double l = sqrt(tau1[0]*tau1[0]+tau1[1]*tau1[1]+tau1[2]*tau1[2]);
   tau1[0] /= l;
   tau1[1] /= l;
   tau1[2] /= l;
   double p = tau1[0]*tau2[0]+tau1[1]*tau2[1]+tau1[2]*tau2[2];
   tau2[0] -= p*tau1[0];
   tau2[1] -= p*tau1[1];
   tau2[2] -= p*tau1[2];
   l = sqrt(tau2[0]*tau2[0]+tau2[1]*tau2[1]+tau2[2]*tau2[2]);
   tau2[0] /= l;
   tau2[1] /= l;
   tau2[2] /= l;
   if (nldir[i]>0) {
     DGMHelm3d_LM *l = dynamic_cast<DGMHelm3d_LM*>(lm[i]);
     for(int j=0;j<lm[i]->nDofs();j++) {
       l->ldir(j,tau1,tau2,cldir+3*cc);
       cc++;
     }
   }
 }
 delete[] fi;

 int arbFlag[6] = { 0,0,0,0,0,0};
 for(int i=0;i<6;i++) if (bc[i]==1) arbFlag[i] = 1;

 complex<double>* kee = new complex<double>[ndir*ndir];
 complex<double>* kel = new complex<double>[ndir*tnldir];
 complex<double>* kpp = new complex<double>[oc*oc];
 complex<double>* kpe = new complex<double>[ndir*oc];
 complex<double>* kpl = new complex<double>[oc*tnldir];
 HelmDEMMatrices3d(o, xyz, ndir, cdir, nldir, cldir,
                    kappa, arbFlag, xlref,  xref, kee, kel, kpp, kpl, kpe); 

 cc = 0;
 for(int i=0;i<6;i++) {
  for(int j=cc;j<cc+nldir[i];j++) {
    for(int k=0;k<ndir;k++) kel[j*ndir+k] *= sign[i];
    for(int k=0;k<oc;k++) kpl[j*oc+k] *= sign[i];
  }
  cc += nldir[i];
 }

 for(int i=0;i<oc;i++) for(int j=0;j<oc;j++) M[j*(oc+tnldir+ndir)+i] = kpp[j*oc+i]/rho;
 for(int i=0;i<ndir;i++) for(int j=0;j<ndir;j++)
   M[(oc+tnldir+j)*(oc+ndir+tnldir)+(oc+tnldir+i)] = kee[j*ndir+i]/rho;
 for(int i=0;i<ndir;i++) for(int j=0;j<tnldir;j++) {
   M[(oc+tnldir+i)*(oc+ndir+tnldir)+oc+j] = kel[j*ndir+i]/rho;
   M[(oc+j)*(oc+ndir+tnldir)+(oc+tnldir+i)] = kel[j*ndir+i]/rho;
 }
 for(int i=0;i<oc;i++) for(int j=0;j<ndir;j++) {
   M[i*(oc+ndir+tnldir)+(oc+tnldir+j)] = kpe[j*oc+i]/rho;
   M[(oc+tnldir+j)*(oc+ndir+tnldir)+i] = kpe[j*oc+i]/rho;
 }
 for(int i=0;i<oc;i++) for(int j=0;j<tnldir;j++) {
   M[i*(oc+ndir+tnldir)+(oc+j)] = kpl[j*oc+i]/rho;
   M[(oc+j)*(oc+ndir+tnldir)+i] = kpl[j*oc+i]/rho;
 }
 for(int i=0;i<tnldir;i++) for(int j=0;j<tnldir;j++)
   M[(oc+j)*(oc+ndir+tnldir)+oc+i] = 0.0;

 delete[] kee;
 delete[] kel;
 delete[] kpp;
 delete[] kpe;
 delete[] kpl;
 delete[] xyz;
 delete[] cldir;
 delete[] cdir;
}

class HelmDEMNeumVFunction : public IntegFunctionA3d {
 int order;
 int ndir;
 complex<double> *dirs;
 complex<double> *incdir;
 double *xc;
 complex<double> *v;
public:
 HelmDEMNeumVFunction(int _order,
                      int _ndir, complex<double> *_dirs,
                      complex<double> *_incdir, double *_xc,
                      complex<double> *_v) {
   order = _order; ndir = _ndir; dirs = _dirs; incdir = _incdir; xc = _xc;
   v = _v;
 }
 void evaluate(double *x, double *N, double *cross, double nsign, double w) {
  
   int oc = order*order*order;
   
   complex<double> ince = w* (nsign* 
                  (incdir[0]*cross[0]+incdir[1]*cross[1]+incdir[2]*cross[2])) *
                  exp(incdir[0]*x[0] + incdir[1]*x[1] + incdir[2]*x[2]);

   for(int i=0;i<oc;i++) v[i] += ince*N[i];

   complex<double> *e = new complex<double>[ndir];
   for(int i=0;i<ndir;i++) e[i] = exp(dirs[i*3+0]*(x[0]-xc[0]) +
                                      dirs[i*3+1]*(x[1]-xc[1])+
                                      dirs[i*3+2]*(x[2]-xc[2]));
   for(int i=0;i<ndir;i++)
     v[oc+i] += ince*e[i];
   delete[] e;
 }
};


void HelmDEMNeumV(int order, double *xy,
                     int ndir, complex<double> *dirs,
                     complex<double> *incdir, int faceindex,
                     double *xc,
                     complex<double>* v) {
 IsoParamUtils ipu(order);

 HelmDEMNeumVFunction f(order,ndir,dirs,incdir,xc,v);
 ipu.surfInt3d(xy, faceindex, f);
}



void DEMHelm3d::createRHS(CoordSet &cs, complex<double>*v) {
 IsoParamUtils ipu(o);
 int oc = ipu.getorderc();
 double *xyz= new double[3*oc];
 cs.getCoordinates(nn,oc,xyz,xyz+oc,xyz+2*oc);

 double kappa = prop ->kappaHelm;
 double rho = prop ->rho;

 double xref[3];
 getRef(xyz,xref);

 complex<double>* cdir = new complex<double>[ndir*3];
 for(int i=0;i<ndir;i++) dir(i,cdir+i*3);

 complex<double> incdir[3] = {complex<double>(0.0,kappa*0.0),
                              complex<double>(0.0,kappa*0.0),
                              complex<double>(0.0,kappa*1.0)};

 complex<double> *vv = new complex<double>[nPolynomialDofs()+ndir];
 for(int i=0;i<nPolynomialDofs()+ndir;i++) vv[i] = 0.0;

 for(int i=0;i<nFaces();i++) {
   if (bc[i]==2)
     HelmDEMNeumV(o, xyz, ndir, cdir, incdir, i+1 , xref, vv);
 }

 delete[] cdir;
 delete[] xyz;

 for(int i=0;i<nPolynomialDofs();i++)
   v[i] = vv[i]/rho;
 int tnldir = nLagrangeDofs();
 for(int i=0;i<tnldir;i++) v[nPolynomialDofs()+i] = 0;
 for(int i=0;i<ndir;i++)
   v[nPolynomialDofs()+tnldir+i] = vv[nPolynomialDofs()+i]/rho;
 delete[] vv;
}
