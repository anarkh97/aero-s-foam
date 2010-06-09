#ifndef _DEMHELM3D_H_
#define _DEMHELM3D_H_

#include <math.h>
#ifndef SANDIA
#include <Element.d/DEM.d/DEMElement.h>
#include <Element.d/Helm.d/IntegFunction.h>
#else
#include "DEMElement.h"
#endif

class DGMHelm3d_LM: public DEMLM {
public:
 virtual int nDofs()=0;
 virtual int type()=0;
 virtual void init() {};
 virtual void ldir(int,double*,double*,complex<double>*)=0;
};


class DGMHelm3d_1_LM: public DGMHelm3d_LM {
public:
 virtual int nDofs() { return 1; }
 virtual int type() { return 51; }
 virtual void ldir(int,double*,double*,complex<double>*);
};

class DGMHelm3d_4_LM: public DGMHelm3d_LM {
public:
 virtual int nDofs() { return 4; }
 virtual int type() { return 52; }
 virtual void ldir(int,double*,double*,complex<double>*);
};

class DGMHelm3d_8_LM: public DGMHelm3d_LM {
public:
 virtual int nDofs() { return 8; }
 virtual int type() { return 53; }
 virtual void ldir(int,double*,double*,complex<double>*);
};

class DGMHelm3d_12_LM: public DGMHelm3d_LM {
public:
 virtual int nDofs() { return 12; }
 virtual int type() { return 54; }
 virtual void ldir(int,double*,double*,complex<double>*);
};




class DGMHelm3d: public DEMElement {
protected:
public:
 void HelmDGMEMatrices3d(double *xyz,
                    int ndir, complex<double> *dirs,
                    int *nldirs, complex<double> *ldirs,
                    double kappa, int *sflags, double *xsc,
                    double *xc,
                    complex<double> *K, complex<double> *L);
 void HelmDGMPMLEMatrices3d(double *xyz,
                    int ndir, complex<double> *dirs,
                    int *nldirs, complex<double> *ldirs,
                    double kappa, int *sflags, double *xsc,
                    double *xc, PMLProps *pml,
                    complex<double> *K, complex<double> *L);
 void HelmDGMENeumV(double *xyz,
                     int ndir, complex<double> *dirs,
                     double kappa, complex<double> *incdir, int faceindex,
                     double *xc,
                     complex<double>* v);
 double * getCubeDir(int n);
 int ndir;
 virtual void dir(int,complex<double>*) {};
 int o;
 void init(int _nnodes, int* nodenums);
 virtual int defaultLMType() { return 51; }
 virtual bool dgmFlag() { return true; }
 virtual bool condensedFlag() {
   for(int i=0;i<nFaces();i++) if (bc[i]==3) return false&&condensedF;
   return true&&condensedF;
 }
 virtual bool storeMatricesFlag() { return storeMatricesF; }
// virtual bool storeMatricesFlag() { return false; }
 virtual int nPolynomialDofs() { return 0; }
 virtual int nEnrichmentDofs() { return ndir; }
 virtual int *faceCornerI(int fi)=0; 
 virtual int *faceCorners(int fi) {
  int *ix = faceCornerI(fi);
  int nf = nFaceCorners(fi);
  for(int i=0;i<nf;i++) {
    int tmp = ix[i]; ix[i] = nn[tmp];
  }
  return ix;
 }

 virtual void externalNormal(double *xyz, int faceindex, double *n)=0;
 virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f)=0;
 virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f)=0;
 virtual void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f)=0;

#ifndef SANDIA
 virtual int polyDofType() { return DofSet::Helm; }
#endif
 virtual int polyDofsPerNode() { return 1; }

 virtual void getRef(double *xyz,double *xy);
 virtual void createM(complex<double>*);
 virtual void interfMatrix(int fi, DEMElement*,complex<double>*);
 virtual void createRHS(complex<double>*);
 virtual void createSol(double *xyz, complex<double>*,
                            complex<double>*);
};

class HexDGMElement3d: public DGMHelm3d {
public:
 HexDGMElement3d(int _n, int* nodenums);
 virtual int nGeomNodes() { return o*o*o; }
 virtual int nFaces() { return 6; }
 virtual int nFaceCorners(int fi) { return 4; }
 virtual int *faceCornerI(int fi);
 virtual void externalNormal(double *xyz, int faceindex, double *n);
 virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f);
 virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f);
 virtual void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f);
};


class TetraDGMElement3d: public DGMHelm3d {
public:
 TetraDGMElement3d(int _n, int* nodenums);
 virtual int nGeomNodes() { return (o*(o+1)*(o+2))/6; }
 virtual int nFaces() { return 4; }
 virtual int nFaceCorners(int fi) { return 3; }
 virtual int *faceCornerI(int fi);
 virtual void externalNormal(double *xyz, int faceindex, double *n);
 virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f);
 virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f);
 virtual void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f);
};


class PrismDGMElement3d: public DGMHelm3d {
public:
 PrismDGMElement3d(int _n, int* nodenums);
 virtual int nGeomNodes() { return o*(o*(o+1))/2; }
 virtual int nFaces() { return 5; }
 virtual int nFaceCorners(int fi) { return (fi<3)?3:4; }
 virtual int *faceCornerI(int fi);
 virtual void externalNormal(double *xyz, int faceindex, double *n);
 virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f);
 virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f);
 virtual void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f);
};


class PyramidDGMElement3d: public DGMHelm3d {
public:
 PyramidDGMElement3d(int _n, int* nodenums);
 virtual int nGeomNodes() { return 5; }
 virtual int nFaces() { return 5; }
 virtual int nFaceCorners(int fi) { return (fi<2)?4:3; }
 virtual int *faceCornerI(int fi);
 virtual void externalNormal(double *xyz, int faceindex, double *n);
 virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f);
 virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f);
 virtual void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f);
};


class DGMHelm3d_6: public HexDGMElement3d {
public:
 DGMHelm3d_6(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 51; }
};


class DGMHelm3d_6t: public TetraDGMElement3d {
public:
 DGMHelm3d_6t(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 51; }
};


class DGMHelm3d_6p: public PrismDGMElement3d {
public:
 DGMHelm3d_6p(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 51; }
};

class DGMHelm3d_6pd: public PyramidDGMElement3d {
public:
 DGMHelm3d_6pd(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 51; }
};

class DGMHelm3d_26: public HexDGMElement3d {
public:
 DGMHelm3d_26(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 52; }
};

class DGMHelm3d_26t: public TetraDGMElement3d {
public:
 DGMHelm3d_26t(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 52; }
};

class DGMHelm3d_26p: public PrismDGMElement3d {
public:
 DGMHelm3d_26p(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 52; }
};

class DGMHelm3d_26pd: public PyramidDGMElement3d {
public:
 DGMHelm3d_26pd(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 52; }
};

class DGMHelm3d_56: public HexDGMElement3d {
public:
 DGMHelm3d_56(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 53; }
};

class DGMHelm3d_56t: public TetraDGMElement3d {
public:
 DGMHelm3d_56t(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 53; }
};

class DGMHelm3d_56p: public PrismDGMElement3d {
public:
 DGMHelm3d_56p(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 53; }
};

class DGMHelm3d_56pd: public PyramidDGMElement3d {
public:
 DGMHelm3d_56pd(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 53; }
};

class DGMHelm3d_98: public HexDGMElement3d {
public:
 DGMHelm3d_98(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 54; }
};


class DEMHelm3d: public HexDGMElement3d {
public:
 DEMHelm3d(int _o, int* nodenums);
 virtual int defaultLMType() { return 51; }
 virtual bool dgmFlag() { return false; }
 virtual int nPolynomialDofs() { return o*o*o; }

 virtual void createM(complex<double>*);
 virtual void interfMatrix(int fi, DEMElement*,complex<double>*);
 virtual void createRHS(complex<double>*);
};


class DEMHelm3d_6: public DEMHelm3d {
public:
 DEMHelm3d_6(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 51; }
};

class DEMHelm3d_26: public DEMHelm3d {
public:
 DEMHelm3d_26(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 52; }
};

class DEMHelm3d_56: public DEMHelm3d {
public:
 DEMHelm3d_56(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 53; }
};

class DEMHelm3d_98: public DEMHelm3d {
public:
 DEMHelm3d_98(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 54; }
};

#endif
