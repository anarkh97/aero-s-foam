#ifndef _DEMHELM3D_H_
#define _DEMHELM3D_H_

#include <cmath>
#ifndef SANDIA
#include <Element.d/DEM.d/DEMElement.h>
#include <Element.d/Helm.d/IntegFunction.h>
#else
#include "DEMElement.h"
#endif

class DGMHelm3d_LM: public DEMLM {
public:
    virtual void init() {};
    virtual void ldir(int,double*,double*,complex<double>*)=0;
};


class DGMHelm3d_1_LM: public DGMHelm3d_LM {
public:
    int nDofs() const override { return 1; }
    int type() const override { return 51; }
    void ldir(int,double*,double*,complex<double>*) override;
};

class DGMHelm3d_4_LM: public DGMHelm3d_LM {
public:
    int nDofs() const override { return 4; }
    int type() const override { return 52; }
    void ldir(int,double*,double*,complex<double>*) override;
};

class DGMHelm3d_8_LM: public DGMHelm3d_LM {
public:
    int nDofs() const override { return 8; }
    int type() const override { return 53; }
    void ldir(int,double*,double*,complex<double>*) override;
};

class DGMHelm3d_12_LM: public DGMHelm3d_LM {
public:
    int nDofs() const override { return 12; }
    int type() const override { return 54; }
    void ldir(int,double*,double*,complex<double>*) override;
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
    void HelmDGMEMatricesExactFace3d(double *xyz,
                                     int ndir, complex<double> *dirs,
                                     int nldir, complex<double> *ldirs,
                                     double kappa,  int sflag,
                                     double *xsc,  double *xc, int faceindex,
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
    int defaultLMType() const override { return 51; }
    bool dgmFlag() const override { return true; }
    bool condensedFlag() const override  {
        for(int i=0;i<nFaces();i++) if (bc[i]==3) return false&&condensedF;
        return true&&condensedF;
    }
    virtual bool storeMatricesFlag() { return storeMatricesF; }
// virtual bool storeMatricesFlag() { return false; }
    int nPolynomialDofs() const  { return 0; }
    virtual int nEnrichmentDofs() const override { return ndir; }
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
    virtual int isFlatAndStraight(double *xyz,int faceindex)=0;

#ifndef SANDIA
    virtual int polyDofType() const override { return DofSet::Helm; }
#endif
    virtual int polyDofsPerNode() const override { return 1; }

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
    int nGeomNodes() const override { return o*o*o; }
    int nFaces() const override { return 6; }
    int nFaceCorners(int fi) const override { return 4; }
    virtual int *faceCornerI(int fi);
    virtual void externalNormal(double *xyz, int faceindex, double *n);
    virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f);
    virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f);
    virtual void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f);
    virtual int isFlatAndStraight(double *xyz,int faceindex);
};


class TetraDGMElement3d: public DGMHelm3d {
public:
    TetraDGMElement3d(int _n, int* nodenums);
    int nGeomNodes() const override { return (o*(o+1)*(o+2))/6; }
    int nFaces() const override { return 4; }
    int nFaceCorners(int fi) const override { return 3; }
    virtual int *faceCornerI(int fi);
    virtual void externalNormal(double *xyz, int faceindex, double *n);
    virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f);
    virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f);
    virtual void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f);
    virtual int isFlatAndStraight(double *xyz,int faceindex);
};


class PrismDGMElement3d: public DGMHelm3d {
public:
    PrismDGMElement3d(int _n, int* nodenums);
    int nGeomNodes() const override { return o*(o*(o+1))/2; }
    int nFaces() const override { return 5; }
    int nFaceCorners(int fi) const override { return (fi<3)?3:4; }
    virtual int *faceCornerI(int fi);
    virtual void externalNormal(double *xyz, int faceindex, double *n);
    virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f);
    virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f);
    virtual void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f);
    virtual int isFlatAndStraight(double *xyz,int faceindex);
};


class PyramidDGMElement3d: public DGMHelm3d {
public:
    PyramidDGMElement3d(int _n, int* nodenums);
    int nGeomNodes() const override { return 5; }
    int nFaces() const override { return 5; }
    int nFaceCorners(int fi) const override { return (fi<2)?4:3; }
    virtual int *faceCornerI(int fi);
    virtual void externalNormal(double *xyz, int faceindex, double *n);
    virtual void volumeInt3d(double *xyz, IntegFunctionV3d &f);
    virtual void surfInt3d(double *xyz, int faceindex, IntegFunctionA3d &f);
    virtual void surftInt3d(double *xyz, int faceindex, IntegFunctionAt3d &f);
    virtual int isFlatAndStraight(double *xyz,int faceindex);
};


class DGMHelm3d_6: public HexDGMElement3d {
public:
    DGMHelm3d_6(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 51; }
};


class DGMHelm3d_6t: public TetraDGMElement3d {
public:
    DGMHelm3d_6t(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 51; }
};


class DGMHelm3d_6p: public PrismDGMElement3d {
public:
    DGMHelm3d_6p(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 51; }
};

class DGMHelm3d_6pd: public PyramidDGMElement3d {
public:
    DGMHelm3d_6pd(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 51; }
};

class DGMHelm3d_26: public HexDGMElement3d {
public:
    DGMHelm3d_26(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 52; }
};

class DGMHelm3d_26t: public TetraDGMElement3d {
public:
    DGMHelm3d_26t(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 52; }
};

class DGMHelm3d_26p: public PrismDGMElement3d {
public:
    DGMHelm3d_26p(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 52; }
};

class DGMHelm3d_26pd: public PyramidDGMElement3d {
public:
    DGMHelm3d_26pd(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 52; }
};

class DGMHelm3d_56: public HexDGMElement3d {
public:
    DGMHelm3d_56(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 53; }
};

class DGMHelm3d_56t: public TetraDGMElement3d {
public:
    DGMHelm3d_56t(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 53; }
};

class DGMHelm3d_56p: public PrismDGMElement3d {
public:
    DGMHelm3d_56p(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 53; }
};

class DGMHelm3d_56pd: public PyramidDGMElement3d {
public:
    DGMHelm3d_56pd(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 53; }
};

class DGMHelm3d_98: public HexDGMElement3d {
public:
    DGMHelm3d_98(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 54; }
};


class DEMHelm3d: public HexDGMElement3d {
public:
    DEMHelm3d(int _o, int* nodenums);
    int defaultLMType() const override { return 51; }
    bool dgmFlag() const override { return false; }
    int nPolynomialDofs() const  { return o*o*o; }

    virtual void createM(complex<double>*);
    virtual void interfMatrix(int fi, DEMElement*,complex<double>*);
    virtual void createRHS(complex<double>*);
};


class DEMHelm3d_6: public DEMHelm3d {
public:
    DEMHelm3d_6(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 51; }
};

class DEMHelm3d_26: public DEMHelm3d {
public:
    DEMHelm3d_26(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 52; }
};

class DEMHelm3d_56: public DEMHelm3d {
public:
    DEMHelm3d_56(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 53; }
};

class DEMHelm3d_98: public DEMHelm3d {
public:
    DEMHelm3d_98(int _o, int* nodenums);
    virtual void dir(int,complex<double>*);
    int defaultLMType() const override { return 54; }
};

#endif
