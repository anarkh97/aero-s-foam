#ifndef _DEMHELM3D_H_
#define _DEMHELM3D_H_

#include <math.h>
#include <Element.d/DEM.d/DEMElement.h>
#include <Math.d/ComplexD.h>
#include <Utils.d/dofset.h>


class DGMHelm3d_LM: public DEMLM {
public:
 virtual int nDofs()=0;
 virtual int type()=0;
 virtual void init(CoordSet &cs) {};
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
 double * getCubeDir(int n);
public:
 int ndir;
 virtual void dir(int,complex<double>*) {};
 int o;
 DGMHelm3d(int _o, int* nodenums);
 virtual int defaultLMType() { return 51; }
 virtual bool dgmFlag() { return true; }
 virtual bool condensedFlag() {
   for(int i=0;i<nFaces();i++) if (bc[i]==3) return false;
   return true;
 }
 virtual bool storeMatricesFlag() { return true; }
// virtual bool storeMatricesFlag() { return false; }
 virtual int nPolynomialDofs() { return 0; }
 virtual int nEnrichmentDofs() { return ndir; }
 virtual int nGeomNodes() { return (o>0)?o*o*o:((-o)*(-o+1)*(-o+2))/6; }
 virtual int nFaces() { return (o>0)?6:4; }
 virtual int nFaceCorners(int fi) { return (o>0)?4:3; }
 virtual int *faceCorners(int fi);
 virtual int polyDofType() { return DofSet::Helm; }
 virtual int polyDofsPerNode() { return 1; }

 virtual void getRef(double *xyz,double *xy);
 virtual void createM(CoordSet &cs, complex<double>*);
 virtual void createRHS(CoordSet &cs, complex<double>*);
 virtual void createSol(double *xyz, complex<double>*,
                            complex<double>*);
};


class DGMHelm3d_6: public DGMHelm3d {
public:
 DGMHelm3d_6(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 51; }
};


class DGMHelm3d_6t: public DGMHelm3d {
public:
 DGMHelm3d_6t(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 51; }
};

class DGMHelm3d_26: public DGMHelm3d {
public:
 DGMHelm3d_26(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 52; }
};

class DGMHelm3d_26t: public DGMHelm3d {
public:
 DGMHelm3d_26t(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 52; }
};

class DGMHelm3d_56: public DGMHelm3d {
public:
 DGMHelm3d_56(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 53; }
};

class DGMHelm3d_56t: public DGMHelm3d {
public:
 DGMHelm3d_56t(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 53; }
};

class DGMHelm3d_98: public DGMHelm3d {
public:
 DGMHelm3d_98(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 54; }
};


class DEMHelm3d: public DGMHelm3d {
public:
 DEMHelm3d(int _o, int* nodenums);
 virtual int defaultLMType() { return 51; }
 virtual bool dgmFlag() { return false; }
 virtual int nPolynomialDofs() { return o*o*o; }

 virtual void createM(CoordSet &cs, complex<double>*);
 virtual void createRHS(CoordSet &cs, complex<double>*);
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

