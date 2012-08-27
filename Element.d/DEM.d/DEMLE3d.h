#ifndef _DEMLE3D_H_
#define _DEMLE3D_H_

#include <cmath>
#include <Element.d/DEM.d/DEMElement.h>

class DGMLE3d_LM: public DEMLM {
public:
 virtual int nDofs()=0;
 virtual int type()=0;
 virtual void init() {};
 virtual void ldir(int,double*,double*,complex<double>*)=0;
};


class DGMLE3d_3_LM: public DGMLE3d_LM {
public:
 virtual int nDofs() { return 3; }
 virtual int type() { return 251; }
 virtual void ldir(int,double*,double*,complex<double>*);
};

class DGMLE3d_15_LM: public DGMLE3d_LM {
public:
 virtual int nDofs() { return 15; }
 virtual int type() { return 252; }
 virtual void ldir(int,double*,double*,complex<double>*);
};

class DGMLE3d_28_LM: public DGMLE3d_LM {
public:
 virtual int nDofs() { return 28; }
 virtual int type() { return 253; }
 virtual void ldir(int,double*,double*,complex<double>*);
};

class DGMLE3d: public DEMElement {
public:
 int ndir;
 virtual void dir(int,complex<double>*) {};
 int o;
 DGMLE3d(int _o, int* nodenums);
 virtual int defaultLMType() { return 251; }
 virtual bool dgmFlag() { return true; }
 virtual bool condensedFlag() {
   for(int i=0;i<nFaces();i++) if (bc[i]==3) return false;
   return true;
 }
 virtual bool storeMatricesFlag() { return true; }
 virtual int nPolynomialDofs() { return 0; }
 virtual int nEnrichmentDofs() { return ndir; }
 virtual int nGeomNodes() { return o*o*o; }
 virtual int nFaces() { return 6; }
 virtual int nFaceCorners(int fi) { return 4; }
 virtual int *faceCorners(int fi) { 
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
 }
 virtual int polyDofType() { return DofSet::Xdisp|DofSet::Ydisp|DofSet::Zdisp; }
 virtual int polyDofsPerNode() { return 3; }

 virtual void getRef(double *xyz,double *xy);
 virtual void createM(complex<double>*);
 virtual void createRHS(complex<double>*);
 virtual void createSol(double *xyz, complex<double>*,
                            complex<double>*);

 virtual void enrichmentF(double *x, complex<double> *f);
 virtual void polynomialF(double *x, double *f);
};


class DGMLE3d_6: public DGMLE3d {
public:
 DGMLE3d_6(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 251; }
};

class DGMLE3d_26: public DGMLE3d {
public:
 DGMLE3d_26(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 252; }
};

class DGMLE3d_50: public DGMLE3d {
public:
 DGMLE3d_50(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 253; }
};


class DEMLE3d: public DGMLE3d {
public:
 DEMLE3d(int _o, int* nodenums);
 virtual bool dgmFlag() { return false; }
 virtual int nPolynomialDofs() { return 3*o*o*o; }

 virtual void createM(complex<double>*);
 virtual void createRHS(complex<double>*);
};

class DEMLE3d_6: public DEMLE3d {
public:
 DEMLE3d_6(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 251; }
};

class DEMLE3d_26: public DEMLE3d {
public:
 DEMLE3d_26(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 252; }
};

class DEMLE3d_50: public DEMLE3d {
public:
 DEMLE3d_50(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 253; }
};

#endif
