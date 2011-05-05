#ifndef _DEMLE2D_H_
#define _DEMLE2D_H_

#include <cmath>
#include <Element.d/DEM.d/DEMElement.h>

class DGMLE2d_LM: public DEMLM {
public:
 virtual int nDofs()=0;
 virtual int type()=0;
 virtual void init() {};
 virtual void ldir(int,double[2],complex<double>*)=0;
};


class DGMLE2d_1_LM: public DGMLE2d_LM {
public:
 virtual int nDofs() { return 2; }
 virtual int type() { return 201; }
 virtual void ldir(int,double[2],complex<double>*);
};

class DGMLE2d_4_LM: public DGMLE2d_LM {
public:
 virtual int nDofs() { return 4; }
 virtual int type() { return 202; }
 virtual void ldir(int,double[2],complex<double>*);
};

class DGMLE2d: public DEMElement {
public:
 int ndir;
 virtual void dir(int,complex<double>*) {};
 int o;
 DGMLE2d(int _o, int* nodenums);
 virtual int defaultLMType() { return 201; }
 virtual bool dgmFlag() { return true; }
 virtual bool condensedFlag() {
   for(int i=0;i<nFaces();i++) if (bc[i]==3) return false;
   return true;
//   return false;
 }
 virtual bool storeMatricesFlag() { return true; }
 virtual int nPolynomialDofs() { return 0; }
 virtual int nEnrichmentDofs() { return ndir; }
 virtual int nGeomNodes() { return o*o; }
 virtual int nFaces() { return 4; }
 virtual int nFaceCorners(int fi) { return 2; }
 virtual int *faceCorners(int fi) { 
   int *fc = new int[2];
   if (fi==1) { fc[0] = nn[0]; fc[1] = nn[o-1]; }
   else if (fi==2) { fc[0] = nn[o-1]; fc[1] = nn[o*o-1]; }
   else if (fi==3) { fc[0] = nn[o*o-1]; fc[1] = nn[o*(o-1)]; }
   else if (fi==4) { fc[0] = nn[o*(o-1)]; fc[1] = nn[0]; }
   return fc;
 }
 virtual int polyDofType() { return DofSet::Xdisp|DofSet::Ydisp; }
 virtual int polyDofsPerNode() { return 2; }

 virtual void getRef(double *xyz,double *xy);
 virtual void createM(complex<double>*);
 virtual void createRHS(complex<double>*);
 virtual void createSol(double *xyz, complex<double>*,
                            complex<double>*);

 virtual void enrichmentF(double *x, complex<double> *f);
 virtual void polynomialF(double *x, double *f);
};


class DGMLE2d_4: public DGMLE2d {
public:
 DGMLE2d_4(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 201; }
};

class DGMLE2d_16: public DGMLE2d {
public:
 DGMLE2d_16(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 202; }
};


class DEMLE2d: public DGMLE2d {
public:
 DEMLE2d(int _o, int* nodenums);
 virtual bool dgmFlag() { return false; }
 virtual int nPolynomialDofs() { return 2*o*o; }

 virtual void createM(complex<double>*);
 virtual void createRHS(complex<double>*);
};

class DEMLE2d_4: public DEMLE2d {
public:
 DEMLE2d_4(int _o, int* nodenums);
 virtual void dir(int,complex<double>*);
 virtual int defaultLMType() { return 201; }
};

#endif

