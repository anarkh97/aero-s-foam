#ifndef _DEMELEMENT_H_
#define _DEMELEMENT_H_

#include <math.h>

#ifndef SANDIA

#include <Element.d/Element.h>
#include <Math.d/ComplexD.h>

#include <Driver.d/GeoSource.h>
extern GeoSource *geoSource;
#include <Driver.d/Domain.h>
extern Domain *domain;

#else
#include <complex>
using std::complex;
#endif

class DEMElement;

class DEMLM {
public:
 DEMElement *e1,*e2;
 int nOffset;
 virtual int type()=0;
 virtual int nDofs()=0;
 virtual void init()=0;
 int setNodeOffset(int j) { nOffset = j; return(0); }
 int getNodeOffset() { return nOffset; }
};

// 
//  | A  B^T | | a | = | d |
//  | B   C  | | c |   | f |
//  C is the condensed part
// default setup condenses out all enrichment dofs

class DEMMatrices {
public:
 DEMMatrices() { B=C=f=0; }
 complex<double>* B;
 complex<double>* C;
 complex<double>* f;
};


class DEMCoreElement {
public:
 int *nn; // geometrical nodes
 bool storeMatricesF;
 bool condensedF;
 DEMLM **lm;
 int *bc;
 virtual ~DEMCoreElement() {
   if (bc) delete[] bc;
//   if (lm) delete[] lm;
   if (nn) delete[] nn;
   if (demm.B) delete[] demm.B;
   if (demm.C) delete[] demm.C;
   if (demm.f) delete[] demm.f;
 }
 virtual int defaultLMType()=0;
 virtual bool condensedFlag()=0;
 virtual bool dgmFlag()=0;
 virtual bool storeMatricesFlag()=0;
 virtual int nPolynomialDofs()=0;
 virtual int nEnrichmentDofs()=0;
 virtual int nGeomNodes()=0;
 virtual int nFaces()=0;
 virtual int nFaceCorners(int fi)=0;
 virtual int * faceCorners(int fi)=0;
 virtual int nLagrangeDofs();
 virtual void createM(complex<double>*)=0;
 virtual void createRHS(complex<double>*)=0;
 virtual void createSol(double *xyz, complex<double>*, complex<double>*)=0;
 virtual void systemMatrix(complex<double>*);
 virtual void interfMatrix(int fi, DEMElement*,complex<double>*);
 virtual void systemRHS(complex<double>*);
 virtual void enrichmentF(double *x, complex<double> *f);
 virtual void polynomialF(double *x, double *f);

 DEMMatrices demm;
 virtual void condensedDofs(int &nc, int *&condensed, int *&kept);
 virtual int createCondensedMatrices(complex<double>* kk, complex<double>* A);

 int *ipiv; // used with one version of static condensation to store pivots
 void staticCondensationLHS(int ne, int nl, complex<double> *kee,
                         complex<double> *kel, complex<double> *kll);
 void staticCondensationRHS(int ne, int nl, int nr,
                         complex<double> *kee, complex<double> *kel,
                         complex<double> *re, complex<double>* rl);
 void staticCondensationSolution(int ne, int nl,
                          complex<double> *kee, complex<double> *kel,
                          complex<double> *pl,
                          complex<double> *re, complex<double> *pe);
};

#ifdef SANDIA

class hex8dem {
public:
    double omega;
    double* nodal_coords;
    double rho;
    double c0;
    complex<double> *forceVector;
    DEMElement *de;

    hex8dem();

    int Initialize(bool condense_enrichment_flag, // indicates whether condensation takes place
                   bool dgm_flag, // indicates whether the element has polynomial field
                   bool store_matrices_flag, // indicates whether element matrices should be stored
                   int element_enrichment_type, // type of enrichment
                   int* lagrange_multiplier_type, // type of lagrange multiplier
                   const double* nodal_coords, // coordinates of element nodes
                   int* connectivity, // global node numbers of element nodes
                   double omega, // omega = 2*pi*f, where f is the frequency
                   double c0, // speed of sound
                   double rho // density of fluid
                   );
    ~hex8dem();
    int NumEhancementDofs() const;
    int NumConstraintDofs() const;
    int NumPolynomialDofs() const;

    void ElemDynamicMtx(complex<double> *matrix) const;
    void ScalarFaceIntegral(int faceid, const double *scalar_field,
                            complex<double> **vector) const;
};


class DEMElement: public DEMCoreElement {
 hex8dem *sinterf;
public:
 virtual int numDofs() {
   return ( (dgmFlag())?0:nPolynomialDofs() ) +
          ( (condensedFlag())?0:nEnrichmentDofs() ) +
           nLagrangeDofs() ; 
 }
 virtual double getOmega() { return sinterf->omega ; }
 virtual double *getWaveDirection() { return 0; }
 virtual complex<double> *getField() { return 0; }
 virtual void getNodalCoord(int n, int *nn, double *xyz) {
   xyz = sinterf->nodal_coords;
 }
 virtual double getRho() { return sinterf->rho; }
 virtual double getSpeedOfSound() {
   return sinterf->c0; }
 virtual double getE() { return 1.0; }
 virtual double getNu() { return 0.3; }
};

#else

class DEMElement: public Element, public DEMCoreElement {
public:
 int ne; // enrichment node offset
 virtual int nodesE(int no);
 virtual void nodalSolution(CoordSet &cs, complex<double>* sol,
                            complex<double>(*nodalSol)[8]);
 complex<double> *forceVector;
// Element functions 
 virtual void renum(int *);
 virtual int numDofs() {
   return ( (dgmFlag())?0:nPolynomialDofs() ) +
          ( (condensedFlag())?0:nEnrichmentDofs() ) +
           nLagrangeDofs() ; 
 }
 virtual int polyDofType()=0;
 virtual int polyDofsPerNode()=0;
 virtual void markDofs(DofSetArray &);
 virtual int* dofs(DofSetArray &, int *p=0);
 virtual int numNodes() {
   return ( (dgmFlag())?0:nGeomNodes() ) +
          ( (condensedFlag())?0:nEnrichmentDofs() )+
          ( nLagrangeDofs() ); 
 }
 virtual int* nodes(int * = 0);
 virtual FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
 virtual FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);

// Interface
 virtual double getOmega() { return geoSource->omega(); }
 virtual double *getWaveDirection() { return domain->getWaveDirection(); }
 virtual complex<double> *getField() { return forceVector; }
 virtual void getNodalCoord(int n, int *nn, double *xyz) {
   (domain->getNodes()).getCoordinates(nn,n,xyz,xyz+n,xyz+2*n);
 }
 virtual double getRho() { return prop->rho; }
 virtual double getSpeedOfSound() {
   return geoSource->omega()/prop->kappaHelm; }
 virtual double getE() { return prop->E; }
 virtual double getNu() { return prop->nu; }
};
#endif

class DEMCoreInterfaceElement {
public:
 int fi;
 DEMElement *deme, *deme2;
 virtual void systemMatrix(complex<double>*);
};

#ifdef SANDIA
class DEMInterfaceElement: public DEMCoreInterfaceElement {
 DEMInterfaceElement(DEMElement *_deme, DEMElement *_deme2, int _fi);
};

#else

class DEMInterfaceElement: public Element, public DEMCoreInterfaceElement {
public:
 DEMInterfaceElement(DEMElement *_deme, DEMElement *_deme2, int _fi);
 virtual void renum(int *);
 virtual int numDofs() { return deme->numDofs()-deme->nLagrangeDofs()+
                                deme2->numDofs()-deme2->nLagrangeDofs(); }
 virtual int* dofs(DofSetArray &, int *p=0);
 virtual void markDofs(DofSetArray &);
 virtual int numNodes() { return deme->numNodes()-deme->nLagrangeDofs()+
                                deme2->numNodes()-deme2->nLagrangeDofs(); }
 virtual int* nodes(int * = 0);
};
#endif

#ifdef _TEMPLATE_FIX_
//#include <Element.d/Helm.d/DEMElementT.C>
#endif

#endif

