#ifndef _DEMELEMENT_H_
#define _DEMELEMENT_H_

#include <math.h>
#include <Element.d/Element.h>
#include <Math.d/ComplexD.h>

class DEMElement;

class DEMLM {
public:
 DEMElement *e1,*e2;
 int nOffset;
 virtual int type()=0;
 virtual int nDofs()=0;
 virtual void init(CoordSet &cs)=0;
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



class DEMElement: public Element {
public:
 int *nn; // geometrical nodes
 int ne; // enrichment node offset
 DEMLM **lm;
 int *bc;
 virtual ~DEMElement() {
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
 virtual int polyDofType()=0;
 virtual int polyDofsPerNode()=0;
 virtual void createM(CoordSet &cs, complex<double>*)=0;
 virtual void createRHS(CoordSet &cs, complex<double>*)=0;
 virtual void createSol(double *xyz, complex<double>*, complex<double>*)=0;
 virtual void systemMatrix(CoordSet &cs, complex<double>*);
 virtual void interfMatrix(CoordSet &cs, int fi, DEMElement*,complex<double>*);
 virtual void systemRHS(CoordSet &cs, complex<double>*);
 virtual void nodalSolution(CoordSet &cs, complex<double>* sol,
                            complex<double>(*nodalSol)[8]);
 virtual void enrichmentF(CoordSet &cs, double *x, complex<double> *f);
 virtual void polynomialF(CoordSet &cs, double *x, double *f);

 virtual int nodesE(int no);

 virtual void condensedDofs(int &nc, int *&condensed, int *&kept);
 virtual int createCondensedMatrices(complex<double>* kk, complex<double>* A);

 DEMMatrices demm;

// Element functions 
 virtual void renum(int *);
 virtual int numDofs() {
   return ( (dgmFlag())?0:nPolynomialDofs() ) +
          ( (condensedFlag())?0:nEnrichmentDofs() ) +
           nLagrangeDofs() ; 
 }
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

 complex<double> DEMSolution(double *xyz, complex<double>* sol,
                                         complex<double> *pe);
};


class DEMInterfaceElement: public Element {
public:
 int fi;
 DEMElement *deme, *deme2;
 DEMInterfaceElement(DEMElement *_deme, DEMElement *_deme2, int _fi);

 virtual void renum(int *);
 virtual int numDofs() { return deme->numDofs()-deme->nLagrangeDofs()+
                                deme2->numDofs()-deme2->nLagrangeDofs(); }
 virtual int* dofs(DofSetArray &, int *p=0);
 virtual void markDofs(DofSetArray &);
 virtual int numNodes() { return deme->numNodes()-deme->nLagrangeDofs()+
                                deme2->numNodes()-deme2->nLagrangeDofs(); }
 virtual int* nodes(int * = 0);

 virtual void systemMatrix(CoordSet &cs, complex<double>*);
};

#ifdef _TEMPLATE_FIX_
//#include <Element.d/Helm.d/DEMElementT.C>
#endif

#endif

