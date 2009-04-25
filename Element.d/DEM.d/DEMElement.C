#include <stdio.h>
#include <stdlib.h>
#include <alloca.h>
#include <math.h>
#include <time.h>
#include <complex>

#include <Element.d/Helm.d/IsoParamUtils.h>
#include <Element.d/Helm.d/GaussRules.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

#include <Element.d/DEM.d/DEMElement.h>

using namespace std;

DEMInterfaceElement::DEMInterfaceElement(DEMElement *_deme, DEMElement *_deme2, int _fi) {
 deme = _deme;
 deme2 = _deme2;
 fi = _fi;
}


// Element functions
void DEMElement::renum(int *table) {
 fprintf(stderr,"DEMElement::renum is not implemented.\n");
}


void DEMInterfaceElement::renum(int *table) {
 fprintf(stderr,"DEMInterfaceElement::renum is not implemented.\n");
}


int* DEMElement::nodes(int *p) {

 if (p == 0) p = new int[numNodes()];
 int i,ii=0;
 if (!dgmFlag())
   for(i=0;i<nGeomNodes();i++) p[ii++] = nn[i];
 for(i=0;i<nFaces();i++) {
   if (lm[i] != 0) {
     int nlm = lm[i]->nDofs();
     for(int j=0;j<nlm;j++) p[ii++] = lm[i]->getNodeOffset()+j;
   }
 }
 if (!condensedFlag()) for(i=0;i<nEnrichmentDofs();i++) p[ii++] = ne+i;

 return p;
}


int* DEMInterfaceElement::nodes(int *p) {

 if (p == 0) p = new int[numNodes()];

 int ii=0;
 if (!deme->dgmFlag()) 
   for(int i=0;i<deme->nGeomNodes();i++) p[ii++] = deme->nn[i];
 for(int i=0;i<deme->nEnrichmentDofs();i++) p[ii++] = deme->ne+i;
 if (!deme2->dgmFlag()) 
   for(int i=0;i<deme2->nGeomNodes();i++) p[ii++] = deme2->nn[i];
 for(int i=0;i<deme2->nEnrichmentDofs();i++) p[ii++] = deme2->ne+i;

 return p;
}


int* DEMElement::dofs(DofSetArray &dsa, int *p) {

 if (p == 0) p = new int[numDofs()];
 int ii=0;
 if (!dgmFlag()) {
   for(int i=0;i<nGeomNodes();i++) { 
     dsa.number(nn[i],polyDofType(),p+ii);
     ii += polyDofsPerNode();
   }
 }

 for(int i=0;i<nFaces();i++) {
   if (lm[i] != 0) {
     int j;
     for(j=0;j<lm[i]->nDofs();j++) {
        dsa.number(lm[i]->getNodeOffset()+j, DofSet::Lagrange1,p+ii++);
       
     }
   }
 }
 
 if (!condensedFlag()) 
// RT: Temporarily put DofSet::Helm, but should create some other new name
   for(int i=0;i<nEnrichmentDofs();i++) dsa.number(ne+i,DofSet::Helm,p+ii++);
 return p;
}


int* DEMInterfaceElement::dofs(DofSetArray &dsa, int *p) {

 if (p == 0) p = new int[numDofs()];
 int ii=0;

 if (!deme->dgmFlag()) {
   for(int i=0;i<deme->nGeomNodes();i++) {
     dsa.number(deme->nn[i], deme->polyDofType(),p+ii);
     ii += deme->polyDofsPerNode();
   }
 }

// RT: Temporarily put DofSet::Helm, but should create some other new name
 for(int i=0;i<deme->nEnrichmentDofs();i++)
   dsa.number(deme->ne+i,DofSet::Helm,p+ii++);

 if (!deme2->dgmFlag()) {
   for(int i=0;i<deme2->nGeomNodes();i++) {
     dsa.number(deme2->nn[i], deme2->polyDofType(),p+ii);
     ii += deme2->polyDofsPerNode();
   }
 }

// RT: Temporarily put DofSet::Helm, but should create some other new name
 for(int i=0;i<deme2->nEnrichmentDofs();i++)
   dsa.number(deme2->ne+i,DofSet::Helm,p+ii++);

 return p;
}


void DEMElement::markDofs(DofSetArray &dsa) {

 if (!dgmFlag())
   for(int i=0;i<nGeomNodes();i++) { dsa.mark(nn[i],polyDofType());
 }

 for(int i=0;i<nFaces();i++) {
   if (lm[i] != 0) {
     int j;
     for(j=0;j<lm[i]->nDofs();j++) {
        dsa.mark(lm[i]->getNodeOffset()+j, DofSet::Lagrange1);
     }
   }
 }

 if (!condensedFlag()) 
// RT: Temporarily put DofSet::Helm, but should create some other new name
   for(int i=0;i<nEnrichmentDofs();i++) dsa.mark(ne+i,DofSet::Helm);
}


void DEMInterfaceElement::markDofs(DofSetArray &dsa) {

 if (!deme->dgmFlag())
   for(int i=0;i<deme->nGeomNodes();i++) dsa.mark(deme->nn[i],deme->polyDofType());

// RT: Temporarily put DofSet::Helm, but should create some other new name
 for(int i=0;i<deme->nEnrichmentDofs();i++) dsa.mark(deme->ne+i,DofSet::Helm);

 if (!deme2->dgmFlag())
   for(int i=0;i<deme2->nGeomNodes();i++) dsa.mark(deme2->nn[i],deme2->polyDofType());

// RT: Temporarily put DofSet::Helm, but should create some other new name
 for(int i=0;i<deme2->nEnrichmentDofs();i++) dsa.mark(deme2->ne+i,DofSet::Helm);
}


FullSquareMatrix DEMElement::massMatrix(CoordSet &cs, double *K, int fl) {
 fprintf(stderr,"DEMElement::massMatrix not implemented.\n");
 FullSquareMatrix ret(0,K);
 return ret;
}


FullSquareMatrix DEMElement::stiffness(CoordSet &cs, double *K, int flg ) {
 fprintf(stderr,"DEMElement::massMatrix not implemented.\n");
 FullSquareMatrix ret(0,K);
 return ret;
}


int DEMElement::nodesE(int no) {
 ne = no;
 return condensedFlag()?0:nEnrichmentDofs();
}

int DEMElement::nLagrangeDofs() {

 int nl = 0;
 for(int i=0;i<nFaces();i++)
   if (lm[i] != 0) nl += lm[i]->nDofs();
 return nl;
}


void DEMElement::condensedDofs(int &nc, int *&condensed, int *&kept) {
 int n = nPolynomialDofs()+nLagrangeDofs()+nEnrichmentDofs();
 nc = nEnrichmentDofs();
 condensed = new int[nc];
 kept = new int[n-nc];
 int i;
 for(i=0;i<nc;i++)
   condensed[i] = nPolynomialDofs()+nLagrangeDofs()+i;
 for(i=0;i<n-nc;i++)
   kept[i] = i;
}


int DEMElement::createCondensedMatrices(complex<double>* kk, 
                                         complex<double>* A) {

 int nc, *condensed, *kept;
 int n = nPolynomialDofs()+nLagrangeDofs()+nEnrichmentDofs();
 condensedDofs(nc,condensed,kept);
 demm.B = new complex<double>[nc*(n-nc)];
 demm.C = new complex<double>[nc*nc];
 int i,j;
 for(i=0;i<nc;i++) for(j=0;j<nc;j++)
   demm.C[i+j*nc] = kk[condensed[i]+n*condensed[j]];
 for(i=0;i<n-nc;i++) for(j=0;j<n-nc;j++)
   A[i+j*(n-nc)] = kk[kept[i]+n*kept[j]];
 for(i=0;i<n-nc;i++) for(j=0;j<nc;j++)
   demm.B[j+i*nc] = kk[kept[i]+n*condensed[j]];
 delete[] condensed;
 delete[] kept;

 staticCondensationLHS(nc,n-nc,demm.C,demm.B,A);
 return nc;
}
  


void DEMElement::systemMatrix(CoordSet &cs,
                              complex<double> *K) {

 int n = nPolynomialDofs()+nLagrangeDofs()+nEnrichmentDofs();
 if (!condensedFlag()) {
   createM(cs, K);
 } else {
   complex<double> *kk;
   kk = new complex<double>[n*n];

   createM(cs, kk);
   int nc = createCondensedMatrices(kk,K);
   delete[] kk;
   
   if (!storeMatricesFlag()) {
     delete[] demm.B;
     delete[] demm.C;
   }
 }
}


void DEMInterfaceElement::systemMatrix(CoordSet &cs,
                              complex<double> *K) {
 deme->interfMatrix(cs,fi,deme2,K);
}


void DEMElement::interfMatrix(CoordSet &cs, int fi, DEMElement* deme2,
                              complex<double> *K) {
 fprintf(stderr,"DEMElement::interfMatrix not implemented.\n");
}



void DEMElement::nodalSolution(CoordSet &cs, complex<double> *sol,
                          complex<double> (*nodalSol)[8]) {


 complex<double> *a;
 if (condensedFlag()) {
   int n = nPolynomialDofs()+nLagrangeDofs()+nEnrichmentDofs();
   int nc, *condensed, *kept;
   condensedDofs(nc,condensed,kept);
   if (!storeMatricesFlag()) {
     complex<double> *kk;
     kk = new complex<double>[n*n];
     createM(cs, kk);
     complex<double> *A = new complex<double>[nc*nc];
     createCondensedMatrices(kk,A);
     delete[] A;
     delete[] kk;
   }
   a = new complex<double>[nc]; 
   complex<double>* c = new complex<double>[n-nc]; 
   for(int i=0;i<n-nc;i++) c[i] = sol[kept[i]];
   delete[] kept;
   delete[] condensed;
   staticCondensationSolution(nc, n-nc, demm.C, demm.B,
                              c, demm.f, a);
   delete[] c;
   if (!storeMatricesFlag()) {
     delete[] demm.B;
     delete[] demm.C;
   }
 } else {
   a = sol + nPolynomialDofs()+nLagrangeDofs();
 }
 int j=0;
 for(int i=0;i<nGeomNodes();i++) {
   double xyz[3];
   cs.getCoordinates(nn+i,1,xyz,xyz+1,xyz+2);
   createSol(xyz, a, nodalSol[i]);
   if (!dgmFlag()) {
      if (polyDofType()& DofSet::Xdisp) { 
        nodalSol[i][1] += sol[j];
        j++;
      }
      if (polyDofType()& DofSet::Ydisp)  {
        nodalSol[i][2] += sol[j];
        j++;
      }
      if (polyDofType()& DofSet::Zdisp) {
        nodalSol[i][3] += sol[j];
        j++;
      }
      if (polyDofType()& DofSet::Helm) {
        nodalSol[i][0] += sol[j];
        j++;
      }
   }
 }
 if (condensedFlag()) delete[] a;
}


void DEMElement::enrichmentF(CoordSet &cs, double *x, complex<double> *f) {
 fprintf(stderr,"DEMElement::enrichmentF is not implemented.\n");
}

void DEMElement::polynomialF(CoordSet &cs, double *x, double *f) {
 fprintf(stderr,"DEMElement::polynomialF is not implemented.\n");
}


void DEMElement::systemRHS(CoordSet &cs, complex<double> *v) {


 if (!condensedFlag()) {
   createRHS(cs,v);
 } else {
   int n = nPolynomialDofs()+nLagrangeDofs()+nEnrichmentDofs();
   complex<double> *vv = new complex<double>[n];
   createRHS(cs,vv);
   int nc, *condensed, *kept;
   condensedDofs(nc,condensed,kept);
   if (!storeMatricesFlag()) {
     complex<double> *kk;
     kk = new complex<double>[n*n];
     createM(cs, kk);
     complex<double> *A = new complex<double>[nc*nc];
     createCondensedMatrices(kk,A);
     delete[] A;
     delete[] kk;
   }

   if (demm.f==0)
     demm.f = new complex<double>[nc];
   int i;
   for(i=0;i<nc;i++) demm.f[i] = vv[condensed[i]];
   for(i=0;i<n-nc;i++) v[i] = vv[kept[i]];
   delete[] condensed;
   delete[] kept;
   delete[] vv;
   staticCondensationRHS(nc, n-nc, 1, demm.C, demm.B, demm.f, v);
   if (!storeMatricesFlag()) {
     delete[] demm.B;
     delete[] demm.C;
   }
 }
}


extern "C" {
void _FORTRAN(zgesvx)( const char &fact, const char & trans, const int &n,
                       const int &nrhs, ComplexD *a, const int &lda,
                       ComplexD *af, const int &ldaf, int *ipiv,
                       char &equed, double *r, double *c,
                       ComplexD *b, const int &ldb,
                       ComplexD *x, const int &ldx,
                       double &rcond,
                       double *ferr, double *berr,
                       ComplexD  *work, double *rwork, int &info);
}

extern "C" {
void _FORTRAN(zgesv)
(const int &n, const int &nrhs,
                            ComplexD *a, const int &lda, int *ipiv,
                            ComplexD *b, const int &ldb, int &info);
}

extern "C" {
void _FORTRAN(zgetrf)
(const int &m, const int &n, ComplexD *a, const int &lda, int *ipiv, int &info);
}

extern "C" {
void _FORTRAN(zgetrs) (const char & trans, const int &n, const int &nrhs,
                            ComplexD *a, const int &lda, int *ipiv,
                            ComplexD *b, const int &ldb, int &info);
}

extern "C" {
void _FORTRAN(zgemv)(const char &trans, const int &m, const int &n,
                          const ComplexD &lapha, ComplexD *a, const int &lda,
                          ComplexD *x, const int &incx,
                          const ComplexD &beta, ComplexD *y, const int &incy);
}

extern "C" {
void _FORTRAN(zgemm)(const char &transa, const char &transb,
                            const int &m,const int &n, const int &k,
                           const ComplexD &alpha, ComplexD *a, const int &lda,
                           ComplexD *b, const int &ldb,
                           const ComplexD &beta, ComplexD *c, const int &ldc);
}

//#define RT_USE_ALLOCA
#define LAFVERSION 1

#if LAFVERSION == 0
void DEMElement::staticCondensationLHS(int ne, int nl, complex<double> *kee,
                         complex<double> *kel, complex<double> *kll) {
#ifdef RT_USE_ALLOCA
 complex<double> *af = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 complex<double> *x =  (complex<double>*)alloca(sizeof(complex<double>)*ne*nl);
#else
 complex<double> *af = new complex<double>[ne*ne];
 complex<double> *x = new complex<double>[ne*nl]; 
#endif

 int info;
 int lda = ne;
 int ldaf = ne;
 int ldb = ne;
 int ldx = ne;
#ifdef RT_USE_ALLOCA
 double *r = (double*)alloca(sizeof(double)*ne);
 double *c = (double*)alloca(sizeof(double)*ne);
 int *ipiv = (int*)alloca(sizeof(int)*ne);
 double *ferr = (double*)alloca(sizeof(double)*nl);
 double *berr = (double*)alloca(sizeof(double)*nl);
 complex<double> *work = (complex<double>*)alloca(sizeof(complex<double>)*2*ne);
 double *rwork = (double*)alloca(sizeof(double)*2*ne);
#else
 double *r = new double[ne];
 double *c = new double[ne];
 int *ipiv = new int[ne];
 double *ferr = new double[nl];
 double *berr = new double[nl];
 complex<double> *work =  new complex<double>[2*ne];
 double *rwork = new double[2*ne];
#endif
// double *r = 0, *c = 0;
 char fact = 'N';
 char trans = 'N';
 char equed = 'N';
 double rcond;
 _FORTRAN(zgesvx)(fact, trans, ne, nl, kee, lda, af, ldaf, ipiv,
                  equed, r, c, kel, ldb, x, ldx,
                  rcond, ferr, berr, work, rwork, info);
//  fprintf(stderr, " rcond = %e \n", rcond);
 if (rcond<1e-15) fprintf(stderr,
   "Nearly singular matrix: rcond = %e in DEMElement::staticCondensationLHS\n",
    rcond);
 if (info!=0) fprintf(stderr,
   "Singular matrix in DEMElement::staticCondensationLHS\n");

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nl, ne, alpha , kel, lda, x, ldx,
                  beta, kll, ldc);
#ifndef RT_USE_ALLOCA
 delete[] af;
 delete[] x;
 delete[] r;
 delete[] c;
 delete[] ipiv;
 delete[] ferr;
 delete[] berr;
 delete[] work;
 delete[] rwork;
#endif
}

#elif LAFVERSION == 1
void DEMElement::staticCondensationLHS(int ne, int nl, complex<double> *kee,
                         complex<double> *kel, complex<double> *kll) {

 ipiv = new int[ne];

 int info;
 int lda = ne;
 _FORTRAN(zgetrf)(ne, ne, kee, lda, ipiv, info);
 if (info>0) fprintf(stderr,
   "Singular matrix in DEMElement::staticCondensationLHS\n");

 char trans = 'N';
 int ldb = ne;
 complex<double> *b = (complex<double>*)alloca(sizeof(complex<double>)*ne*nl);
 for(int i=0;i<ne*nl;i++) b[i] = kel[i];
 _FORTRAN(zgetrs)(trans, ne, nl, kee, lda, ipiv, b, ldb, info);
 if (info!=0) fprintf(stderr,
   "Error by zgetrs in DEMElement::staticCondensationLHS\n");

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nl, ne, alpha , kel, lda, b, ldb,
                  beta, kll, ldc);
}
#else
void DEMElement::staticCondensationLHS(int ne, int nl, complex<double> *kee,
                         complex<double> *kel, complex<double> *kll) {

 int *ipiv = (int*)alloca(sizeof(int)*ne);
 complex<double> *a = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 int i;
 for(i=0;i<ne*ne;i++) a[i] = kee[i];
 complex<double> *b = (complex<double>*)alloca(sizeof(complex<double>)*ne*nl);
 for(i=0;i<ne*nl;i++) b[i] = kel[i];

 int info;
 int lda = ne;
 int ldb = ne;
 _FORTRAN(zgesv)(ne, nl, a, lda, ipiv, b, ldb, info);

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nl, ne, alpha , kel, lda, b, ldb,
                  beta, kll, ldc);
}
#endif


#if LAFVERSION == 0
void DEMElement::staticCondensationRHS(int ne, int nl, int nr,
                         complex<double> *kee, complex<double> *kel,
                         complex<double> *re, complex<double>* rl) {

 complex<double> *af = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 complex<double> *x = (complex<double>*)alloca(sizeof(complex<double>)*ne*nr);

 int info;
 int lda = ne;
 int ldaf = ne;
 int ldb = ne;
 int ldx = ne;
 // double *r = (double*)alloca(sizeof(double)*ne);
 // double *c = (double*)alloca(sizeof(double)*ne);
 double *r = 0, *c = 0;
 char fact = 'N';
 char trans = 'N';
 char equed = 'N';
 double rcond;
 int *ipiv = (int*)alloca(sizeof(int)*ne);
 double *ferr = (double*)alloca(sizeof(double)*nr);
 double *berr = (double*)alloca(sizeof(double)*nr);
 complex<double> *work = (complex<double>*)alloca(sizeof(complex<double>)*2*ne);
 double *rwork = (double*)alloca(sizeof(double)*2*ne);
 _FORTRAN(zgesvx)(fact, trans, ne, nr, kee, lda, af, ldaf, ipiv,
                  equed, r, c, re, ldb, x, ldx,
                  rcond, ferr, berr, work, rwork, info);

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nr, ne, alpha , kel, lda, x, ldx,
                  beta, rl, ldc);
}

#elif LAFVERSION == 1
void DEMElement::staticCondensationRHS(int ne, int nl, int nr,
                         complex<double> *kee, complex<double> *kel,
                         complex<double> *re, complex<double>* rl) {

 int info;
 char trans = 'N';
 int lda = ne;
 int ldb = ne;
 complex<double> *b = (complex<double>*)alloca(sizeof(complex<double>)*ne*nr);
 for(int i=0;i<ne*nr;i++) b[i] = re[i];
 _FORTRAN(zgetrs)(trans, ne, nr, kee, lda, ipiv, b, ldb, info);
 if (info!=0) fprintf(stderr,
   "Error by zgetrs in DEMElement::staticCondensationRHS\n");

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nr, ne, alpha , kel, lda, b, ldb,
                  beta, rl, ldc);
}

#else
void DEMElement::staticCondensationRHS(int ne, int nl, int nr,
                         complex<double> *kee, complex<double> *kel,
                         complex<double> *re, complex<double>* rl) {

 int *ipiv = (int*)alloca(sizeof(int)*ne);
 complex<double> *a = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 int i;
 for(i=0;i<ne*ne;i++) a[i] = kee[i];
 complex<double> *b = (complex<double>*)alloca(sizeof(complex<double>)*ne*nr);
 for(i=0;i<ne*nr;i++) b[i] = re[i];

 int info;
 int lda = ne;
 int ldb = ne;
 _FORTRAN(zgesv)(ne, nr, a, lda, ipiv, b, ldb, info);

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nr, ne, alpha , kel, lda, b, ldb,
                  beta, rl, ldc);
}
#endif


#if LAFVERSION == 0
void DEMElement::staticCondensationSolution(int ne, int nl,
                          complex<double> *kee, complex<double> *kel,
                          complex<double> *pl,
                          complex<double> *re, complex<double> *pe) {
 int i;
 for(i=0;i<ne;i++) pe[i] = re[i];

 char trans = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int incx = 1, incy = 1;
 _FORTRAN(zgemv) ( trans, ne, nl, alpha, kel, ne, pl, incx, beta, pe, incy);

 int nr = 1;

 complex<double> *af = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 complex<double> *x = (complex<double>*)alloca(sizeof(complex<double>)*ne*nr);
 for(int i=0;i<ne;i++) x[i] = pe[i];

 int info;
 int lda = ne;
 int ldaf = ne;
 int ldb = ne;
 int ldx = ne;
 // double *r = (double*)alloca(sizeof(double)*ne);
 // double *c = (double*)alloca(sizeof(double)*ne);
 double *r = 0, *c = 0;
 char fact = 'N';
 trans = 'N';
 char equed = 'N';
 double rcond;
 int *ipiv = (int*)alloca(sizeof(int)*ne);
 double *ferr = (double*)alloca(sizeof(double)*nr);
 double *berr = (double*)alloca(sizeof(double)*nr);
 complex<double> *work = (complex<double>*)alloca(sizeof(complex<double>)*2*ne);
 double *rwork = (double*)alloca(sizeof(double)*2*ne);
 _FORTRAN(zgesvx)(fact, trans, ne, nr, kee, lda, af, ldaf, ipiv,
                  equed, r, c, x, ldb, pe, ldx,
                  rcond, ferr, berr, work, rwork, info);
 if (rcond<1e-15) fprintf(stderr,
   "Nearly singular matrix: rcond = %e in DEMElement::staticCondensationSolution\n",
    rcond);
 if (info!=0) fprintf(stderr,
    "Singular matrix in DEMElement::staticCondensationSolution\n");
}

#elif LAFVERSION == 1
void DEMElement::staticCondensationSolution(int ne, int nl,
                          complex<double> *kee, complex<double> *kel,
                          complex<double> *pl,
                          complex<double> *re, complex<double> *pe) {

 int i;
 for(i=0;i<ne;i++) pe[i] = re[i];

 char trans = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int incx = 1, incy = 1;
 _FORTRAN(zgemv) ( trans, ne, nl, alpha, kel, ne, pl, incx, beta, pe, incy);

 int info;
 int lda = ne;
 int ldb = ne;
 int nr = 1;

 _FORTRAN(zgetrs)(trans, ne, nr, kee, lda, ipiv, pe, ldb, info);
}

#else
void DEMElement::staticCondensationSolution(int ne, int nl,
                          complex<double> *kee, complex<double> *kel,
                          complex<double> *pl,
                          complex<double> *re, complex<double> *pe) {

 int i;
 for(i=0;i<ne;i++) pe[i] = re[i];

 char trans = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int incx = 1, incy = 1;
 _FORTRAN(zgemv) ( trans, ne, nl, alpha, kel, ne, pl, incx, beta, pe, incy);

 int *ipiv = (int*)alloca(sizeof(int)*ne);
 complex<double> *a = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 for(i=0;i<ne*ne;i++) a[i] = kee[i];

 int info;
 int lda = ne;
 int ldb = ne;
 int nr = 1;
 _FORTRAN(zgesv)(ne, nr, a, lda, ipiv, pe, ldb, info);
}
#endif

