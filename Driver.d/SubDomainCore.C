#include <stdio.h>
#include <Driver.d/SubDomain.h>
#include <Driver.d/CornerMaker.h>

extern int salinasFlag;

extern "C"      {
void _FORTRAN(cfjacobi)(double *,double *,double *, double *,int&,double&,int &);
}

void
getJacobi(double *kappa, double * mu, FullSquareMatrix &xx,
          double *eigVal,int nsmax, int subSpaceSize, double tolJac)
{
  int i,j;

  _FORTRAN(cfjacobi)(kappa,mu,xx[0],eigVal,nsmax,tolJac,subSpaceSize);

  // sort eigenvalues.

  int is = 1;
  while(is != 0) {
    is = 0;
    for(i=1; i<subSpaceSize; ++i) {
      if(eigVal[i] < eigVal[i-1] ) {
        is = 1;
        mySwap( eigVal[i-1], eigVal[i] );
        for(j=0; j<subSpaceSize; ++j) {
          mySwap( xx[i][j], xx[i-1][j] );
        }
      }
    }
  }
}

void
getJacobi(DComplex *kappa, DComplex *mu, GenFullSquareMatrix<DComplex> &xx,
          DComplex *eigVal, int nsmax, int subSpaceSize, double tolJac)
{
  fprintf(stderr, " *** WARNING: getJacobi(...) not implemented for DComplex type \n");
}

template<>
void
GenSubDomain<DComplex>::getSRMult(DComplex *lvec, DComplex *interfvec, int nRBM, double *locRBMs, DComplex *alpha)
{
  fprintf(stderr, " *** WARNING: GenSubDomain<DComplex>::getSRMult(...) not implemented \n");
}

template<>
void
GenSubDomain<double>::getSRMult(double *lvec, double *interfvec, int nRBM, double *locRBMs, double *alpha)
{
 double *localvec = (double *) dbg_alloca(sizeof(double)*localLen());

 // Add the interface vector (interfvec) contribution to localvec
 int iDof;
 for(iDof = 0; iDof <localLen(); ++iDof)
   localvec[iDof] = lvec[iDof];

 for(iDof = 0; iDof < scomm->lenT(SComm::std); ++iDof)
   localvec[scomm->boundDofT(SComm::std,iDof)] += interfvec[iDof];

 Tgemv('T',localLen(), nRBM, -1.0, locRBMs,
       localLen(),localvec, 1, 0.0, alpha, 1);
}

template<>
void
GenSubDomain<DComplex>::precondGrbm()
{
  fprintf(stderr, " *** WARNING: GenSubDomain<DComplex>::precondGrbm() not implemented \n");
}

template<>
void
GenSubDomain<double>::precondGrbm()
{
 // WARNING: this routine has to be modified to test any new FETI-DP
 //          coarse grid ideas, if this return is uncommented, then
 //          the augmented coarse grid used the geometric rbms to enforce
 //          G^t Bu=0
 // return;
 // Just checking

 // 1. Augmenting Kcc by the Geometric Gs
 if(nGrbm == 0)
   return;

 int *boundDofs = scomm->boundDofsT(SComm::std);
 int iRBM, iDof;
 int iLen = scomm->lenT(SComm::std);
 int locLen = localRLen();

 //double *v = (double *) dbg_alloca(iLen*sizeof(double));
 double *y = (double *) dbg_alloca(nGrbm*locLen*sizeof(double));

// 2. Preconditioned Gs
/*
 for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
   for(iDof = 0; iDof < iLen; ++iDof)
   v[iDof]=interfaceRBMs[iRBM*iLen+iDof];
   multKbb(v,interfaceRBMs+iRBM*iLen);
 }
 return;
*/

/*
// 3. Try multipling by Kbb's diagonal values?
 for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
   for(iDof = 0; iDof < iLen; ++iDof)
      v[iDof]=interfaceRBMs[iRBM*iLen+iDof];
   multDiagKbb(v,interfaceRBMs+iRBM*iLen);
 }
 return;
*/

 int j;

// 4. just the local solver eigen vectors
/*
 for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
   for(iDof=0; iDof<locLen; ++iDof)
     y[iDof+iRBM*locLen]=0.0;
   for(iDof=0; iDof<iLen; ++iDof)
     y[boundDofs[iDof]+iRBM*locLen] += interfaceRBMs[iRBM*iLen+iDof];
 }
 for(j=0; j<5; ++j) {
   for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
     if(Krr) Krr->reSolve(y+iRBM*locLen);
   }
   // orthonormalize interfaceRBMs
   ortho(y, nGrbm, locLen);
 }

 for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
   for(iDof=0; iDof<iLen; ++iDof)
     interfaceRBMs[iRBM*iLen+iDof] = y[boundDofs[iDof]+iRBM*locLen];
 }
*/

 // Michel's new idea
 // 5. generalized eigen value problem
 // orthogonalize the RBMS
  // ortho(interfaceRBMs,nGrbm,iLen);
 // copy RBMS into some vectors
  //double *origR = (double *) dbg_alloca(nGrbm*locLen*sizeof(double));
  //for(iRBM = 0; iRBM < nGrbm; ++iRBM)
  // for(iDof=0; iDof<iLen; ++iDof)
  //   origR[iDof+iRBM*iLen] = interfaceRBMs[iRBM*iLen+iDof];

 // within the iteration, orthogonalize wrt original RBMS

 // Order 1: Krr^-1, Preconditioner
 double *mu    = (double *) dbg_alloca(sizeof(double)*nGrbm*(nGrbm+1)/2);
 double *kappa = (double *) dbg_alloca(sizeof(double)*nGrbm*(nGrbm+1)/2);
 double *x = (double *) dbg_alloca(sizeof(double)*iLen*nGrbm);

/*
 for(iRBM=nGrbm/2; iRBM<nGrbm; ++iRBM)
 for(iDof=0; iDof<iLen; ++iDof) {
   interfaceRBMs[iRBM*iLen+iDof]=(iDof+iRBM)%nGrbm;
 }
*/

 int numIterations=6;
 int i;
 GenFullSquareMatrix<double> xx(nGrbm);
 double *eval = (double *) dbg_alloca(sizeof(double)*nGrbm);

 for(int jj=0; jj<numIterations; ++jj) {
 int count=0;
 for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
   for(iDof=0; iDof<locLen; ++iDof)
     y[iDof+iRBM*locLen]=0.0;
   for(iDof=0; iDof<iLen; ++iDof) {
     y[boundDofs[iDof]+iRBM*locLen] += interfaceRBMs[iRBM*iLen+iDof];
//     y[boundDofs[iDof]+iRBM*locLen] +=
//              interfaceRBMs[iRBM*iLen+iDof]/scaling[iDof];
     x[iRBM*iLen+iDof] = interfaceRBMs[iRBM*iLen+iDof];
   }
 }

   if(Krr) Krr->reSolve(nGrbm, y);

 for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
   for(iDof=0; iDof<iLen; ++iDof)
     interfaceRBMs[iRBM*iLen+iDof] =
                           y[boundDofs[iDof]+iRBM*locLen];
////                       y[boundDofs[iDof]+iRBM*locLen]/scaling[iDof];

   for(j=0; j<=iRBM; ++j)
     { GenStackVector<double> Z(x+iRBM*iLen,iLen);
       GenStackVector<double> Q(interfaceRBMs+j*iLen,iLen);
       kappa[count++] = Z*Q; }
  }
  count = 0;
  for(iRBM = 0; iRBM < nGrbm; ++iRBM) {

   for(iDof = 0; iDof < iLen; ++iDof)
     x[iRBM*iLen+iDof] = interfaceRBMs[iRBM*iLen+iDof];
   multKbb(x+iRBM*iLen,interfaceRBMs+iRBM*iLen);

   for(j=0; j<=iRBM; ++j)
      { GenStackVector<double> Z(x+iRBM*iLen,iLen);
        GenStackVector<double> Q(interfaceRBMs+j*iLen,iLen);
        mu[count++] = Z*Q;
      }
  
 }
 for(iRBM = 0; iRBM < nGrbm; ++iRBM)
   for(iDof = 0; iDof < iLen; ++iDof)
     x[iRBM*iLen+iDof] = interfaceRBMs[iRBM*iLen+iDof];
 //orthoWithOrigR(origR,interfaceRBMs, nGrbm, iLen);
 //ortho(interfaceRBMs,nGrbm,iLen);
 getJacobi(kappa, mu, xx, eval, 20, nGrbm, 1e-4);

 for(j = 0; j < nGrbm; ++j) {
   GenStackVector<double> Q(interfaceRBMs+j*iLen,iLen);
   Q.zero();
   for(i = 0; i < nGrbm; ++i) {
     GenStackVector<double> Z(x+i*iLen,iLen);
     Q.linAdd(xx[j][i], Z);
   }
 }
 }

/*
 // remultiply by scaling and then Kbb
   for(iRBM = 0; iRBM < nGrbm; ++iRBM) {

     for(iDof = 0; iDof < iLen; ++iDof)
       x[iRBM*iLen+iDof] = interfaceRBMs[iRBM*iLen+iDof];
       //x[iRBM*iLen+iDof] = interfaceRBMs[iRBM*iLen+iDof]*scaling[iDof];
     multKbb(x+iRBM*iLen,interfaceRBMs+iRBM*iLen);
   }
*/
// add Krr^-1 here
 for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
   for(iDof=0; iDof<locLen; ++iDof)
     y[iDof+iRBM*locLen]=0.0;
   for(iDof=0; iDof<iLen; ++iDof) {
     y[boundDofs[iDof]+iRBM*locLen] += interfaceRBMs[iRBM*iLen+iDof];
   }

   Krr->reSolve(y+iRBM*locLen);

   for(iDof=0; iDof<iLen; ++iDof)
     interfaceRBMs[iRBM*iLen+iDof] =
                           y[boundDofs[iDof]+iRBM*locLen];
 }
}

template<> double GenSubDomain<double>::Bcx(int i) { return (bcx) ? bcx[i] : bcxC[i].real(); }
template<> DComplex GenSubDomain<DComplex>::Bcx(int i) { return (bcxC) ? bcxC[i] : DComplex(bcx[i],0.0); }

