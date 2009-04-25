#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Utils.d/Connectivity.h>
#include <Solvers.d/Rbm.h>
#include <Corotational.d/GeomState.h>
#include <Timers.d/GetTime.h>

// Non-linear skyline matrix class

template<class Scalar>
GenNonLinSkyMatrix<Scalar>::GenNonLinSkyMatrix(Connectivity * cn,
                 DofSetArray *c_dsa, double _trbm,int _numele,
                 Connectivity *_dofs, Rbm *rigidBodyModes) :
                 GenSkyMatrix<Scalar>(cn, c_dsa, _trbm, rigidBodyModes)
{
 numele  = _numele;
 allDofs = _dofs;
 if(this->rbm)
   rbmflg = 1;
 else
   rbmflg = 0;

 timeZero     = 0.0;
 timeAssemble = 0.0;
 timeFactor   = 0.0;
}

template<class Scalar>
void
GenNonLinSkyMatrix<Scalar>::reBuildGeometricRbms(GeomState *gs)
{
  // if we are using geometric this->rbm method, rebuild

  timeRbm -= getTime();
  if(rbmflg)
    this->rbm->reBuildGeometricRbms(gs);
  timeRbm += getTime();
}

template<class Scalar>
void
GenNonLinSkyMatrix<Scalar>::reBuild(FullSquareMatrix *kel, int, int)
{
 // First zero Skyline stiffness array
 timeZero -= getTime();
 this->zeroAll();
 timeZero += getTime();

 // Assemble next stiffness
 timeAssemble -= getTime();
 int iele;
 for(iele=0; iele<numele; ++iele)
   this->add(kel[iele],(*allDofs)[iele]);
 timeAssemble += getTime();

 // this->Factor next stiffness
 timeFactor -= getTime();
 if(rbmflg) {
   this->Factor(this->rbm);
 }
 else {
   this->Factor();
 }
 timeFactor += getTime();

}

template<class Scalar>
void
GenNonLinSkyMatrix<Scalar>::reBuild(FullSquareMatrix *kel, 
                                    FullSquareMatrix *mel, 
                                    double delta)
{
 // First zero Skyline stiffness array
 timeZero -= getTime();
 this->zeroAll();
 timeZero += getTime();

 double delta2 = delta*delta;

 // Assemble nonlinear dynamic tangnent stiffness matrix
 timeAssemble -= getTime();
 int iele;
 for( iele=0; iele<numele; ++iele) {
   int dim = kel[iele].dim();

   int i,j;
   for(i = 0; i < dim; ++i)
     for(j = 0; j < dim; ++j) {
       double m = mel[iele][i][j];
       double k = kel[iele][i][j];

       kel[iele][i][j] = delta2*k + m;
     }

   this->add( kel[iele], (*allDofs)[iele] );
 }
 timeAssemble += getTime();

 timeFactor -= getTime();
 this->Factor();
 timeFactor += getTime();
}

// KHP: not needed currently
//      nlGeometricFactor is the hybrid method proposed by Charbel
//      to take the number of rigid body modes computed by GRBM
//      and use svbu4gm to compute this number of rbms.
//      we are now just rebuilding the geometric rigid body modes at
//      each nonlinear iteration
//
// This code is commented out on purpose in case it is needed in the future.
// NOTE: not templated yet
/*
void
SkyMatrix::nlGeometricFactor(Rbm *rigid)
{
   int i, j;

   double *w1 = (double *) dbg_alloca(5*numUncon*sizeof(double));

   double *w2 = w1 + numUncon;
   double *w3 = w2 + numUncon;
   double *w4 = w3 + numUncon;

   double *dummyZEM = w4 + numUncon;
   double *rhs = 0;

   // ... INITIALIZE NUMBER OF OPERATIONS (NOPS)
   double nops = 0.0;

   // ... INITIALIZE # of GEOMETRIC RIGID BODY MODES
   int ngrbm = rigid->numRBM();

   // ... DECLARE A POINTER TO ZEM
   double *zem = NULL;

   // ... GET THE NUMBER OF COMPONENTS
   int numComp = rigid->numComponents();

   // ... INITIALIZE # OF ZERO ENERGY MODES (ZEM)
   int nTotZem = 0;
   int *nZemPerComp = (int *) dbg_alloca(numComp * sizeof(int));

   // ... ALLOCATE MEMORY FOR S1 AND S2
   // int sizeS = (numUncon - kstop)*(numUncon - kstop);

   double *s1 =  (double *) dbg_alloca(400*sizeof(double));
   double *s2 =  (double *) dbg_alloca(400*sizeof(double));

   // ... CALL WITH FLAG = 1 FOR PARTIAL FACTORIZATION
   int flag = 1;

   // ... LOOP OVER COMPONENTS
   int n;
   for(n=0; n<numComp; ++n) {

     // Get number of this->rbm per component
     int locNgrbm = rigid->numRBM(n); // modify this function

     // Get number of dof per component
     int numDofPerComp  = rigid->numDof(n);

     // Get first dof of the component
     int firstDofOfComp = rigid->firstDof(n);

     int kstop = 4*((numDofPerComp - 10)/4);
     if(kstop < 4) kstop = 4;

     int k;
     if(firstDofOfComp != 0)
       for(k = 0; k < numDofPerComp; ++k)
          lacol[firstDofOfComp+k] -= firstDofOfComp;

     _FORTRAN(svbu4gm)(skyA,dlp+firstDofOfComp,lacol+firstDofOfComp,rhs,
      w1+firstDofOfComp, w2+firstDofOfComp,w3+firstDofOfComp,w4+firstDofOfComp,
      pivot+firstDofOfComp,TOLERANCE,numDofPerComp,flag,nops,nzem,
      dummyZEM,kstop,NULL,NULL,locNgrbm);

     if(nzem != 0) fprintf(stderr,"Found %d mechanisms\n",nzem);
     nTotZem += nzem;
     nZemPerComp[n] = nzem;

   }

   Vector *allrbms = new Vector[nTotZem+ngrbm];
   Vector v(numUncon,0.0);

   for(i=0; i<nTotZem+ngrbm; ++i)
     allrbms[i] = v;

   // ... CALL WITH FLAG = 2 TO COMPLETE FACTORIZATION
   flag = 2;

   nTotZem = 0;
   for(n=0; n<numComp; ++n) {

     // Get number of this->rbm per component
     int locNgrbm = rigid->numRBM(n); // modify this function

     // Get number of dof per component
     int numDofPerComp  = rigid->numDof(n);

     // Get first dof of the component
     int firstDofOfComp = rigid->firstDof(n);

     int kstop = 4*((numDofPerComp - 10)/4);
     if(kstop < 4) kstop = 4;

     // ... UPDATE nzem
     int ZemPlusGrbm = nZemPerComp[n] + locNgrbm;

     // ... ALLOCATE MEMORY FOR nzem RIGID BODY MODES (even if we don't care)
     zem = (ZemPlusGrbm) ? (double *) dbg_alloca(numDofPerComp*ZemPlusGrbm*sizeof(double)) : NULL;

     _FORTRAN(svbu4gm)(skyA,dlp+firstDofOfComp,lacol+firstDofOfComp,
        rhs,w1+firstDofOfComp, w2+firstDofOfComp,w3+firstDofOfComp,w4+firstDofOfComp,
         pivot+firstDofOfComp,TOLERANCE,numDofPerComp,
         flag,nops,ZemPlusGrbm,zem,kstop,s1,s2,locNgrbm);

     for(i = 0; i < ZemPlusGrbm; ++i) {
       for(j = 0 ; j < numDofPerComp; ++j)
         allrbms[i][j+firstDofOfComp] = zem[i*numDofPerComp+j];
     }

    nTotZem += nZemPerComp[n];
    int k;
    if(firstDofOfComp != 0)
       for(k = 0; k < numDofPerComp; ++k)
          lacol[firstDofOfComp+k] += firstDofOfComp;
   }

   nzem = nTotZem + ngrbm;

   this->rbm->setGrbm(allrbms); 
}
*/
