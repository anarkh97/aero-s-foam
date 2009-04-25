#include <Utils.d/Connectivity.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Corotational.d/GeomState.h>
#include <Timers.d/GetTime.h>

template<class Scalar>
GenNonLinBLKSparseMatrix<Scalar>::GenNonLinBLKSparseMatrix(Connectivity * _cn,
                 DofSetArray *_dsa,
                 DofSetArray *_c_dsa, int _numele,
                 Connectivity *_dofs, double _tolerance, int spRenum, Rbm *rigidBodyMode) :
                 GenBLKSparseMatrix<Scalar>( _cn, _dsa, _c_dsa, _tolerance, spRenum, rigidBodyMode)
{
 numele  = _numele;
 allDofs = _dofs;
 tolerance = _tolerance;

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
GenNonLinBLKSparseMatrix<Scalar>::reBuildGeometricRbms(GeomState *)
{
  // if(rbmflg) this->rbm->reBuildGeometricRbms(gs);
  fprintf(stderr, "WARNING: GenNonLinBLKSparseMatrix<Scalar>::reBuildGeometricRbms(GeomState *) is not implemented \n");
}


template<class Scalar>
void
GenNonLinBLKSparseMatrix<Scalar>::reBuild(FullSquareMatrix *kel, int, int)
{
 // First zero Block Sparse stiffness array
 timeZero -= getTime();
 this->zeroAll();
 timeZero += getTime();

 // Assemble next stiffness
 timeAssemble -= getTime();
 int iele;
 for(iele=0; iele<numele; ++iele)
   this->add(kel[iele],(*allDofs)[iele]);
 timeAssemble += getTime();

 // Factor next stiffness
 timeFactor -= getTime();
 this->factor();
 timeFactor += getTime();
}

template<class Scalar>
void
GenNonLinBLKSparseMatrix<Scalar>::reBuild(FullSquareMatrix *kel, 
                                          FullSquareMatrix *mel, 
                                          double delta)
{
 // First zero Block Sparse stiffness array
 timeZero -= getTime();
 this->zeroAll();
 timeZero += getTime();

 double delta2 = delta*delta;

 // Assemble nonlinear dynamic tangnent stiffness matrix
 timeAssemble -= getTime();
 int iele;
 for(iele=0; iele<numele; ++iele) {
   int dim = kel[iele].dim();

   int i,j;
   for(i = 0; i < dim; ++i)
     for(j = 0; j < dim; ++j) {
       double m = mel[iele][i][j];
       double k = kel[iele][i][j];

       kel[iele][i][j] = delta2*k + m;
     }

   this->add(kel[iele],(*allDofs)[iele]);
 }
 timeAssemble += getTime();

 timeFactor -= getTime();
 this->factor();
 timeFactor += getTime();
}

