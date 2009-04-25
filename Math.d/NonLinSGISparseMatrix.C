#include <Utils.d/Connectivity.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/SGISparseMatrix.h>
#include <Corotational.d/GeomState.h>
#include <Timers.d/GetTime.h>

template<class Scalar>
GenNonLinSGISparseMatrix<Scalar>::GenNonLinSGISparseMatrix(Connectivity * _cn,
                 DofSetArray *_dsa,
                 ConstrainedDSA *_c_dsa, int _numele,
                 Connectivity *_dofs, Rbm *rigidBodyMode) :
                 GenSGISparseMatrix<Scalar>( _cn, _dsa, _c_dsa, rigidBodyMode)
{
 numele  = _numele;

 allDofs = _dofs;

 if(this->rbm)
   rbmflg = 1;
 else
   rbmflg = 0;

 timeZero     = 0.0;
 timeFactor   = 0.0;
 timeAssemble = 0.0;
}


template<class Scalar>
void
GenNonLinSGISparseMatrix<Scalar>::reBuildGeometricRbms(GeomState *)
{
  // if(rbmflg) this->rbm->reBuildGeometricRbms(gs);
  fprintf(stderr, "WARNING: GenNonLinSGISparseMatrix<Scalar>::reBuildGeometricRbms(GeomState *) is not implemented \n");
}


template<class Scalar>
void
GenNonLinSGISparseMatrix<Scalar>::reBuild(FullSquareMatrix *kel, int, int)
{
 // First zero SGI sparse stiffness array
 timeZero -= getTime();
 this->zeroAll();
 timeZero += getTime();

 // Assemble next stiffness
 timeAssemble -= getTime();
 int iele;
 for(iele=0; iele<numele; ++iele)
   this->add(kel[iele],(*allDofs)[iele]);
 timeAssemble += getTime();

 // Factor next stiffness matrix
 timeFactor -= getTime();
 this->factor();
 timeFactor += getTime();
}

template<class Scalar>
void
GenNonLinSGISparseMatrix<Scalar>::reBuild(FullSquareMatrix *kel, 
                                          FullSquareMatrix *mel, 
                                          double delta)
{
 // First zero SGI sparse stiffness array
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

 // Factor next stiffness matrix
 timeFactor -= getTime();
 this->factor();
 timeFactor += getTime();
}

