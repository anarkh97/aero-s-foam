#include <Math.d/Skyline.d/BlockSky.h>
#include <Utils.d/Connectivity.h>
#include <Solvers.d/Rbm.h>
#include <Corotational.d/GeomState.h>

// Non-linear skyline matrix class

template<class Scalar>
GenNonLinBlockSky<Scalar>::GenNonLinBlockSky(Connectivity * _cn,
                 DofSetArray *_c_dsa, double _trbm,int _numele,
                 Connectivity *_dofs) :
                 GenBlockSky<Scalar>( _cn, _c_dsa, _trbm)
{
 numele  = _numele;
 allDofs = _dofs;
}

template<class Scalar>
void
GenNonLinBlockSky<Scalar>::reBuild(FullSquareMatrix *kel, int, int)
{
 // First zero Skyline stiffness array
 this->zeroAll();

 // Assemble next stiffness
 int iele;
 for(iele=0; iele<numele; ++iele)
   this->add(kel[iele],(*allDofs)[iele]);

 // Factor next stiffness
 this->factor();
}

template<class Scalar>
void
GenNonLinBlockSky<Scalar>::reBuild(FullSquareMatrix *kel, 
                                   FullSquareMatrix *mel, 
                                   double delta)
{
 // First zero Skyline stiffness array
 this->zeroAll();

 double delta2 = delta*delta;

 // Assemble nonlinear dynamic tangnent stiffness matrix
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

 // Factor next stiffness
 this->factor();
}

