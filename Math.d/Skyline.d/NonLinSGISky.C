#include <Utils.d/Connectivity.h>
#include <Solvers.d/Rbm.h>
#include <Math.d/Skyline.d/NonLinSGISky.h>
#include <Corotational.d/GeomState.h>

// Non-linear skyline matrix class

NonLinSGISky::NonLinSGISky(Connectivity * _cn,DofSetArray *_dsa,
                 ConstrainedDSA *_c_dsa, double _trbm,int _numele,
                 Connectivity *_dofs, Rbm *rbm) :
                 SGISky( _cn, _dsa, _c_dsa, _trbm, rbm)
{
 numele  = _numele;
 allDofs = _dofs;

 if(rbm)
   rbmflg = 1;
 else
   rbmflg = 0;
}

void
NonLinSGISky::reBuildGeometricRbms(GeomState *gs)
{
  // if we are using geometric rbm method, rebuild
  //if(rbmflg)
  //  rbm->reBuildGeometricRbms(gs);
}

void
NonLinSGISky::reBuild(FullSquareMatrix *kel, int, int)
{
 // First zero Skyline stiffness array
 zeroAll();

 // Assemble next stiffness
 int iele;
 for(iele=0; iele<numele; ++iele)
   add(kel[iele],(*allDofs)[iele]);

 // Factor next stiffness
 factor();

}

void
NonLinSGISky::reBuild(FullSquareMatrix *kel, FullSquareMatrix *mel, 
                         double delta)
{
 // First zero Skyline stiffness array
 zeroAll();

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

   add(kel[iele],(*allDofs)[iele]);
 }

 // Factor next stiffness
 factor();

}

