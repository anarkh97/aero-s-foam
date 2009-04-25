#ifndef NONLINSKYMATRIX_H_
#define NONLINSKYMATRIX_H_

#include <Utils.d/MyComplex.h>

class Connectivity;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class GeomState;
class DofSetArray;
class DofSetArray;
class Rbm;
template<class Scalar> class GenSkyMatrix;

template<class Scalar>
class GenNonLinSkyMatrix : public GenSkyMatrix<Scalar> 
{
 int numele; // number of elements
 Connectivity *allDofs;	
 int rbmflg; // rbmflg = 0 - tolerance
             // rbmflg = 1 - geometric

 double timeZero;
 double timeAssemble;
 double timeFactor;
 double timeRbm;

public:
  void reBuild(FullSquareMatrix *kel, int iter=0, int step = 1);
  void reBuild(FullSquareMatrix *kel, FullSquareMatrix *mel, double delta);
  void reBuildGeometricRbms(GeomState *gs);

  GenNonLinSkyMatrix(Connectivity *, DofSetArray *,
                     double trbm, int numele, Connectivity*, Rbm *rbm=0);

  Connectivity *getAllDofs()  { return allDofs; }
};

typedef GenNonLinSkyMatrix<double> NonLinSkyMatrix;
typedef GenNonLinSkyMatrix<DComplex> ComplexNonLinSkyMatirx;

#ifdef _TEMPLATE_FIX_
#include <Math.d/Skyline.d/NonLinSkyMatrix.C>
#endif

#endif
