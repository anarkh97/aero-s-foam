#ifndef NONLIN_BLK_SPARSEMATRIX_H_
#define NONLIN_BLK_SPARSEMATRIX_H_

#include<Utils.d/MyComplex.h>

class Connectivity;
class GeomState;
class DofSetArray;
class DofSetArray;
class Rbm;
template <class Scalar> class GenBLKSparseMatrix;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;

template<class Scalar>
class GenNonLinBLKSparseMatrix : public GenBLKSparseMatrix<Scalar>
{
 int numele; // number of elements
 Connectivity *allDofs;	
 double tolerance;
 int rbmflg; // rbmflg = 0 - tolerance
             // rbmflg = 1 - geometric

 double timeZero;
 double timeAssemble;
 double timeFactor;

public:
  void reBuild(FullSquareMatrix *kel, int iter=0, int step=1);
  void reBuild(FullSquareMatrix *kel, FullSquareMatrix *mel, double delta);
  void reBuildGeometricRbms(GeomState *gs);
  Connectivity *getAllDofs()  { return allDofs; }

  // Constructor
  GenNonLinBLKSparseMatrix(Connectivity *, DofSetArray *, DofSetArray *,
                           int numele, Connectivity*, double tol, int spRenum, Rbm *rbm=0);
  virtual ~GenNonLinBLKSparseMatrix() { /* nothing to delete */ };
};

typedef GenNonLinBLKSparseMatrix<double> NonLinBLKSparseMatrix;
typedef GenNonLinBLKSparseMatrix<DComplex> ComplexNonLinBLKSparseMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/NonLinBLKSparseMatrix.C>
#endif


#endif
