#ifndef NONLIN_SGISPARSEMATRIX_H_
#define NONLIN_SGISPARSEMATRIX_H_

#include<Utils.d/MyComplex.h>

class Connectivity;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class GeomState;
class DofSetArray;
class ConstrainedDSA;
class Rbm;
template <class Scalar> class GenSGISparseMatrix;

template<class Scalar>
class GenNonLinSGISparseMatrix : public GenSGISparseMatrix<Scalar> 
{
  int numele;            // number of elements
  Connectivity *allDofs;	
  int rbmflg; // rbmflg = 0 - tolerance
              // rbmflg = 1 - geometric
  double timeZero;
  double timeFactor;
  double timeAssemble;
  
public:
  void reBuild(FullSquareMatrix *kel, int iter=0, int step=1);
  void reBuild(FullSquareMatrix *kel, FullSquareMatrix *mel, double delta);
  void reBuildGeometricRbms(GeomState *gs);

  // Constructor
  GenNonLinSGISparseMatrix(Connectivity *, DofSetArray *, ConstrainedDSA *,
                           int numele, Connectivity*, Rbm *rbm=0);
  Connectivity *getAllDofs()  { return allDofs; }
};

typedef GenNonLinSGISparseMatrix<double> NonLinSGISparseMatrix;
typedef GenNonLinSGISparseMatrix<DComplex> ComplexNonLinSGISparseMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/NonLinSGISparseMatrix.C>
#endif

#endif
