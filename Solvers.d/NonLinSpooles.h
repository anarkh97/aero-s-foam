#ifndef NONLIN_SPOOLES_H_
#define NONLIN_SPOOLES_H_

#include<Utils.d/MyComplex.h>

class Connectivity;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class GeomState;
class DofSetArray;
class Rbm;
template <class Scalar> class GenSpoolesSolver;

template<class Scalar>
class GenNonLinSpoolesSolver : public GenSpoolesSolver<Scalar> {

 int numele;            // number of elements
 Connectivity *allDofs;
 // double tolerance;
 // int rbmflg; // rbmflg = 0 - tolerance
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
  GenNonLinSpoolesSolver(Connectivity *, DofSetArray *, ConstrainedDSA *,
                         int numele, Connectivity *, Rbm *rbm = 0);
  virtual ~GenNonLinSpoolesSolver() { /* nothing to delete */ };
};

typedef GenNonLinSpoolesSolver<double> NonLinSpoolesSolver;
typedef GenNonLinSpoolesSolver<DComplex> ComplexNonLinSpoolesSolver;

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/NonLinSpooles.C>
#endif

#endif

