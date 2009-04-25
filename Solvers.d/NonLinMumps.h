#ifndef NONLIN_MUMPS_H_
#define NONLIN_MUMPS_H_

#include<Utils.d/MyComplex.h>

class Connectivity;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class GeomState;
class DofSetArray;
class Rbm;
template <class Scalar> class GenMumpsSolver;

template<class Scalar>
class GenNonLinMumpsSolver : public GenMumpsSolver<Scalar> {

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
  GenNonLinMumpsSolver(Connectivity *, DofSetArray *, ConstrainedDSA *,
                         int numele, Connectivity *, Rbm *rbm = 0);
  virtual ~GenNonLinMumpsSolver() { /* nothing to delete */ };
};

typedef GenNonLinMumpsSolver<double> NonLinMumpsSolver;
typedef GenNonLinMumpsSolver<DComplex> ComplexNonLinMumpsSolver;

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/NonLinMumps.C>
#endif

#endif

