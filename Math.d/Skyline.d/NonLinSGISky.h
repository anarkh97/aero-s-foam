#ifndef NONLINSGISKY_H_
#define NONLINSGISKY_H_

#include <Math.d/Skyline.d/SGISky.h>

class Connectivity;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class GeomState;
class DofSetArray;
class ConstrainedDSA;
class Rbm;

class NonLinSGISky : public SGISky {

 int numele;            // number of elements
 Connectivity *allDofs;	
 int rbmflg; // rbmflg = 0 - tolerance
             // rbmflg = 1 - geometric

public:
  void reBuild(FullSquareMatrix *kel, int iter=0, int step = 1);
  void reBuild(FullSquareMatrix *kel, FullSquareMatrix *mel, double delta);
  void reBuildGeometricRbms(GeomState *gs);

  // Constructor
  NonLinSGISky(Connectivity *, DofSetArray *, ConstrainedDSA *,
                  double trbm,int numele,Connectivity*,Rbm *rbm=0);
};

#endif
