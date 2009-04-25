#ifndef NONLINBLOCKSKY_H_
#define NONLINBLOCKSKY_H_

#include <Utils.d/MyComplex.h>

class Connectivity;
class GeomState;
class DofSetArray;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenBlockSky;

template<class Scalar>
class GenNonLinBlockSky : public GenBlockSky<Scalar> {

 int numele;            // number of elements
 Connectivity *allDofs;	
 int rbmflg; // rbmflg = 0 - tolerance
             // rbmflg = 1 - geometric

public:
  void reBuild(FullSquareMatrix *kel, int iter=0, int step = 1);
  void reBuild(FullSquareMatrix *kel, FullSquareMatrix *mel, double delta);

  // Constructor
  GenNonLinBlockSky(Connectivity *, DofSetArray *,
                    double trbm, int numele, Connectivity*);
};

typedef GenNonLinBlockSky<double> NonLinBlockSky;
typedef GenNonLinBlockSky<DComplex> ComplexNonLinBlockSky;

#ifdef _TEMPLATE_FIX_
#include <Math.d/Skyline.d/NonLinBlockSky.C>
#endif

#endif
