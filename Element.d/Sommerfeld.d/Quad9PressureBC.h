#ifndef _QUAD9PRESSUREBC_H_
#define _QUAD9PRESSUREBC_H_
#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)

#include <Element.d/Function.d/ExternalForce.d/SurfacePressureForceFunction.h>
#include <Element.d/Function.d/Shape.d/Quad9LagrangePolynomial.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Sommerfeld.d/PressureElement.h>

template <typename S>
using Quad9LagrangePolynomialSurfacePressureForceFunction = SurfacePressureForceFunction<S, Quad9LagrangePolynomialShapeFunction, GaussLegendre2d>;

class Quad9PressureBC : public PressureElement<Quad9LagrangePolynomialSurfacePressureForceFunction>
{
  public:
    Quad9PressureBC(int* _nn, PressureBCond* _pbc); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,36,1>&, Eigen::Array<int,2,1>&);
};

#endif
#endif
