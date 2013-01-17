#ifndef _QUAD8PRESSUREBC_H_
#define _QUAD8PRESSUREBC_H_
#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)

#include <Element.d/Sommerfeld.d/SurfacePressureForceFunction.h>
#include <Element.d/Function.d/Shape.d/Quad8LagrangePolynomial.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Sommerfeld.d/PressureElement.h>

template <typename S>
using Quad8LagrangePolynomialSurfacePressureForceFunction = SurfacePressureForceFunction<S, Quad8LagrangePolynomialShapeFunction, GaussLegendre2d>;

class Quad8PressureBC : public PressureElement<Quad8LagrangePolynomialSurfacePressureForceFunction>
{
  public:
    Quad8PressureBC(int* _nn, double _pressure); 

  protected:
    double pressure;
    void getConstants(CoordSet& cs, Eigen::Array<double,25,1>&, Eigen::Array<int,1,1>&);
};

#endif
#endif
