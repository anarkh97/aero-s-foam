#ifndef _QUAD12PRESSUREBC_H_
#define _QUAD12PRESSUREBC_H_
#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)

#include <Element.d/Function.d/ExternalForce.d/SurfacePressureForceFunction.h>
#include <Element.d/Function.d/Shape.d/Quad12LagrangePolynomial.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Sommerfeld.d/PressureElement.h>

template <typename S>
using Quad12LagrangePolynomialSurfacePressureForceFunction = SurfacePressureForceFunction<S, Quad12LagrangePolynomialShapeFunction, GaussLegendre2d>;

class Quad12PressureBC : public PressureElement<Quad12LagrangePolynomialSurfacePressureForceFunction>
{
  public:
    Quad12PressureBC(int* _nn, PressureBCond* _pbc); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,45,1>&, Eigen::Array<int,2,1>&);
};

#endif
#endif
