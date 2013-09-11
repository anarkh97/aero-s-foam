#ifndef _TRIANGLE10PRESSUREBC_H_
#define _TRIANGLE10PRESSUREBC_H_
#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)

#include <Element.d/Function.d/ExternalForce.d/SurfacePressureForceFunction.h>
#include <Element.d/Function.d/Shape.d/Tri10LagrangePolynomial.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Sommerfeld.d/PressureElement.h>

template <typename S>
using Tri10LagrangePolynomialSurfacePressureForceFunction = SurfacePressureForceFunction<S, Tri10LagrangePolynomialShapeFunction,
                                                                                         TriangleQuadratureRule<double,Eigen::Vector2d> >;

class Triangle10PressureBC : public PressureElement<Tri10LagrangePolynomialSurfacePressureForceFunction>
{
  public:
    Triangle10PressureBC(int* _nn, PressureBCond* _pbc); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,39,1>&, Eigen::Array<int,2,1>&);
};

#endif
#endif
