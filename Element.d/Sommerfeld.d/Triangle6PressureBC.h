#ifndef _TRIANGLE6PRESSUREBC_H_
#define _TRIANGLE6PRESSUREBC_H_
#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)

#include <Element.d/Sommerfeld.d/SurfacePressureForceFunction.h>
#include <Element.d/Function.d/Shape.d/Tri6LagrangePolynomial.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Element.d/Sommerfeld.d/PressureElement.h>

template <typename S>
using Tri6LagrangePolynomialSurfacePressureForceFunction = SurfacePressureForceFunction<S, Tri6LagrangePolynomialShapeFunction,
                                                                                        TriangleQuadratureRule<double,Eigen::Vector2d> >;

class Triangle6PressureBC : public PressureElement<Tri6LagrangePolynomialSurfacePressureForceFunction>
{
  public:
    Triangle6PressureBC(int* _nn, PressureBCond* _pbc); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,27,1>&, Eigen::Array<int,2,1>&);
};

#endif
#endif
