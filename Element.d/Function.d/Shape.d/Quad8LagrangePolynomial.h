#ifndef _QUAD8LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _QUAD8LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Quad8LagrangePolynomialShapeFunction : public VectorValuedFunction<2,8,Scalar,0,0,double>
{
  public:
    Quad8LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,8,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate

      // outputs:
      // Shape functions for eight-noded quad element

      Eigen::Matrix<Scalar,8,1> y;
      const Scalar &xi = q[0], &eta = q[1];
      Scalar etam = 0.25*(1.0-eta);
      Scalar etap = 0.25*(1.0+eta);
      Scalar xim = 0.25*(1.0-xi);
      Scalar xip = 0.25*(1.0+xi);

      y(0) = 4.0*etam*xim*(-xi-eta-1.0);
      y(1) = 32.0*xim*etam*etap;
      y(2) = 4.0*etap*xim*(-xi+eta-1.0);
      y(3) = 32.0*xim*xip*etap;
      y(4) = 4.0*xip*etap*(xi+eta-1.0);
      y(5) = 32.0*xip*etap*etam;
      y(6) = 4.0*xip*etam*(xi-eta-1.0);
      y(7) = 32.0*xim*xip*etam;

      return y;
    }
};

#endif
