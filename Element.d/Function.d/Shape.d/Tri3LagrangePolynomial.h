#ifndef _TRI3LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _TRI3LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Tri3LagrangePolynomialShapeFunction : public VectorValuedFunction<2,3,Scalar,0,0,double>
{
  public:
    Tri3LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate

      // outputs:
      // Shape functions for three-noded triangle element

      Eigen::Matrix<Scalar,3,1> y;
      const Scalar &xi = q[0], &eta = q[1];
      y[0] = xi;
      y[1] = eta;
      y[2] = 1-xi-eta;

      return y;
    }
};

#endif
