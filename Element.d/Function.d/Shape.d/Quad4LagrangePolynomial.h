#ifndef _QUAD4LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _QUAD4LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Quad4LagrangePolynomialShapeFunction : public VectorValuedFunction<2,4,Scalar,0,0,double>
{
  public:
    Quad4LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,4,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate

      // outputs:
      // Shape functions for four-noded quad element

      Eigen::Matrix<Scalar,4,1> y;
      const Scalar &xi = q[0], &eta = q[1];
      y[0] = 1/4.*(1-xi)*(1-eta);
      y[1] = 1/4.*(1-xi)*(1+eta);
      y[2] = 1/4.*(1+xi)*(1+eta);
      y[3] = 1/4.*(1+xi)*(1-eta);

      return y;
    }
};

#endif
