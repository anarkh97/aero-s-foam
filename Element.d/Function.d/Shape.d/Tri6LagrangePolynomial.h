#ifndef _TRI6LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _TRI6LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Tri6LagrangePolynomialShapeFunction : public VectorValuedFunction<2,6,Scalar,0,0,double>
{
  public:
    Tri6LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,6,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate

      // outputs:
      // Shape functions for six-noded triangle element

      Eigen::Matrix<Scalar,6,1> y;
      const Scalar &xi = q[0], &eta = q[1];

      using std::sqrt;
      Scalar l1 = 1/3.*(1+2*xi);
      Scalar l2 = 1/3.*(1-xi-sqrt(3)*eta);
      Scalar l3 = 1/3.*(1-xi+sqrt(3)*eta);

      y(1) = (2*l1-1)*l1;
      y(2) = 4*l1*l2;
      y(3) = (2*l2-1)*l2;
      y(4) = 4*l2*l3;
      y(5) = (2*l3-1)*l3;
      y(6) = 4*l3*l1;

      return y;
    }
};

#endif
