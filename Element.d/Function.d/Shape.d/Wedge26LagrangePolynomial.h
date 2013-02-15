#ifndef _WEDGE26LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _WEDGE26LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Wedge26LagrangePolynomialShapeFunction : public VectorValuedFunction<3,26,Scalar,0,0,double>
{
  public:
    Wedge26LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,26,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate
      // q[2] = z local coordinate

      // outputs:
      // Shape functions for 26-noded wedge element

      Eigen::Matrix<Scalar,26,1> y;
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];
      y[0] = 1/8.*(1-xi)*(1-eta)*(1-zeta);
      y[1] = 1/8.*(1+xi)*(1-eta)*(1-zeta);
      y[2] = 1/8.*(1+xi)*(1+eta)*(1-zeta);
      y[3] = 1/8.*(1-xi)*(1+eta)*(1-zeta);
      y[4] = 1/8.*(1-xi)*(1-eta)*(1+zeta);
      y[5] = 1/8.*(1+xi)*(1-eta)*(1+zeta);
      y[6] = 1/8.*(1+xi)*(1+eta)*(1+zeta);
      y[7] = 1/8.*(1-xi)*(1+eta)*(1+zeta);

      return y;
    }
};

#endif