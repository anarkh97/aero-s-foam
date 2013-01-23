#ifndef _HEX20LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _HEX20LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Hex20LagrangePolynomialShapeFunction : public VectorValuedFunction<3,20,Scalar,0,0,double>
{
  public:
    Hex20LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,20,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate
      // q[2] = z local coordinate

      // outputs:
      // Shape functions for 20-noded brick element

      Eigen::Matrix<Scalar,20,1> y;
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];
      y(0)  = 1/8.*(1-xi)*(1-eta)*(1-zeta)*(-xi-eta-zeta-2);
      y(1)  = 1/8.*(1+xi)*(1-eta)*(1-zeta)*( xi-eta-zeta-2);
      y(2)  = 1/8.*(1+xi)*(1+eta)*(1-zeta)*( xi+eta-zeta-2);
      y(3)  = 1/8.*(1-xi)*(1+eta)*(1-zeta)*(-xi+eta-zeta-2);
      y(4)  = 1/8.*(1-xi)*(1-eta)*(1+zeta)*(-xi-eta+zeta-2);
      y(5)  = 1/8.*(1+xi)*(1-eta)*(1+zeta)*( xi-eta+zeta-2);
      y(6)  = 1/8.*(1+xi)*(1+eta)*(1+zeta)*( xi+eta+zeta-2);
      y(7)  = 1/8.*(1-xi)*(1+eta)*(1+zeta)*(-xi+eta+zeta-2);
      y(8)  = 1/4.*(1-xi*xi)*(1-eta)*(1-zeta);
      y(9)  = 1/4.*(1+xi)*(1-eta*eta)*(1-zeta);
      y(10) = 1/4.*(1-xi*xi)*(1+eta)*(1-zeta);
      y(11) = 1/4.*(1-xi)*(1-eta*eta)*(1-zeta);
      y(12) = 1/4.*(1-xi*xi)*(1-eta)*(1+zeta);
      y(13) = 1/4.*(1+xi)*(1-eta*eta)*(1+zeta);
      y(14) = 1/4.*(1-xi*xi)*(1+eta)*(1+zeta);
      y(15) = 1/4.*(1-xi)*(1-eta*eta)*(1+zeta);
      y(16) = 1/4.*(1-xi)*(1-eta)*(1-zeta*zeta);
      y(17) = 1/4.*(1+xi)*(1-eta)*(1-zeta*zeta);
      y(18) = 1/4.*(1+xi)*(1+eta)*(1-zeta*zeta);
      y(19) = 1/4.*(1-xi)*(1+eta)*(1-zeta*zeta);

      return y;
    }
};

#endif
