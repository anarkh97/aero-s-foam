#ifndef _WEDGE6LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _WEDGE6LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Wedge6LagrangePolynomialShapeFunction : public VectorValuedFunction<3,6,Scalar,0,0,double>
{
  public:
    Wedge6LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,6,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate
      // q[2] = z local coordinate

      // outputs:
      // Shape functions for six-noded wedge element

      Eigen::Matrix<Scalar,6,1> y;
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];
      y[0] = 1/2.*(1-xi-eta)*(1-zeta);
      y[1] = 1/2.*xi*(1-zeta);
      y[2] = 1/2.*eta*(1-zeta);
      y[3] = 1/2.*(1-xi-eta)*(1+zeta);
      y[4] = 1/2.*xi*(1+zeta);
      y[5] = 1/2.*eta*(1+zeta);

      return y;
    }
};

#endif
