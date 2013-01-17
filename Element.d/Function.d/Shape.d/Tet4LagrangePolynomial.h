#ifndef _TET4LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _TET4LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Tet4LagrangePolynomialShapeFunction : public VectorValuedFunction<3,4,Scalar,0,0,double>
{
  public:
    Tet4LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,4,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate
      // q[2] = z local coordinate

      // outputs:
      // Shape functions for four-noded tet element

      Eigen::Matrix<Scalar,4,1> y;
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];
      y[0] = 1-xi-eta-zeta;
      y[1] = xi;
      y[2] = eta;
      y[3] = zeta;

      return y;
    }
};

#endif
