#ifndef _QUAD9LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _QUAD9LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Quad9LagrangePolynomialShapeFunction : public VectorValuedFunction<2,9,Scalar,0,0,double>
{
  public:
    Quad9LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,9,1> operator() (const Eigen::Matrix<Scalar,2,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate

      // outputs:
      // Shape functions for nine-noded quad element

      Eigen::Matrix<Scalar,9,1> y;
      const Scalar &xi = q[0], &eta = q[1];

      y[0] = 0.25*(xi-1.0)*(eta-1.0)*xi*eta;
      y[1] = 0.25*(xi+1.0)*(eta-1.0)*xi*eta;
      y[2] = 0.25*(xi+1.0)*(eta+1.0)*xi*eta;
      y[3] = 0.25*(xi-1.0)*(eta+1.0)*xi*eta;
      y[4] = 0.50*(1.0-xi*xi)*eta*(eta-1.0);
      y[5] = 0.50*(1.0-eta*eta)*xi*(xi+1.0);
      y[6] = 0.50*(1.0-xi*xi)*eta*(eta+1.0);
      y[7] = 0.50*(1.0-eta*eta)*xi*(xi-1.0);
      y[8] = (1.0-xi*xi)*(1.0-eta*eta);

      return y;
    }
};

#endif
