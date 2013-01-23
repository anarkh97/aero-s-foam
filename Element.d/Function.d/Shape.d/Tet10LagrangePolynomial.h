#ifndef _TET10LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _TET10LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Tet10LagrangePolynomialShapeFunction : public VectorValuedFunction<3,10,Scalar,0,0,double>
{
  public:
    Tet10LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,10,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate
      // q[2] = z local coordinate

      // outputs:
      // Shape functions for 10-noded tet element

      Eigen::Matrix<Scalar,10,1> y;
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];
      y[0] = (1-xi-eta-zeta)*(2*(1-xi-eta-zeta)-1);
      y[1] = xi*(2*xi-1);
      y[2] = eta*(2*eta-1);
      y[3] = zeta*(2*zeta-1);
      y[4] = 4*xi*(1-xi-eta-zeta);
      y[5] = 4*xi*eta;
      y[6] = 4*eta*(1-xi-eta-zeta);
      y[7] = 4*zeta*(1-xi-eta-zeta);
      y[8] = 4*zeta*xi;
      y[9] = 4*zeta*eta;

      return y;
    }
};

#endif
