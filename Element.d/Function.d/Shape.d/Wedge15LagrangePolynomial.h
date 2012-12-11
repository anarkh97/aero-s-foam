#ifndef _WEDGE15LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _WEDGE15LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Wedge15LagrangePolynomialShapeFunction : public VectorValuedFunction<3,15,Scalar,0,0,double>
{
  public:
    Wedge15LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,15,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate
      // q[2] = z local coordinate

      // outputs:
      // Shape functions for 15-noded wedge element

      Eigen::Matrix<Scalar,15,1> y;
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];
      y[0] = -0.5*(1-xi-eta)*(1.0-zeta)*(2.0*xi+2.0*eta+zeta);
      y[1] = 0.5*xi*(1.0-zeta)*(2.0*xi-2.0-zeta);
      y[2] = 0.5*eta*(1.0-zeta)*(2.0*eta-2.0-zeta);
      y[3] = -0.5*(1-xi-eta)*(1.0+zeta)*(2.0*xi+2.0*eta-zeta);
      y[4] = 0.5*xi*(1.0+zeta)*(2.0*xi-2.0+zeta);
      y[5] = 0.5*eta*(1.0+zeta)*(2.0*eta-2.0+zeta);
      y[6] = 2.0*xi*(1-xi-eta)*(1.0-zeta);
      y[7] = 2.0*xi*eta*(1.0-zeta);
      y[8] = 2.0*eta*(1-xi-eta)*(1.0-zeta);
      y[9] = 2.0*xi*(1-xi-eta)*(1.0+zeta);
      y[10] = 2.0*xi*eta*(1.0+zeta);
      y[11] = 2.0*eta*(1-xi-eta)*(1.0+zeta);
      y[12] = (1-xi-eta)*(1.0-zeta*zeta);
      y[13] = xi*(1.0-zeta*zeta);
      y[14] = eta*(1.0-zeta*zeta);

      return y;
    }
};

#endif
