#ifndef _HEX32LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_
#define _HEX32LAGRANGEPOLYNOMIALSHAPEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>

template<typename Scalar>
class Hex32LagrangePolynomialShapeFunction : public VectorValuedFunction<3,32,Scalar,0,0,double>
{
  public:
    Hex32LagrangePolynomialShapeFunction(const Eigen::Array<double,0,1>&, const Eigen::Array<int,0,1>&)
    {}

    Eigen::Matrix<Scalar,32,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x local coordinate
      // q[1] = y local coordinate
      // q[2] = z local coordinate

      // outputs:
      // Shape functions for 32-noded brick element

      Eigen::Matrix<Scalar,32,1> y;
      const Scalar &xi = q[0], &eta = q[1], &zeta = q[2];
      y[0 ] = 1/64.*(1-xi)*(1-eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      y[1 ] = 1/64.*(1+xi)*(1-eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      y[2 ] = 1/64.*(1+xi)*(1+eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      y[3 ] = 1/64.*(1-xi)*(1+eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      y[4 ] = 1/64.*(1-xi)*(1-eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      y[5 ] = 1/64.*(1+xi)*(1-eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      y[6 ] = 1/64.*(1+xi)*(1+eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      y[7 ] = 1/64.*(1-xi)*(1+eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
      y[8 ] = 9/64.*(1-xi*xi)*(1-3*xi)*(1-eta)*(1-zeta);
      y[9 ] = 9/64.*(1-xi*xi)*(1+3*xi)*(1-eta)*(1-zeta);
      y[10] = 9/64.*(1-eta*eta)*(1-3*eta)*(1+xi)*(1-zeta);
      y[11] = 9/64.*(1-eta*eta)*(1+3*eta)*(1+xi)*(1-zeta);
      y[12] = 9/64.*(1-xi*xi)*(1+3*xi)*(1+eta)*(1-zeta);
      y[13] = 9/64.*(1-xi*xi)*(1-3*xi)*(1+eta)*(1-zeta);
      y[14] = 9/64.*(1-eta*eta)*(1+3*eta)*(1-xi)*(1-zeta);
      y[15] = 9/64.*(1-eta*eta)*(1-3*eta)*(1-xi)*(1-zeta);
      y[16] = 9/64.*(1-xi*xi)*(1-3*xi)*(1-eta)*(1+zeta);
      y[17] = 9/64.*(1-xi*xi)*(1+3*xi)*(1-eta)*(1+zeta);
      y[18] = 9/64.*(1-eta*eta)*(1-3*eta)*(1+xi)*(1+zeta);
      y[19] = 9/64.*(1-eta*eta)*(1+3*eta)*(1+xi)*(1+zeta);
      y[20] = 9/64.*(1-xi*xi)*(1+3*xi)*(1+eta)*(1+zeta);
      y[21] = 9/64.*(1-xi*xi)*(1-3*xi)*(1+eta)*(1+zeta);
      y[22] = 9/64.*(1-eta*eta)*(1+3*eta)*(1-xi)*(1+zeta);
      y[23] = 9/64.*(1-eta*eta)*(1-3*eta)*(1-xi)*(1+zeta);
      y[24] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1-xi)*(1-eta);
      y[25] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1-xi)*(1-eta);
      y[26] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1+xi)*(1-eta);
      y[27] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1+xi)*(1-eta);
      y[28] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1+xi)*(1+eta);
      y[29] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1+xi)*(1+eta);
      y[30] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1-xi)*(1+eta);
      y[31] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1-xi)*(1+eta);
    
      return y;
    }
};

#endif
