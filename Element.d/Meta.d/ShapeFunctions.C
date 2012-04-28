#include <Element.d/Meta.d/ShapeFunctions.h>
#include <iostream>

// shape function online resources:
//   http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/Home.html
//   http://www.kxcad.net/ansys/ANSYS/ansyshelp/thy_shp.html
//   http://web.mit.edu/calculix_v1.6/CalculiX/ccx_1.6/doc/ccx/node11.html
//   http://www.softeng.rl.ac.uk/st/projects/felib4/Software/library

#ifdef USE_EIGEN3

template<>
template<typename T>
int
ShapeFunctions<Hexahedron,8,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,3,1>& x, Eigen::Matrix<T,8,1> *y) const
{
  // Shape functions for eight-noded brick element
  const T &xi = x[0], &eta = x[1], &zeta = x[2];
  (*y)[0] = 1/8.*(1-xi)*(1-eta)*(1-zeta);
  (*y)[1] = 1/8.*(1+xi)*(1-eta)*(1-zeta);
  (*y)[2] = 1/8.*(1+xi)*(1+eta)*(1-zeta);
  (*y)[3] = 1/8.*(1-xi)*(1+eta)*(1-zeta);
  (*y)[4] = 1/8.*(1-xi)*(1-eta)*(1+zeta);
  (*y)[5] = 1/8.*(1+xi)*(1-eta)*(1+zeta);
  (*y)[6] = 1/8.*(1+xi)*(1+eta)*(1+zeta);
  (*y)[7] = 1/8.*(1-xi)*(1+eta)*(1+zeta);

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<Hexahedron,20,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,3,1>& x, Eigen::Matrix<T,20,1> *y) const
{
  // Shape functions for twenty-noded brick element
  const T &xi = x[0], &eta = x[1], &zeta = x[2];
  (*y)(0)  = 1/8.*(1-xi)*(1-eta)*(1-zeta)*(-xi-eta-zeta-2);
  (*y)(1)  = 1/8.*(1+xi)*(1-eta)*(1-zeta)*( xi-eta-zeta-2);
  (*y)(2)  = 1/8.*(1+xi)*(1+eta)*(1-zeta)*( xi+eta-zeta-2);
  (*y)(3)  = 1/8.*(1-xi)*(1+eta)*(1-zeta)*(-xi+eta-zeta-2);
  (*y)(4)  = 1/8.*(1-xi)*(1-eta)*(1+zeta)*(-xi-eta+zeta-2);
  (*y)(5)  = 1/8.*(1+xi)*(1-eta)*(1+zeta)*( xi-eta+zeta-2);
  (*y)(6)  = 1/8.*(1+xi)*(1+eta)*(1+zeta)*( xi+eta+zeta-2);
  (*y)(7)  = 1/8.*(1-xi)*(1+eta)*(1+zeta)*(-xi+eta+zeta-2);
  (*y)(8)  = 1/4.*(1-xi*xi)*(1-eta)*(1-zeta);
  (*y)(9)  = 1/4.*(1+xi)*(1-eta*eta)*(1-zeta);
  (*y)(10) = 1/4.*(1-xi*xi)*(1+eta)*(1-zeta);
  (*y)(11) = 1/4.*(1-xi)*(1-eta*eta)*(1-zeta);
  (*y)(12) = 1/4.*(1-xi*xi)*(1-eta)*(1+zeta);
  (*y)(13) = 1/4.*(1+xi)*(1-eta*eta)*(1+zeta);
  (*y)(14) = 1/4.*(1-xi*xi)*(1+eta)*(1+zeta);
  (*y)(15) = 1/4.*(1-xi)*(1-eta*eta)*(1+zeta);
  (*y)(16) = 1/4.*(1-xi)*(1-eta)*(1-zeta*zeta);
  (*y)(17) = 1/4.*(1+xi)*(1-eta)*(1-zeta*zeta);
  (*y)(18) = 1/4.*(1+xi)*(1+eta)*(1-zeta*zeta);
  (*y)(19) = 1/4.*(1-xi)*(1+eta)*(1-zeta*zeta);

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<Hexahedron,32,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,3,1>& x, Eigen::Matrix<T,32,1> *y) const
{
  // Shape functions for thirty-two-noded brick element
  const T &xi = x[0], &eta = x[1], &zeta = x[2];
  (*y)[0 ] = 1/64.*(1-xi)*(1-eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
  (*y)[1 ] = 1/64.*(1+xi)*(1-eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
  (*y)[2 ] = 1/64.*(1+xi)*(1+eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19);
  (*y)[3 ] = 1/64.*(1-xi)*(1+eta)*(1-zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19); 
  (*y)[4 ] = 1/64.*(1-xi)*(1-eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19); 
  (*y)[5 ] = 1/64.*(1+xi)*(1-eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19); 
  (*y)[6 ] = 1/64.*(1+xi)*(1+eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19); 
  (*y)[7 ] = 1/64.*(1-xi)*(1+eta)*(1+zeta)*(9*(xi*xi+eta*eta+zeta*zeta)-19); 
  (*y)[8 ] = 9/64.*(1-xi*xi)*(1-3*xi)*(1-eta)*(1-zeta);
  (*y)[9 ] = 9/64.*(1-xi*xi)*(1+3*xi)*(1-eta)*(1-zeta);
  (*y)[10] = 9/64.*(1-eta*eta)*(1-3*eta)*(1+xi)*(1-zeta);
  (*y)[11] = 9/64.*(1-eta*eta)*(1+3*eta)*(1+xi)*(1-zeta);
  (*y)[12] = 9/64.*(1-xi*xi)*(1+3*xi)*(1+eta)*(1-zeta);
  (*y)[13] = 9/64.*(1-xi*xi)*(1-3*xi)*(1+eta)*(1-zeta);
  (*y)[14] = 9/64.*(1-eta*eta)*(1+3*eta)*(1-xi)*(1-zeta);
  (*y)[15] = 9/64.*(1-eta*eta)*(1-3*eta)*(1-xi)*(1-zeta);
  (*y)[16] = 9/64.*(1-xi*xi)*(1-3*xi)*(1-eta)*(1+zeta);
  (*y)[17] = 9/64.*(1-xi*xi)*(1+3*xi)*(1-eta)*(1+zeta);
  (*y)[18] = 9/64.*(1-eta*eta)*(1-3*eta)*(1+xi)*(1+zeta);
  (*y)[19] = 9/64.*(1-eta*eta)*(1+3*eta)*(1+xi)*(1+zeta);
  (*y)[20] = 9/64.*(1-xi*xi)*(1+3*xi)*(1+eta)*(1+zeta);
  (*y)[21] = 9/64.*(1-xi*xi)*(1-3*xi)*(1+eta)*(1+zeta);
  (*y)[22] = 9/64.*(1-eta*eta)*(1+3*eta)*(1-xi)*(1+zeta);
  (*y)[23] = 9/64.*(1-eta*eta)*(1-3*eta)*(1-xi)*(1+zeta);
  (*y)[24] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1-xi)*(1-eta);
  (*y)[25] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1-xi)*(1-eta);
  (*y)[26] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1+xi)*(1-eta);
  (*y)[27] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1+xi)*(1-eta); 
  (*y)[28] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1+xi)*(1+eta);
  (*y)[29] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1+xi)*(1+eta);
  (*y)[30] = 9/64.*(1-zeta*zeta)*(1-3*zeta)*(1-xi)*(1+eta);
  (*y)[31] = 9/64.*(1-zeta*zeta)*(1+3*zeta)*(1-xi)*(1+eta);

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,2,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,2,1> *y) const
{
  // Shape functions for two-noded line element
  const T &xi = x[0];
  (*y)[0] = 1/2.*(1-xi);
  (*y)[1] = 1/2.*(1+xi);

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,3,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,3,1> *y) const
{
  // Shape functions for three-noded line element
  const T &xi = x[0];
  (*y)[0] = 1/2.*xi*(xi-1);
  (*y)[1] = 1-xi*xi;
  (*y)[2] = 1/2.*xi*(xi+1);

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<Tetrahedron,4,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,3,1>& x, Eigen::Matrix<T,4,1> *y) const
{
  // Shape functions for four-noded tetrahedral element
  const T &xi = x[0], &eta = x[1], &zeta = x[2];
  (*y)[0] = 1-xi-eta-zeta;
  (*y)[1] = xi;
  (*y)[2] = eta;
  (*y)[3] = zeta;

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<Tetrahedron,10,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,3,1>& x, Eigen::Matrix<T,10,1> *y) const
{
  // Shape functions for ten-noded tetrahedral element
  const T &xi = x[0], &eta = x[1], &zeta = x[2];
  (*y)[0] = (1-xi-eta-zeta)*(2*(1-xi-eta-zeta)-1);
  (*y)[1] = xi*(2*xi-1);
  (*y)[2] = eta*(2*eta-1);
  (*y)[3] = zeta*(2*zeta-1);
  (*y)[4] = 4*xi*(1-xi-eta-zeta);
  (*y)[5] = 4*xi*eta;
  (*y)[6] = 4*eta*(1-xi-eta-zeta);
  (*y)[7] = 4*zeta*(1-xi-eta-zeta);
  (*y)[8] = 4*zeta*xi;
  (*y)[9] = 4*zeta*eta;

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<Wedge,6,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,3,1>& x, Eigen::Matrix<T,6,1> *y) const
{
  // Shape functions for six-noded wedge element
  const T &xi = x[0], &eta = x[1], &zeta = x[2];
  (*y)[0] = 1/2.*(1-xi-eta)*(1-zeta);
  (*y)[1] = 1/2.*xi*(1-zeta);
  (*y)[2] = 1/2.*eta*(1-zeta);
  (*y)[3] = 1/2.*(1-xi-eta)*(1+zeta);
  (*y)[4] = 1/2.*xi*(1+zeta);
  (*y)[5] = 1/2.*eta*(1+zeta);

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<Wedge,15,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,3,1>& x, Eigen::Matrix<T,15,1> *y) const
{
  // Shape functions for fifteen-noded wedge element
  const T &xi = x[0], &eta = x[1], &zeta = x[2];
  (*y)[0] = -0.5*(1-xi-eta)*(1.0-zeta)*(2.0*xi+2.0*eta+zeta);
  (*y)[1] = 0.5*xi*(1.0-zeta)*(2.0*xi-2.0-zeta);
  (*y)[2] = 0.5*eta*(1.0-zeta)*(2.0*eta-2.0-zeta);
  (*y)[3] = -0.5*(1-xi-eta)*(1.0+zeta)*(2.0*xi+2.0*eta-zeta);
  (*y)[4] = 0.5*xi*(1.0+zeta)*(2.0*xi-2.0+zeta);
  (*y)[5] = 0.5*eta*(1.0+zeta)*(2.0*eta-2.0+zeta);
  (*y)[6] = 2.0*xi*(1-xi-eta)*(1.0-zeta);
  (*y)[7] = 2.0*xi*eta*(1.0-zeta);
  (*y)[8] = 2.0*eta*(1-xi-eta)*(1.0-zeta);
  (*y)[9] = 2.0*xi*(1-xi-eta)*(1.0+zeta);
  (*y)[10] = 2.0*xi*eta*(1.0+zeta);
  (*y)[11] = 2.0*eta*(1-xi-eta)*(1.0+zeta);
  (*y)[12] = (1-xi-eta)*(1.0-zeta*zeta);
  (*y)[13] = xi*(1.0-zeta*zeta);
  (*y)[14] = eta*(1.0-zeta*zeta);

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<Wedge,26,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,3,1>& x, Eigen::Matrix<T,26,1> *y) const
{
  // Shape functions for twenty-six-noded wedge element
  const T &xi = x[0], &eta = x[1], &zeta = x[2];
  (*y)[ 0] = (3*xi-1)*(3*xi-2)*xi*(1-zeta)/4.-9/8.*(1+zeta)*(zeta-1/3.)*(zeta-1)*xi+9/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*xi;
  (*y)[ 1] = (3*eta-1)*(3*eta-2)*eta*(1-zeta)/4.-9/8.*(1+zeta)*(zeta-1/3.)*(zeta-1)*eta+9/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*eta;
  (*y)[ 2] = (2-3*xi-3*eta)*(1-3*xi-3*eta)*(1-xi-eta)*(1-zeta)/4.-9/8.*(1+zeta)*(zeta-1/3.)*(zeta-1)*(1-xi-eta)+9/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*(1-xi-eta);
  (*y)[ 3] = (3*xi-1)*(3*xi-2)*xi*(1+zeta)/4.-9/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*xi+9/8.*(1+zeta)*(zeta+1/3.)*(zeta-1)*xi;
  (*y)[ 4] = (3*eta-1)*(3*eta-2)*eta*(1+zeta)/4.-9/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*eta+9/8.*(1+zeta)*(zeta+1/3.)*(zeta-1)*eta;
  (*y)[ 5] = (2-3*xi-3*eta)*(1-3*xi-3*eta)*(1-xi-eta)*(1+zeta)/4.-9/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*(1-xi-eta)+9/8.*(1+zeta)*(zeta+1/3.)*(zeta-1)*(1-xi-eta);
  (*y)[ 6] = 9/4.*xi*eta*(3*xi-1)*(1-zeta);
  (*y)[ 7] = 9/4.*xi*eta*(3*eta-1)*(1-zeta);
  (*y)[ 8] = 9/4.*(1-xi-eta)*eta*(3*eta-1)*(1-zeta);
  (*y)[ 9] = 9/4.*(1-xi-eta)*eta*(2-3*xi-3*eta)*(1-zeta);
  (*y)[10] = 9/4.*xi*(1-xi-eta)*(2-3*xi-3*eta)*(1-zeta);
  (*y)[11] = 9/4.*xi*(1-xi-eta)*(3*xi-1)*(1-zeta);
  (*y)[12] = 9/4.*xi*eta*(3*xi-1)*(1+zeta);
  (*y)[13] = 9/4.*xi*eta*(3*eta-1)*(1+zeta);
  (*y)[14] = 9/4.*(1-xi-eta)*eta*(3*eta-1)*(1+zeta);
  (*y)[15] = 9/4.*(1-xi-eta)*eta*(2-3*xi-3*eta)*(1+zeta);
  (*y)[16] = 9/4.*xi*(1-xi-eta)*(2-3*xi-3*eta)*(1+zeta);
  (*y)[17] = 9/4.*xi*(1-xi-eta)*(3*xi-1)*(1+zeta);
  (*y)[18] = 27/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*xi;
  (*y)[19] = -27/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*xi;
  (*y)[20] = 27/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*eta;
  (*y)[21] = -27/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*eta;
  (*y)[22] = 27/16.*(1+zeta)*(zeta-1/3.)*(zeta-1)*(1-xi-eta);
  (*y)[23] = -27/16.*(1+zeta)*(zeta+1/3.)*(zeta-1)*(1-xi-eta);
  (*y)[24] = 27/2.*xi*eta*(1-xi-eta)*(1-zeta);
  (*y)[25] = 27/2.*xi*eta*(1-xi-eta)*(1+zeta);

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<Triangle,3,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,2,1>& x, Eigen::Matrix<T,3,1> *y) const
{
  // Shape functions for three-noded triangle element
  const T &xi = x[0], &eta = x[1];
  (*y)[0] = xi;
  (*y)[1] = eta;
  (*y)[2] = 1-xi-eta;

/* from:  http://www.softeng.rl.ac.uk/st/projects/felib4/Software/library
  (*y)[0] = 1/3.*(1+2*xi);
  (*y)[1] = 1/3.*(1-xi-std::sqrt(3.)*eta);
  (*y)[2] = 1/3.*(1-xi+std::sqrt(3.)*eta);
*/

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<Quadrilateral,4,double,LagrangePolynomial>
::operator() (const Eigen::Matrix<T,2,1>& x, Eigen::Matrix<T,4,1> *y) const
{
  // Shape functions for four-noded quadrilateral element
  const T &xi = x[0], &eta = x[1];
  (*y)[0] = 1/4.*(1-xi)*(1-eta);
  (*y)[1] = 1/4.*(1-xi)*(1+eta);
  (*y)[2] = 1/4.*(1+xi)*(1+eta);
  (*y)[3] = 1/4.*(1+xi)*(1-eta);

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,3,double,HermiteBirkhoffType1>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,3,1> *y) const
{
  // Shape functions for two-noded quadratic hermite birkhoff line element
  // with function value and it's derivative at nodeA, and function value only at nodeB
  const T &xi = x[0];
  (*y)[0] = -0.25*xi*xi - 0.5*xi + 0.75;
  (*y)[1] = -0.5*xi*xi + 0.5;
  (*y)[2] = 0.25*xi*xi + 0.5*xi + 0.25;

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,4,double,HermiteBirkhoffType1>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,4,1> *y) const
{
  // compute shape functions y(x) for two-noded cubic hermite interpolation: phi(x) = sum y[i](x)*a[i]
  // with degrees of freedom (a): the function value and it's first and second derivatives at node Ai (xi = -1),
  // and function value only at node B (xi = +1)
  // note: the derivatives must be w.r.t. xi (local coordinate)

  const T &xi = x[0];
  (*y)[0] = -(0.125*xi*xi*xi + 0.375*xi*xi + 0.375*xi - 0.875);
  (*y)[1] = -0.25*xi*xi*xi - 0.75*xi*xi + 0.25*xi + 0.75;
  (*y)[2] = -0.25*xi*xi*xi - 0.25*xi*xi + 0.25*xi + 0.25;
  (*y)[3] = 0.125*xi*xi*xi + 0.375*xi*xi + 0.375*xi + 0.125;

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,4,double,HermiteBirkhoffType1LS>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,4,1> *y) const
{
  // (- 0.1*xx^2 - 0.4*xx + 0.6)*u1 + (- 0.2*xx^2 + 0.2*xx + 0.2)*v1 + (0.3*xx^2 + 0.2*xx - 0.3)*a1 + (0.1*xx^2 + 0.4*xx + 0.4)*u2
  const T &xi = x[0];
  (*y)[0] = (-0.1*xi*xi - 0.4*xi + 0.6); // y[0](-1) = -0.1 + 0.4 + 0.6 = 0.9
  (*y)[1] = (-0.2*xi*xi + 0.2*xi + 0.2);
  (*y)[2] = (0.3*xi*xi + 0.2*xi - 0.3);
  (*y)[3] = (0.1*xi*xi + 0.4*xi + 0.4); // y[3](1) = 0.1 + 0.4 + 0.4 = 0.9
                                        // y[3](-1) = 0.1 - 0.4 + 0.4 = 0.1

}
template<>
template<typename T>
int
ShapeFunctions<LineSegment,4,double,HermiteBirkhoffType1WLS>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,4,1> *y) const
{
  double w1 = 1.0, w2 = 1e6, w3 = 10, w4 = 1.0;
  const T &xx = x[0];
  (*y)[0] = (1.0*xx*((1.0*w1*(w1*w2 + 2.0*w1*w3 - 3.0*w2*w4 - 2.0*w3*w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4) + (1.0*w1*w2*(w1 + w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4) - (1.0*w1*(w1 + w4)*(w2 + w3))/(w1*w2*w3 + 4.0*w1*w2*w4 + 4.0*w1*w3*w4 + w2*w3*w4)) - 1.0*xx*xx*(- (1.0*w1*(w1*w2 + 4.0*w1*w4 + w2*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) + (1.0*w1*(4.0*w1*w4 - 1.0*w1*w2 + 3.0*w2*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) + (1.0*w1*w2*(w1 + w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4)) + (1.0*w1*(w1*w2 + 4.0*w1*w3 + 4.0*w1*w4 + 4.0*w2*w3 + 9.0*w2*w4 + 4.0*w3*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) - (1.0*w1*(4.0*w1*w4 - 1.0*w1*w2 + 3.0*w2*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) - (1.0*w1*(w1*w2 + 2.0*w1*w3 - 3.0*w2*w4 - 2.0*w3*w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4));
  (*y)[1] = (- 1.0*xx*xx*((2.0*w2*(w1*w2 + 4.0*w1*w4 + w2*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) - (1.0*w2*w2*(w1 + w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4)) - 1.0*xx*((2.0*w2*w2*(w1 + w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4) - (1.0*w2*(w1 + w4)*(w2 + w3))/(w1*w2*w3 + 4.0*w1*w2*w4 + 4.0*w1*w3*w4 + w2*w3*w4)) + (2.0*w2*(4.0*w1*w4 - 1.0*w1*w2 + 3.0*w2*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) + (1.0*w2*(w1*w2 + 2.0*w1*w3 - 3.0*w2*w4 - 2.0*w3*w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4));
  (*y)[2] = (- (2.0*w3*(4.0*w1*w4 - 1.0*w1*w2 + 3.0*w2*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) + (2.0*w3*xx*xx*(w1*w2 + 4.0*w1*w4 + w2*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) + (2.0*w2*w3*xx*(w1 + w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4));
  (*y)[3] = (1.0*xx*((1.0*w4*(w1*w2 + 2.0*w1*w3 - 3.0*w2*w4 - 2.0*w3*w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4) + (1.0*w2*w4*(w1 + w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4) + (1.0*w4*(w1 + w4)*(w2 + w3))/(w1*w2*w3 + 4.0*w1*w2*w4 + 4.0*w1*w3*w4 + w2*w3*w4)) + 1.0*xx*xx*((1.0*w4*(w1*w2 + 4.0*w1*w4 + w2*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) - (1.0*w4*(4.0*w1*w4 - 1.0*w1*w2 + 3.0*w2*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) + (1.0*w2*w4*(w1 + w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4)) + (1.0*w4*(w1*w2 + 4.0*w1*w3 + 4.0*w1*w4 + 4.0*w2*w3 + 9.0*w2*w4 + 4.0*w3*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) - (1.0*w4*(4.0*w1*w4 - 1.0*w1*w2 + 3.0*w2*w4))/(4.0*w1*w2*w3 + 16.0*w1*w2*w4 + 16.0*w1*w3*w4 + 4.0*w2*w3*w4) + (1.0*w4*(w1*w2 + 2.0*w1*w3 - 3.0*w2*w4 - 2.0*w3*w4))/(2.0*w1*w2*w3 + 8.0*w1*w2*w4 + 8.0*w1*w3*w4 + 2.0*w2*w3*w4));

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,5,double,HermiteBirkhoffType1LS>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,5,1> *y) const
{
  // Shape functions for two-noded cubic hermite-birkhoff least squares line element
  // with function falue and it's first, second and third derivatives at nodeA, and function value only at nodeB
  const T &xi = x[0];
 // no weighting
  (*y)[0] = (- 0.018867924528301886792452830188679*xi*xi*xi - 0.14150943396226415094339622641509*xi*xi - 0.39622641509433962264150943396226*xi + 0.64150943396226415094339622641509);
  (*y)[1] = (- 0.037735849056603773584905660377358*xi*xi*xi - 0.28301886792452830188679245283019*xi*xi + 0.20754716981132075471698113207547*xi + 0.28301886792452830188679245283019);
  (*y)[2] = (- 0.037735849056603773584905660377358*xi*xi*xi + 0.21698113207547169811320754716981*xi*xi + 0.20754716981132075471698113207547*xi - 0.21698113207547169811320754716981);
  (*y)[3] = (0.14150943396226415094339622641509*xi*xi*xi + 0.31132075471698113207547169811321*xi*xi - 0.028301886792452830188679245283019*xi - 0.31132075471698113207547169811321);
  (*y)[4] = (0.018867924528301886792452830188679*xi*xi*xi + 0.14150943396226415094339622641509*xi*xi + 0.39622641509433962264150943396226*xi + 0.35849056603773584905660377358491);
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,3,double,HermiteBirkhoffType2>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,3,1> *y) const
{
  // Shape functions for two-noded quadratic hermite-birkhoff line element
  // with function value and it's first derivative at nodeA, and first derivative only at nodeB
  // 1.0*u1 + (- 0.25*xx^2 + 0.5*xx + 0.75)*v1 + (0.25*xx^2 + 0.5*xx + 0.25)*v2
  const T &xi = x[0];
  (*y)[0] = 1.0;
  (*y)[1] = (-0.25*xi*xi + 0.5*xi + 0.75);
  (*y)[2] = (0.25*xi*xi + 0.5*xi + 0.25);

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,4,double,HermiteBirkhoffType2>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,4,1> *y) const
{
  // compute shape functions y(x) for two-noded cubic hermite-birkoff interpolation: phi(x) = sum y[i](x)*a[i]
  // with degrees of freedom (a): the function value and it's first and second derivatives at node A (xi = -1),
  // and first derivative only at node B (xi = +1)
  // note: the derivatives must be w.r.t. xi (local coordinate)
  // 1.0*u1 + (-0.08333*xx^3 - 0.25*xx^2 + 0.75*xx + 0.9167)*v1 + (-0.1667*xx^3 + 0.5*xx + 0.333)*a1 + (0.0833*xx^3 + 0.25*xx^2 + 0.25*xx + 0.0833)*v2

  const T &xi = x[0];
  (*y)[0] = 1.0;
  (*y)[1] = (-1/12.*xi*xi*xi - 0.25*xi*xi + 0.75*xi + 11/12.);
  (*y)[2] = (-1/6.*xi*xi*xi + 0.5*xi + 1/3.);
  (*y)[3] = (1/12.*xi*xi*xi + 0.25*xi*xi + 0.25*xi + 1/12.); // 1/12 + 0.25 + 0.25 +1/12. = 4/6 != 1

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,4,double,HermiteBirkhoffType2LS>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,4,1> *y) const
{
  const T &xx = x[0];
  (*y)[0] = 1.0; // *u1
  (*y)[1] = (-1/6.*xx*xx + 0.5*xx + 2/3.); // *v1
  (*y)[2] = (1/6.*xx*xx - 1/6.); // *a1
  (*y)[3] = (1/6.*xx*xx + 0.5*xx + 1/3.); // *v2
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,4,double,HermiteBirkhoffType2WLS>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,4,1> *y) const
{
  double w1 = 1, w2 = 1, w3 = 1.0, w4 = 1;
  const T &xx = x[0]; 
  (*y)[0] = (1.0*xx*((1.0*w1*(w2 + 2.0*w3 + 3.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4) - (1.0*w1*(w2 + w3 + w4))/(w2*w3 + 4.0*w2*w4 + w3*w4) + (1.0*w1*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) + (1.0*(w1*w2 + 4.0*w1*w3 + 9.0*w1*w4 + 4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) + 1.0*xx*xx*((1.0*w1*(w2 - 3.0*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) + (1.0*w1*(w2 + w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) - (1.0*w1*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) - (1.0*w1*(w2 + 2.0*w3 + 3.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4) + (1.0*w1*(w2 - 3.0*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4)); // *u1
  (*y)[1] = (1.0*xx*((1.0*w2*(w2 + w3 + w4))/(w2*w3 + 4.0*w2*w4 + w3*w4) - (2.0*w2*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) - 1.0*xx*xx*((2.0*w2*(w2 + w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) - (1.0*w2*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) + (1.0*w2*(w2 + 2.0*w3 + 3.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4) - (2.0*w2*(w2 - 3.0*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4)); // *v1
  (*y)[2] = ((2.0*w3*(w2 - 3.0*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) + (2.0*w3*xx*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4) + (2.0*w3*xx*xx*(w2 + w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4)); // *a1
  (*y)[3] = (1.0*xx*xx*((2.0*w4*(w2 + w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) + (1.0*w4*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) + 1.0*xx*((1.0*w4*(w2 + w3 + w4))/(w2*w3 + 4.0*w2*w4 + w3*w4) + (2.0*w4*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) + (1.0*w4*(w2 + 2.0*w3 + 3.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4) + (2.0*w4*(w2 - 3.0*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4)); // *v2
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,4,double,HermiteBirkhoffType3>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,4,1> *y) const
{
  // compute shape functions y(x) for two-noded cubic hermite-birkoff interpolation: phi(x) = sum y[i](x)*a[i]
  // with degrees of freedom (a): the function value and it's first and second derivatives at node A (xi = -1),
  // and second derivative only at node B (xi = +1)
  // note: the derivatives must be w.r.t. xi (local coordinate)
  //1.0*u1 + (1.0*xx + 1.0)*v1 + (- 0.0833*xx^3 + 0.25*xx^2 + 0.75*xx + 0.4167)*a1 + (0.0833*xx^3 + 0.25*xx^2 + 0.25*xx + 0.0833)*a2
  const T &xi = x[0];
  (*y)[0] = 1.0;
  (*y)[1] = (xi + 1.0);
  (*y)[2] = (-1/12.*xi*xi*xi + 0.25*xi*xi + 0.75*xi + 5/12.);
  (*y)[3] = (1/12.*xi*xi*xi + 0.25*xi*xi + 0.25*xi + 1/12.); // 1/12 + 0.25 + 0.25 +1/12. = 4/6 != 1

  return 0;
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,4,double,HermiteBirkhoffType3LS>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,4,1> *y) const
{
  // 1.0*u1 + (1.0*xx + 1.0)*v1 + (0.25*xx^2 + 0.5*xx + 0.25)*a1 + (0.25*xx^2 + 0.5*xx + 0.25)*a2
  const T &xx = x[0];
  (*y)[0] = 1.0; // *u1
  (*y)[1] = xx + 1.0; // *v1
  (*y)[2] = 0.25*xx*xx + 0.5*xx + 0.25; // *a1
  (*y)[3] = 0.25*xx*xx + 0.5*xx + 0.25; // *a2
}

template<>
template<typename T>
int
ShapeFunctions<LineSegment,4,double,HermiteBirkhoffType3WLS>
::operator() (const Eigen::Matrix<T,1,1>& x, Eigen::Matrix<T,4,1> *y) const
{
  double w1 = 1, w2 = 1, w3 = 1, w4 = 1;
  const T &xx = x[0];
  (*y)[0] = (xx*(w1/(2.0*w3 + 2.0*w4) - (w1*(w2 + w3 + w4))/(w2*(w3 + w4)) + (0.5*w1*(w2 + 2.0*w3 + 2.0*w4))/(w2*(w3 + w4))) - xx*xx*(w1/(2.0*w3 + 2.0*w4) - (2.0*w1)/(4.0*w3 + 4.0*w4)) + w1/(4.0*w3 + 4.0*w4) + (0.25*(w1*w2 + 4.0*w1*w3 + 4.0*w1*w4 + 4.0*w2*w3 + 4.0*w2*w4))/(w2*(w3 + w4)) - (0.5*w1*(w2 + 2.0*w3 + 2.0*w4))/(w2*(w3 + w4)));  // *u1
  (*y)[1] = ((0.5*(w2 + 2.0*w3 + 2.0*w4))/(w3 + w4) + xx*xx*(w2/(2.0*w3 + 2.0*w4) - (2.0*w2)/(4.0*w3 + 4.0*w4)) - (2.0*w2)/(4.0*w3 + 4.0*w4) - xx*((2.0*w2)/(2.0*w3 + 2.0*w4) - (w2 + w3 + w4)/(w3 + w4))); // *v1
  (*y)[2] = ((2.0*w3)/(4.0*w3 + 4.0*w4) + (2.0*w3*xx)/(2.0*w3 + 2.0*w4) + (2.0*w3*xx*xx)/(4.0*w3 + 4.0*w4)); // *a1
  (*y)[3] = ((2.0*w4)/(4.0*w3 + 4.0*w4) + (2.0*w4*xx)/(2.0*w3 + 2.0*w4) + (2.0*w4*xx*xx)/(4.0*w3 + 4.0*w4)); // *a2
}


// least squares type 2
// 1.0*u1 + (- 0.16666666666666666666666666666667*xx^2 + 0.5*xx + 0.66666666666666666666666666666667)*v1 + (0.16666666666666666666666666666667*xx^2 - 0.16666666666666666666666666666667)*a1 + (0.16666666666666666666666666666667*xx^2 + 0.5*xx + 0.33333333333333333333333333333333)*v2
// weighted version
// (1.0*xx*((1.0*w1*(w2 + 2.0*w3 + 3.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4) - (1.0*w1*(w2 + w3 + w4))/(w2*w3 + 4.0*w2*w4 + w3*w4) + (1.0*w1*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) + (1.0*(w1*w2 + 4.0*w1*w3 + 9.0*w1*w4 + 4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) + 1.0*xx^2*((1.0*w1*(w2 - 3.0*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) + (1.0*w1*(w2 + w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) - (1.0*w1*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) - (1.0*w1*(w2 + 2.0*w3 + 3.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4) + (1.0*w1*(w2 - 3.0*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4))*u1 + (1.0*xx*((1.0*w2*(w2 + w3 + w4))/(w2*w3 + 4.0*w2*w4 + w3*w4) - (2.0*w2*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) - 1.0*xx^2*((2.0*w2*(w2 + w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) - (1.0*w2*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) + (1.0*w2*(w2 + 2.0*w3 + 3.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4) - (2.0*w2*(w2 - 3.0*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4))*v1 + ((2.0*w3*(w2 - 3.0*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) + (2.0*w3*xx*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4) + (2.0*w3*xx^2*(w2 + w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4))*a1 + (1.0*xx^2*((2.0*w4*(w2 + w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4) + (1.0*w4*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) + 1.0*xx*((1.0*w4*(w2 + w3 + w4))/(w2*w3 + 4.0*w2*w4 + w3*w4) + (2.0*w4*(w2 - 1.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4)) + (1.0*w4*(w2 + 2.0*w3 + 3.0*w4))/(2.0*w2*w3 + 8.0*w2*w4 + 2.0*w3*w4) + (2.0*w4*(w2 - 3.0*w4))/(4.0*w2*w3 + 16.0*w2*w4 + 4.0*w3*w4))*v2


// least squares type 3
// 1.0*u1 + (1.0*xx + 1.0)*v1 + (0.25*xx^2 + 0.5*xx + 0.25)*a1 + (0.25*xx^2 + 0.5*xx + 0.25)*a2
// weighted version
// (1.0*xx*((1.0*w1)/(2.0*w3 + 2.0*w4) - (1.0*w1*(w2 + w3 + w4))/(w2*(w3 + w4)) + (0.5*w1*(w2 + 2.0*w3 + 2.0*w4))/(w2*(w3 + w4))) - 1.0*xx^2*((1.0*w1)/(2.0*w3 + 2.0*w4) - (2.0*w1)/(4.0*w3 + 4.0*w4)) + (1.0*w1)/(4.0*w3 + 4.0*w4) + (0.25*(w1*w2 + 4.0*w1*w3 + 4.0*w1*w4 + 4.0*w2*w3 + 4.0*w2*w4))/(w2*(w3 + w4)) - (0.5*w1*(w2 + 2.0*w3 + 2.0*w4))/(w2*(w3 + w4)))*u1 + ((0.5*(w2 + 2.0*w3 + 2.0*w4))/(w3 + w4) + 1.0*xx^2*((1.0*w2)/(2.0*w3 + 2.0*w4) - (2.0*w2)/(4.0*w3 + 4.0*w4)) - (2.0*w2)/(4.0*w3 + 4.0*w4) - 1.0*xx*((2.0*w2)/(2.0*w3 + 2.0*w4) - (1.0*(w2 + w3 + w4))/(w3 + w4)))*v1 + ((2.0*w3)/(4.0*w3 + 4.0*w4) + (2.0*w3*xx)/(2.0*w3 + 2.0*w4) + (2.0*w3*xx^2)/(4.0*w3 + 4.0*w4))*a1 + ((2.0*w4)/(4.0*w3 + 4.0*w4) + (2.0*w4*xx)/(2.0*w3 + 2.0*w4) + (2.0*w4*xx^2)/(4.0*w3 + 4.0*w4))*a2

#endif
