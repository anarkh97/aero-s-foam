#ifdef USE_EIGEN3
#include <Element.d/Function.d/Constraint.d/PointVariPlaneDistanceConstraintFunction.h>

namespace Simo {

// specializing the member function template of constraint jacobian operator for point-
// plane distance constraint function with double precision scalar
template<>
Eigen::Matrix<double,1,12>
Jacobian<double,PointVariPlaneDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,12,1>& q, double t)
{
  Eigen::Matrix<double,12,1> J;

  Eigen::Vector3d x0 = sconst.segment<3>(0);                                       
  Eigen::Vector3d x1 = sconst.segment<3>(3);
  Eigen::Vector3d x2 = sconst.segment<3>(6);
  Eigen::Vector3d x3 = sconst.segment<3>(9);


  return J.transpose();
}

// specializing the member function template of constraint hessian operator for point-
// plane distance constraint function with double precision scalar
template<>
Eigen::Matrix<double,12,12>
Hessian<double,PointVariPlaneDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,12,1>& q, double t)
{
  Eigen::Matrix<double,12,12> H;
  H.setZero();

  return H;
}

} // namespace Simo

#endif
