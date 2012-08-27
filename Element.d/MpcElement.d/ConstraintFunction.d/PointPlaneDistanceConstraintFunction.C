#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/ConstraintFunction.d/PointPlaneDistanceConstraintFunction.h>

// specializing the member function template of constraint jacobian operator for point-
// plane distance constraint function with double precision scalar
template<> template<>
int
ConstraintJacobian<double,PointPlaneDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, Eigen::Matrix<double,3,1>& J) const
{
  Eigen::Vector3d x0 = sconst.segment<3>(0);                                       
  Eigen::Vector3d x1 = sconst.segment<3>(3);
  Eigen::Vector3d x2 = sconst.segment<3>(6);
  Eigen::Vector3d x3 = sconst.segment<3>(9);

  // unit normal to plane
  Eigen::Vector3d nhat = (x2-x1).cross(x3-x1).normalized();

  if(iconst[0]) 
    J = -nhat;
  else
    J = nhat;

  return 1;
}

// specializing the member function template of constraint hessian operator for point-
// plane distance constraint function with double precision scalar
template<> template<>
int
SacadoReverseJacobian<ConstraintJacobian<double,PointPlaneDistanceConstraintFunction> >
::operator() (const Eigen::Matrix<double,3,1>& q, Eigen::Matrix<double,3,3>& H) const
{
  H.setZero();

  return 1;
}
#endif
