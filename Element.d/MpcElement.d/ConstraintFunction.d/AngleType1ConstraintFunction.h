#ifndef _ANGLETYPE1CONSTRAINTFUNCTION_H_
#define _ANGLETYPE1CONSTRAINTFUNCTION_H_

#include <cmath>
#include <Element.d/MpcElement.d/ConstraintFunction.d/ConstraintFunction.h>

template<typename Scalar>
class AngleType1ConstraintFunction : public RheonomicConstraintFunction<6,Scalar,7,0,double>
{
   // a0 and b0 are vectors which can be interpreted as one selected axis of each of the body-attached frames attached 
   // to two nodes in some specified configuration.
   // note that the specified configuration can be either the undeformed configuration (for total lagrangian),
   // or the reference configuration at t_n (for updated lagrangian), or the current configuration at t_{n+1-alphaf}^k
   // (where k is the newton iteration) for eulerian.
   Eigen::Matrix<double,3,1> a0hat, b0hat;
   double theta0;

  public:
    AngleType1ConstraintFunction(const Eigen::Array<double,7,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      Eigen::Matrix<double,3,1> a0,b0;
      a0 << sconst(0), sconst(1), sconst(2);
      b0 << sconst(3), sconst(4), sconst(5);
      theta0 = sconst(6);
      a0hat = a0.normalized();
      b0hat = b0.normalized();
    }

    Scalar operator() (const Eigen::Matrix<Scalar,6,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = 1st component of axis/angle rotation vector of node 1
      // q[1] = 2nd component of axis/angle rotation vector of node 1
      // q[2] = 3rd component of axis/angle rotation vector of node 1
      // q[3] = 1st component of axis/angle rotation vector of node 2
      // q[4] = 2nd component of axis/angle rotation vector of node 2
      // q[5] = 3rd component of axis/angle rotation vector of node 2
      // note that the rotation vector can either be the total increment from the undeformed configuration (for total
      // total lagrangian), or the increment from the reference configuration at t_n (for updated lagrangian), or the 
      // spin (i.e. zero/infintesimal increment) from the current configuration at t_{n+1-alphaf}^k (where k is the
      // newton iteration) for eulerian.

      // return value:
      // -theta0 plus the angle between the rotated axes 
      using std::acos;

      Eigen::Quaternion<Scalar> q1;
      q1.setFromOneVector(q.template segment<3>(0));

      Eigen::Quaternion<Scalar> q2;
      q2.setFromOneVector(q.template segment<3>(3));

      return -theta0 + acos((q1.toRotationMatrix()*a0hat.template cast<Scalar>()).dot
                            (q2.toRotationMatrix()*b0hat.template cast<Scalar>()));
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<> template<>
int
ConstraintJacobian<double,AngleType1ConstraintFunction>
::operator() (const Eigen::Matrix<double,6,1>& q, Eigen::Matrix<double,6,1>& J) const;

template<> template<>
int
SacadoReverseJacobian<ConstraintJacobian<double,AngleType1ConstraintFunction> >
::operator() (const Eigen::Matrix<double,6,1>& q, Eigen::Matrix<double,6,6>& H) const;

#endif
