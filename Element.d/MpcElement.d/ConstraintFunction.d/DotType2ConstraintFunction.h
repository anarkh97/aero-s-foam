#ifndef _DOTTYPE2CONSTRAINTFUNCTION_H_
#define _DOTTYPE2CONSTRAINTFUNCTION_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/ConstraintFunction.h>

template<typename Scalar>
class DotType2ConstraintFunction : public RheonomicConstraintFunction<9,Scalar,6,0,double>
{
   Eigen::Matrix<double,3,1> a0, b0hat;

  public:
    DotType2ConstraintFunction(const Eigen::Array<double,6,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      Eigen::Matrix<double,3,1> b0;
      a0 << sconst(0), sconst(1), sconst(2);
      b0 << sconst(3), sconst(4), sconst(5);
      b0hat = b0.normalized();
    }

    Scalar operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar) const
    {
      // inputs:
      // q[0] = x translation of node 1
      // q[1] = y translation of node 1
      // q[2] = z translation of node 1
      // q[3] = 1st component of axis-angle rotation vector of node 1
      // q[4] = 2nd component of axis-angle rotation vector of node 1
      // q[5] = 3rd component of axis-angle rotation vector of node 1
      // q[6] = x translation of node 2
      // q[7] = y translation of node 2
      // q[8] = z translation of node 2

      // return value:
      // scalar projection of a onto b, where a = a0+u (i.e. translated), and b = R*b0 (i.e. rotated)
      // see: http://en.wikipedia.org/wiki/Scalar_projection

      // typical usage: a0 is the vector from node 1 to node 2 in the undeformed configuration
      //                b0 is a selected axis of the body attached frame of node 1 in some specified configuration

      Eigen::Matrix<Scalar,3,1> u1 = q.template segment<3>(0);
      Eigen::Quaternion<Scalar> q1;
      q1.setFromOneVector(q.template segment<3>(3));
      Eigen::Matrix<Scalar,3,1> u2 = q.template segment<3>(6);

      return (a0.template cast<Scalar>() + u2 - u1).dot(q1.toRotationMatrix()*b0hat.template cast<Scalar>());
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<> template<>
int
ConstraintJacobian<double,DotType2ConstraintFunction>
::operator() (const Eigen::Matrix<double,9,1>& q, Eigen::Matrix<double,9,1>& J) const;

template<> template<>
int
SacadoReverseJacobian<ConstraintJacobian<double,DotType2ConstraintFunction> >
::operator() (const Eigen::Matrix<double,9,1>& q, Eigen::Matrix<double,9,9>& H) const;

#endif
