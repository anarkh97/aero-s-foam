#ifndef _DOTTYPE2V2CONSTRAINTFUNCTION_H_
#define _DOTTYPE2V2CONSTRAINTFUNCTION_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/ConstraintFunction.h>

template<typename Scalar>
class DotType2v2ConstraintFunction : public RheonomicConstraintFunction<12,Scalar,9,0,double>
{
   Eigen::Matrix<double,3,1> a0,b0hat,c0hat;

  public:
    DotType2v2ConstraintFunction(const Eigen::Array<double,9,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      Eigen::Matrix<double,3,1> b0,c0;
      a0 << sconst(0), sconst(1), sconst(2);
      b0 << sconst(3), sconst(4), sconst(5);
      c0 << sconst(6), sconst(7), sconst(8);
      b0hat = b0.normalized();
      c0hat = c0.normalized();
    }

    Scalar operator() (const Eigen::Matrix<Scalar,12,1>& q, Scalar t) const
    {
      // q[0] = x translation of node 1
      // q[1] = y translation of node 1
      // q[2] = z translation of node 1
      // q[3] = 1st component of axis-angle rotation vector of node 1
      // q[4] = 2nd component of axis-angle rotation vector of node 1
      // q[5] = 3rd component of axis-angle rotation vector of node 1
      // q[6] = x translation of node 2
      // q[7] = y translation of node 2
      // q[8] = z translation of node 2
      // q[9] = 1st component of axis-angle rotation vector of node 2
      // q[10] = 2nd component of axis-angle rotation vector of node 2
      // q[11] = 3rd component of axis-angle rotation vector of node 2

      // return value:
      // sum of scalar projection of a onto b, and the scalar projection of a onto c
      // where a = a0+u (i.e. translated), b = R_1*b0 (i.e. rotated) and c = R_2*c0
      // see: http://en.wikipedia.org/wiki/Scalar_projection

      Eigen::Matrix<Scalar,3,1> u1 = q.template segment<3>(0);
      Eigen::Quaternion<Scalar> q1;
      q1.setFromOneVector(q.template segment<3>(3));

      Eigen::Matrix<Scalar,3,1> u2 = q.template segment<3>(6);
      Eigen::Quaternion<Scalar> q2;
      q2.setFromOneVector(q.template segment<3>(9));

      // "unbiased" alternative to dot constraint type 2: gives consistency with linear elasticity for beam
      return (a0.template cast<Scalar>() + u2 - u1).dot
             (q1.toRotationMatrix()*b0hat.template cast<Scalar>()+q2.toRotationMatrix()*c0hat.template cast<Scalar>());
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<> template<>
int
ConstraintJacobian<double,DotType2v2ConstraintFunction>
::operator() (const Eigen::Matrix<double,12,1>& q, Eigen::Matrix<double,12,1>& J) const;

template<> template<>
int
SacadoReverseJacobian<ConstraintJacobian<double,DotType2v2ConstraintFunction> >
::operator() (const Eigen::Matrix<double,12,1>& q, Eigen::Matrix<double,12,12>& H) const;

#endif
