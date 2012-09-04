#ifndef _DISTANCECONSTRAINTFUNCTION_H_
#define _DISTANCECONSTRAINTFUNCTION_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/ConstraintFunction.h>
#include <iostream>

template<typename Scalar>
class DistanceConstraintFunction : public RheonomicConstraintFunction<6,Scalar,4,0,double>
{
    Eigen::Matrix<double,3,1> a0; // initial vector from node 1 to node 2
    double len;                   // constant length

  public:
    DistanceConstraintFunction(const Eigen::Array<double,4,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      a0 << sconst[0], sconst[1], sconst[2];
      len = sconst[3];
    }

    Scalar operator() (const Eigen::Matrix<Scalar,6,1>& q, Scalar t) const
    {
      // inputs:
      // q[0] = x translation of node 1
      // q[1] = y translation of node 1
      // q[2] = z translation of node 1
      // q[3] = x translation of node 2
      // q[4] = y translation of node 2
      // q[5] = z translation of node 2

      Eigen::Matrix<Scalar,3,1> a = a0.template cast<Scalar>() + q.template segment<3>(3) - q.template segment<3>(0);
      return -len + a.norm();
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<> template<>
int
ConstraintJacobian<double,DistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,6,1>& q, Eigen::Matrix<double,6,1>& J) const;

template<> template<>
int
SacadoReverseJacobian<ConstraintJacobian<double,DistanceConstraintFunction> >
::operator() (const Eigen::Matrix<double,6,1>& q, Eigen::Matrix<double,6,6>& H) const;

#endif
