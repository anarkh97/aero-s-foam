#ifndef _POINTVARIPLANEDISTANCECONSTRAINT_H_
#define _POINTVARIPLANEDISTANCECONSTRAINT_H_

#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/ConstraintFunction.h>

template<typename Scalar>
class PointVariPlaneDistanceConstraintFunction : public RheonomicConstraintFunction<12,Scalar,17,1,double>
{
    // constrains the distance (d) between a point and a variable plane (defined by three points x1, x2 and x3) according to
    // d - (A*sin(omega*t+phi) + (B-C*t)*d0) = 0, <= 0 or >= 0
    // see: http://mathworld.wolfram.com/Point-PlaneDistance.html

    Eigen::Matrix<double,3,1> x0;
    Eigen::Matrix<double,3,1> x1;
    Eigen::Matrix<double,3,1> x2;
    Eigen::Matrix<double,3,1> x3;
    double A, omega, phase, B, C, d0;
    bool negate;

  public:
    PointVariPlaneDistanceConstraintFunction(const Eigen::Array<double,17,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      x0 = sconst.segment(0,3);
      x1 = sconst.segment(3,3);
      x2 = sconst.segment(6,3);
      x3 = sconst.segment(9,3);
      A = sconst[12];
      omega = sconst[13];
      phase = sconst[14];
      B = sconst[15];
      C = sconst[16];
      negate = bool(iconst[0]);

      d0 = (x2-x1).cross(x3-x1).normalized().dot(x0-x1); // initial distance
    }

    Scalar operator() (const Eigen::Matrix<Scalar,12,1>& q, Scalar t) const
    {
      // q(0) = x translation of point 0
      // q(1) = y translation of point 0
      // q(2) = z translation of point 0
      // q(3) = x translation of point 1
      // q(4) = y translation of point 1
      // q(5) = z translation of point 1
      // q(6) = x translation of point 2
      // q(7) = y translation of point 2
      // q(8) = z translation of point 2
      // q(9) = x translation of point 3
      // q(10) = y translation of point 3
      // q(11) = z translation of point 3
      Eigen::Matrix<Scalar,3,1> x0 = PointVariPlaneDistanceConstraintFunction::x0.template cast<Scalar>() + q.segment(0,3);
      Eigen::Matrix<Scalar,3,1> x1 = PointVariPlaneDistanceConstraintFunction::x1.template cast<Scalar>() + q.segment(3,3);
      Eigen::Matrix<Scalar,3,1> x2 = PointVariPlaneDistanceConstraintFunction::x1.template cast<Scalar>() + q.segment(6,3);
      Eigen::Matrix<Scalar,3,1> x3 = PointVariPlaneDistanceConstraintFunction::x1.template cast<Scalar>() + q.segment(9,3);

      Scalar d = (x2-x1).cross(x3-x1).normalized().dot(x0-x1);
      // note: the sign of d is positive if x0 is on the same side of the plane as the normal
      //       and negative if it is on the opposite side
      using std::sin;
      Scalar f = d -(A*sin(omega*t+phase) + (B-C*t)*d0);
      if(negate) return -f; else return f;
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointVariPlaneDistanceConstraint : public ConstraintFunctionElement<PointVariPlaneDistanceConstraintFunction>
{
  public:
    PointVariPlaneDistanceConstraint(int* _nn); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst);
};
#endif
#endif
