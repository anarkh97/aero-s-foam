#ifndef _POINTPLANEDISTANCECONSTRAINT_H_
#define _POINTPLANEDISTANCECONSTRAINT_H_

#include <Element.d/MpcElement.d/ConstraintFunction.h>

template<typename Scalar>
class PointPlaneDistanceConstraintFunction : public RheonomicConstraintFunction<3,Scalar,17,1,double>
{
    // constrains the distance (d) between a point (x0) and a fixed plane (defined by three points x1, x2 and x3) according to
    // d - (A*sin(omega*t+phi) + (BC*t)*d0) = 0, <= 0 or >= 0
    // see: http://mathworld.wolfram.com/Point-PlaneDistance.html

    Eigen::Matrix<double,3,1> x0;
    Eigen::Matrix<double,3,1> x1;
    Eigen::Matrix<double,3,1> x2;
    Eigen::Matrix<double,3,1> x3;
    double A, omega, phase, B, C, d0;
    bool negate;

    Eigen::Matrix<double,3,1> nhat; // unit normal

  public:
    PointPlaneDistanceConstraintFunction(const Eigen::Array<double,17,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      x0 = sconst.segment<3>(0);
      x1 = sconst.segment<3>(3);
      x2 = sconst.segment<3>(6);
      x3 = sconst.segment<3>(9);
      A = sconst[12];
      omega = sconst[13];
      phase = sconst[14];
      B = sconst[15];
      C = sconst[16];
      negate = bool(iconst[0]);

      nhat = (x2-x1).cross(x3-x1).normalized();
      d0 = nhat.dot(x0-x1);
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar t) const
    {
      // q(0) = x translation of point 0
      // q(1) = y translation of point 0
      // q(2) = z translation of point 0
      Eigen::Matrix<Scalar,3,1> x0 = PointPlaneDistanceConstraintFunction::x0.template cast<Scalar>() + q;
      Eigen::Matrix<Scalar,3,1> x1 = PointPlaneDistanceConstraintFunction::x1.template cast<Scalar>();
      Eigen::Matrix<Scalar,3,1> nhat = PointPlaneDistanceConstraintFunction::nhat.template cast<Scalar>();

      Scalar d = nhat.dot(x0-x1);
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

class PointPlaneDistanceConstraint : public ConstraintFunctionElement<PointPlaneDistanceConstraintFunction>
{
    double x1[3], x2[3], x3[3]; // coordinates of the 3 points defining the plane

  public:
    PointPlaneDistanceConstraint(int* _nn); 
    void setFrame(EFrame *);

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst);
};

#endif
