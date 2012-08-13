#ifndef _POINTVARILINEDISTANCECONSTRAINT_H_
#define _POINTVARILINEDISTANCECONSTRAINT_H_

#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/ConstraintFunction.h>

template<typename Scalar>
class PointVariLineDistanceConstraintFunction : public RheonomicConstraintFunction<9,Scalar,14,1,double>
{
    // constrains the distance (d) between a point (x0) and a variable line (defined by 2 points x1 and x2) according to
    // d - (A*sin(omega*t+phi) + (B-C*t)*d0) = 0, <= 0 or >= 0
    // see: http://mathworld.wolfram.com/Point-VariLineDistance3-Dimensional.html
    
    Eigen::Matrix<double,3,1> x0;
    Eigen::Matrix<double,3,1> x1;
    Eigen::Matrix<double,3,1> x2;
    double A, omega, phase, B, C, d0;
    bool negate;

  public:
    PointVariLineDistanceConstraintFunction(const Eigen::Array<double,14,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      x0 = sconst.segment(0,3);
      x1 = sconst.segment(3,3);
      x2 = sconst.segment(6,3);
      A = sconst[9];
      omega = sconst[10];
      phase = sconst[11];
      B = sconst[12];
      C = sconst[13];
      negate = bool(iconst[0]);

      d0 = (x0-x1).cross(x0-x2).norm()/(x2-x1).norm();
    }

    Scalar operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar t) const
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
      Eigen::Matrix<Scalar,3,1> x0 = PointVariLineDistanceConstraintFunction::x0.template cast<Scalar>() + q.segment(0,3);
      Eigen::Matrix<Scalar,3,1> x1 = PointVariLineDistanceConstraintFunction::x1.template cast<Scalar>() + q.segment(3,3);
      Eigen::Matrix<Scalar,3,1> x2 = PointVariLineDistanceConstraintFunction::x2.template cast<Scalar>() + q.segment(6,3);

      Scalar d = (x0-x1).cross(x0-x2).norm()/(x2-x1).norm();
      using std::sin;
      Scalar f = d -(A*sin(omega*t+phase) + (B-C*t)*d0);
      if(negate) return -f; else return f;
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointVariLineDistanceConstraint : public ConstraintFunctionElement<PointVariLineDistanceConstraintFunction>
{
  public:
    PointVariLineDistanceConstraint(int* _nn); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,14,1>& sconst, Eigen::Array<int,1,1>& iconst);
};
#endif
#endif
