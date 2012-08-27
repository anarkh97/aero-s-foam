#ifndef _LINELINEDISTANCECONSTRAINTFUNCTION_H_
#define _LINELINEDISTANCECONSTRAINTFUNCTION_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/ConstraintFunction.h>

template<typename Scalar>
class LineLineDistanceConstraintFunction : public RheonomicConstraintFunction<6,Scalar,17,1,double>
{
    // constrains the distance (d) between two lines (defined by points p0, p1, and q0,q1 respectively)
    // according to d - (A*sin(omega*t+phi) + (B-C*t)*d0) = 0, <= 0 or >= 0
    // see: http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment
    // see: http://mathworld.wolfram.com/Line-LineDistance.html 
    Eigen::Matrix<double,3,1> p0;
    Eigen::Matrix<double,3,1> p1;
    Eigen::Matrix<double,3,1> q0;
    Eigen::Matrix<double,3,1> q1;
    double A, omega, phase, B, C, d0;
    bool negate;

  public:
    LineLineDistanceConstraintFunction(const Eigen::Array<double,17,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      p0 = sconst.segment<3>(0);
      p1 = sconst.segment<3>(3);
      q0 = sconst.segment<3>(6);
      q1 = sconst.segment<3>(9);
      A = sconst[12];
      omega = sconst[13];
      phase = sconst[14];
      B = sconst[15];
      C = sconst[16];
      negate = bool(iconst[0]);
      /*
      Eigen::Matrix<double,3,1> w0 = p0 - q0; 
      Eigen::Matrix<double,3,1> u = p1 - p0; 
      Eigen::Matrix<double,3,1> v = q1 - q0; 
      double a = u.dot(u);
      double b = u.dot(v);
      double c = v.dot(v);
      double D = u.dot(w0);
      double e = v.dot(w0);
      //d0 = (w0 + ((b*e-c*D)*u - (a*e-b*D)*v)/(a*c-b*b)).norm();      
      */
	using std::abs;
	//Check to see if lines are not nearly parallel
	if (1.0 - abs((p1-p0).normalized().dot((q1-q0).normalized())) > 1.0e-6){
    		d0 = abs((q0-p0).dot(((p1-p0).cross(q1-q0)).normalized()));
	}else{
		//formula for point to line using p0 as initial point
		 d0 = (p0-q0).cross(p0-q1).norm()/(q1-q0).norm();
	}
    }

    Scalar operator() (const Eigen::Matrix<Scalar,6,1>& q, Scalar t) const
    {
      // q[0] = x translation of point 0
      // q[1] = y translation of point 0
      // q[2] = z translation of point 0
      // q[3] = x translation of point 1
      // q[4] = y translation of point 1
      // q[5] = z translation of point 1
      Eigen::Matrix<Scalar,3,1> p0 = LineLineDistanceConstraintFunction::p0.template cast<Scalar>() + q.template segment<3>(0);
      Eigen::Matrix<Scalar,3,1> p1 = LineLineDistanceConstraintFunction::p1.template cast<Scalar>() + q.template segment<3>(3);
      Eigen::Matrix<Scalar,3,1> q0 = LineLineDistanceConstraintFunction::q0.template cast<Scalar>();
      Eigen::Matrix<Scalar,3,1> q1 = LineLineDistanceConstraintFunction::q1.template cast<Scalar>();
     /*
      Eigen::Matrix<double,3,1> w0 = p0 - q0; 
      Eigen::Matrix<double,3,1> u = p1 - p0; 
      Eigen::Matrix<double,3,1> v = q1 - q0; 
      double a = u.dot(u);
      double b = u.dot(v);
      double c = v.dot(v);
      double D = u.dot(w0);
      double e = v.dot(w0);
      Scalar d = (w0 + ((b*e-c*D)*u - (a*e-b*D)*v)/(a*c-b*b)).norm();
     */
       using std::abs;
	Scalar d; 
       	if (1.0 - abs((p1-p0).normalized().dot((q1-q0).normalized())) > 1.0e-6){
      		 d = abs((q0-p0).dot(((p1-p0).cross(q1-q0)).normalized()));
	}else{
		 d = (p0-q0).cross(p0-q1).norm()/(q1-q0).norm();
	}
      using std::sin;
      Scalar f = d -(A*sin(omega*t+phase) + (B-C*t)*d0);
      if(negate) return -f; else return f;
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/*template<> template<>
int
ConstraintJacobian<double,LineLineDistanceConstraintFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, Eigen::Matrix<double,3,1>& J) const;

template<> template<>
int
SacadoReverseJacobian<ConstraintJacobian<double,LineLineDistanceConstraintFunction> >
::operator() (const Eigen::Matrix<double,3,1>& q, Eigen::Matrix<double,3,3>& H) const;
*/
#endif
