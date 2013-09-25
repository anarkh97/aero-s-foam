#ifndef _INERTIALFORCEFUNCTIONEXP_H_
#define _INERTIALFORCEFUNCTIONEXP_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/utilities.hpp>

// Compute the inertial moment for a rigid body with 3 rotational degrees of freedom.
// according to the classical Euler's equation for the rigid body using the extended Newmark 
// method (ALGO_1 from reference below)
// Ref: UNCONDITIONALLY STABLE ALGORITHMS FOR RIGID BODY DYNAMICS THAT EXACTLY PRESERVE ENERGY
//      AND MOMENTUM, J. C. SIMO AND K. K. WONG, INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN
//      ENGINEERING, VOL. 31, 19-52 (1991)

template<typename Scalar>
class InertialForceFunctionExp : public VectorValuedFunction<3,3,Scalar,24,0,double>
{
  public:
    Eigen::Matrix<double,3,3> J;         // inertial tensor and damping matrix
    Eigen::Matrix<double,3,1> A_n, V_n;  // angular velocity and acceleration vectors
                                         // note: according to current convention, A_n and V_n store the
                                         // "convected" angular velocities and accelerations
    Eigen::Matrix<double,3,1> Psi_n;
    double beta, gamma, alphaf, alpham, dt; // time integration scheme parameters
    double alphaDamp; // mass-proportional damping parameter

  public:
    InertialForceFunctionExp(const Eigen::Array<double,24,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      J = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+0);
      A_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+9);
      V_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+12);
      Psi_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+15);
      beta   = sconst[18];
      gamma  = sconst[19];
      alphaf = sconst[20];
      alpham = sconst[21];
      dt     = sconst[22];
      alphaDamp = sconst[23];
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar t)
    {
      // inputs:
      // q[0] = x component of total axis-angle rotation vector of node 1
      // q[1] = y component of total axis-angle rotation vector of node 1
      // q[2] = z component of total axis-angle rotation vector of node 1

      // compute the increment in the total rotation vector
      Eigen::Matrix<Scalar,3,1> Psi = q;
      Eigen::Matrix<Scalar,3,1> inc_displacement = Psi - Psi_n.template cast<Scalar>();

      // compute the updated total angular velocity and acceleration
      Eigen::Matrix<Scalar,3,1> V, A;
      V = gamma/(dt*beta)*inc_displacement
          + ((1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n).template cast<Scalar>();
      A = (1-alpham)/(dt*dt*beta*(1-alphaf))*inc_displacement
          + (-(1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n).template cast<Scalar>();
/*
      if(Psi.squaredNorm() != 0 && Psi.norm() > M_PI) {
        Eigen::Matrix<Scalar,3,3> B = complement_transf(Psi);
        Eigen::Matrix<Scalar,3,3> Bdot = complement_transf_dot(Psi,V);
        A = (B*A + Bdot*V).eval();
        V = (B*V).eval();
        Psi = complement_rot_vec(Psi);
      }
*/
      // compute the tangential transformation matrix T
      Eigen::Matrix<Scalar,3,3> T;
      tangential_transf(Psi, T);

      Eigen::Matrix<Scalar,3,3> Tdot;
      tangential_transf_dot(Psi, V, Tdot);

      // convected description of the angular momentum balance equation (see Simo & Wong eq. 29)
      // premultiplied by transformation matrix
      if(alphaDamp == 0) {
        return T.transpose()*
               (J.template cast<Scalar>()*(T*A + Tdot*V) + (T*V).cross(J.template cast<Scalar>()*T*V));
      }
      else {
        return T.transpose()*
               (J.template cast<Scalar>()*(T*A + Tdot*V + alphaDamp*V) + (T*V).cross(J.template cast<Scalar>()*T*V));
      }
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
