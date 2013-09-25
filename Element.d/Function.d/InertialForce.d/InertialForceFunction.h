#ifndef _INERTIALFORCEFUNCTION_H_
#define _INERTIALFORCEFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/utilities.hpp>

// Compute the inertial moment for a rigid body with 3 rotational degrees of freedom.
// according to the classical Euler's equation for the rigid body using the extended Newmark 
// method (ALGO_1 from reference below)
// Ref: UNCONDITIONALLY STABLE ALGORITHMS FOR RIGID BODY DYNAMICS THAT EXACTLY PRESERVE ENERGY
//      AND MOMENTUM, J. C. SIMO AND K. K. WONG, INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN
//      ENGINEERING, VOL. 31, 19-52 (1991)

template<typename Scalar>
class InertialForceFunction : public VectorValuedFunction<3,3,Scalar,39,0,double>
{
  public:
    Eigen::Matrix<double,3,3> J;         // inertial tensor and damping matrix
    Eigen::Matrix<double,3,1> A_n, V_n;  // angular velocity and acceleration vectors
                                         // note: according to current convention, A_n and V_n store the
                                         // "convected" angular velocities and accelerations
    Eigen::Matrix<double,3,3> R_n;
    Eigen::Matrix<double,3,3> Rref;
    double beta, gamma, alphaf, alpham, dt; // time integration scheme parameters
    double alphaDamp; // mass-proportional damping parameter

  public:
    InertialForceFunction(const Eigen::Array<double,39,1>& sconst, const Eigen::Array<int,0,1>&)
    {
      J = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+0);
      A_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+9);
      V_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+12);
      R_n = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+15);
      Rref = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+24);
      beta   = sconst[33];
      gamma  = sconst[34];
      alphaf = sconst[35];
      alpham = sconst[36];
      dt     = sconst[37];
      alphaDamp = sconst[38];
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar t)
    {
      // inputs:
      // q[0] = x component of spatial incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1
      // q[1] = y component of spatial incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1
      // q[2] = z component of spatial incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1

      // compute incRref, the spatial (left) incremental rotation matrix w.r.t Rref
      Eigen::Matrix<Scalar,3,3> incRref;
      vec_to_mat<Scalar>(q, incRref); 

      // compute incR_n, the material (right) incremental rotation matrix w.r.t R_n, and it's corresponding rotation vector
      Eigen::Matrix<Scalar,3,3> incR_n;
      incR_n = R_n.template cast<Scalar>().transpose()*incRref*Rref.template cast<Scalar>();
      Eigen::Matrix<Scalar,3,1> inc_displacement;
      mat_to_vec<Scalar>(incR_n, inc_displacement);

      // compute the tangential transformation matrix T
      Eigen::Matrix<Scalar,3,3> T;
      tangential_transf(q, T);

      // compute the updated convected angular velocity acceleration
      Eigen::Matrix<Scalar,3,1> V, A;
      V = gamma/(dt*beta)*inc_displacement
          + ((1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n).template cast<Scalar>();
      A = (1-alpham)/(dt*dt*beta*(1-alphaf))*inc_displacement
          + (-(1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n).template cast<Scalar>();

      // convected description of the angular momentum balance equation (see Simo & Wong eq. 29)
      // premultiplied by transformation to fixed reference frame
      // note: Even though T(0) = I we still multiply by T.transpose() so that the Jacobian will correctly evaluated
      //       when this function is automatically or numerically differentiated
      if(alphaDamp == 0) {
        return T.transpose()*Rref.template cast<Scalar>()*
               (J.template cast<Scalar>()*A + V.cross(J.template cast<Scalar>()*V));
      }
      else {
        return T.transpose()*Rref.template cast<Scalar>()*
               (J.template cast<Scalar>()*(A + alphaDamp*V) + V.cross(J.template cast<Scalar>()*V));
      }
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
