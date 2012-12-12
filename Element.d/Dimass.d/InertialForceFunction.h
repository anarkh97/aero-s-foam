#ifndef _INERTIALFORCEFUNCTION_H_
#define _INERTIALFORCEFUNCTION_H_

#include <Element.d/Function.d/VectorValuedFunction.h>
#include <Element.d/Function.d/utilities.hpp>

// Compute the inertial force for a rigid body with 3 rotational degrees of freedom.

template<typename Scalar>
class InertialForceFunction : public VectorValuedFunction<3,3,Scalar,39,1,double>
{
  public:
    Eigen::Matrix<double,3,3> elM, elC;  // inertial tensor and damping matrix
    Eigen::Matrix<double,3,1> a_n, v_n;  // velocity and acceleration vectors
                                         // note: according to current convention, a_n and v_n store the "spatial" angular velocities and accelerations
    Eigen::Matrix<double,3,3> R_n;
    Eigen::Matrix<double,3,3> Rref;
    double beta, gamma, alphaf, alpham, dt; // time integration scheme parameters
    double alphaDamp; // mass-proportional damping parameter
    int algo;

  public:
    InertialForceFunction(const Eigen::Array<double,39,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      elM = Eigen::Map<Eigen::Matrix<double,3,3>,Eigen::RowMajor >(const_cast<double*>(sconst.data())+0);
      a_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+9);
      v_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+12);
      R_n = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+15);
      Rref = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+24);
      beta   = sconst[33];
      gamma  = sconst[34];
      alphaf = sconst[35];
      alpham = sconst[36];
      dt     = sconst[37];
      alphaDamp = sconst[38];
      algo = iconst[0];

      elC = alphaDamp*elM;
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,3,1>& q, Scalar t) const
    {
      // inputs:
      // q[0] = x component of incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1
      // q[1] = y component of incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1
      // q[2] = z component of incremental axis-angle rotation vector (w.r.t. reference configuration) of node 1

      Eigen::Matrix<Scalar,3,1> inc_displacement;
      Eigen::Matrix<Scalar,3,3> R, T;

      get_inc_displacement(q, inc_displacement, R);
      tangential_transf(q, T);

      if(algo == 1) {
        Eigen::Matrix<Scalar,3,1> V, A;
        V = gamma/(dt*beta)*inc_displacement
            + (1-(1-alphaf)*gamma/beta)*(R_n.transpose()*v_n).template cast<Scalar>()
            + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*(R_n.transpose()*a_n).template cast<Scalar>();

        A = (1-alpham)/(dt*dt*beta*(1-alphaf))*inc_displacement
            - (1-alpham)/(dt*beta)*(R_n.transpose()*v_n).template cast<Scalar>()
            + ((alpham-1)/(2*beta)+1)*(R_n.transpose()*a_n).template cast<Scalar>();

        // convected description of the angular momentum balance equation (see Simo & Wong eq. 29)
        // note: when SO3param == 2, T = I. However, we still multiply by T so that the Jacobian will correctly evaluated
        //       when this function is automatically or numerically differentiated
        return T.transpose()*Rref.template cast<Scalar>()*
               (elM.template cast<Scalar>()*A + elC.template cast<Scalar>()*V + V.cross(elM.template cast<Scalar>()*V));
      }
      else if(algo == 2) {
        Eigen::Matrix<Scalar,3,1> v, a;
        v = gamma/(dt*beta)*inc_displacement 
            + (1-(1-alphaf)*gamma/beta)*v_n.template cast<Scalar>() 
            + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*a_n.template cast<Scalar>();

        a = (1-alpham)/(dt*dt*beta*(1-alphaf))*inc_displacement 
            - (1-alpham)/(dt*beta)*v_n.template cast<Scalar>()
            + ((alpham-1)/(2*beta)+1)*a_n.template cast<Scalar>();

        Eigen::Matrix<Scalar,3,3> M = R*elM.template cast<Scalar>()*R.transpose();

        // spatial description of the angular momentum balance equation (see Simo & Wong eq. 29)
        // note: when SO3param == 2, T = I. However, we still multiply by T so that the Jacobian will correctly evaluated
        //       when this function is automatically or numerically differentiated
        return T*(M*a + elC.template cast<Scalar>()*v + v.cross(M*v));
      }
    }

  public:
    void get_inc_displacement(const Eigen::Matrix<Scalar,3,1>& q, Eigen::Matrix<Scalar,3,1>& inc_displacement,
                              Eigen::Matrix<Scalar,3,3>& R) const
    {
      Eigen::Matrix<Scalar,3,3> dR, incR;

      vec_to_mat<Scalar>(q, dR); // dR is either the spatial (left) or incremental rotation matrix w.r.t Rref
      R = dR*Rref.template cast<Scalar>(); // R is the current total rotation i.e. the R_{n+1-alphaf} at Newton iteration k

      if(algo == 1)
        incR = R_n.template cast<Scalar>().transpose()*R; // here incR is the material (right) incremental rotation matrix w.r.t R_n
                                                          // i.e. R = R_n*incR --> incR = R_n.transpose()*R
      else if(algo == 2) 
        incR = R*R_n.template cast<Scalar>().transpose(); // here incR is the spatial (left) incremental rotation matrix w.r.t R_n
                                                          // i.e. R = incR*R_n --> incR = R*R_n.transpose()

      mat_to_vec<Scalar>(incR, inc_displacement); 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
