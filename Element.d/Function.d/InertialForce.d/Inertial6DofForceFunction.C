#ifdef USE_EIGEN3
#include <Element.d/Function.d/InertialForce.d/Inertial6DofForceFunction.h>
#include <Element.d/Function.d/Rotation.d/IncrementalRotationVector.h>
#include <iostream>

namespace Simo {

template<>
Eigen::Matrix<double,6,6>
Jacobian<double,Inertial6DofForceFunction>
::operator() (const Eigen::Matrix<double,6,1>& q, double t)
{
  Eigen::Vector3d c, u_n, F, G, P, Q;
  Eigen::Matrix3d J, J0, C, E, R_n, Rref, Fx, Gx, Omegax, uddotx, Px, Qx;
  Eigen::Matrix<double,6,1> j, a_n, v_n, a, v, inc_displacement;

  const double &m = sconst[0];
  j = Eigen::Map<Eigen::Matrix<double,6,1> >(const_cast<double*>(sconst.data())+1);
  c = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+7);
  a_n = Eigen::Map<Eigen::Matrix<double,6,1> >(const_cast<double*>(sconst.data())+10);
  v_n = Eigen::Map<Eigen::Matrix<double,6,1> >(const_cast<double*>(sconst.data())+16);
  u_n = Eigen::Map<Eigen::Matrix<double,3,1> >(const_cast<double*>(sconst.data())+22);
  R_n = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+25);
  Rref = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+34);
  const double &beta   = sconst[43];
  const double &gamma  = sconst[44];
  const double &alphaf = sconst[45];
  const double &alpham = sconst[46];
  const double &dt     = sconst[47];

  J << j[0],  j[5],  j[4],
       j[5],  j[1],  j[3],
       j[4],  j[3],  j[2];
  C <<  0.0, -c[2],  c[1],
       c[2],   0.0, -c[0],
      -c[1],  c[0],   0.0;
  J0 = J - m*C*C;

  double s1 = gamma/(dt*beta);
  double s2 = (1-alpham)/(dt*dt*beta*(1-alphaf));

  inc_displacement.head<3>() = q.head<3>() - u_n;
  if( (q.tail<3>().array() == 0).all() ) {
    mat_to_vec<double>(R_n.transpose()*Rref, inc_displacement.tail<3>());
  }
  else {
    Eigen::Matrix3d dR;
    vec_to_mat<double>(q.tail<3>(), dR);
    mat_to_vec<double>(R_n.transpose()*dR*Rref, inc_displacement.tail<3>());
  }

  v = s1*inc_displacement + (1-(1-alphaf)*gamma/beta)*v_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*a_n;
  a = s2*inc_displacement - (1-alpham)/(dt*beta)*v_n + ((alpham-1)/(2*beta)+1)*a_n;

  Eigen::VectorBlock<Eigen::Matrix<double,6,1>,3> udot  = v.head<3>(), Omega = v.tail<3>(), 
                                                  uddot = a.head<3>(), Alpha = a.tail<3>(); 
  E = Rref*C*Rref.transpose();
  F = m*E*uddot + Rref*(J0*Alpha + Omega.cross(J0*Omega));
  G = Rref*(C*Alpha + Omega.cross(C*Omega));
  P = J0*Omega;
  Q = C*Omega;
  Fx <<     0, -F[2],  F[1],
         F[2],     0, -F[0],
        -F[1],  F[0],     0;
  Gx <<     0, -G[2],  G[1],
         G[2],     0, -G[0],
        -G[1],  G[0],     0;
  Omegax <<         0, -Omega[2],  Omega[1],
             Omega[2],         0, -Omega[0],
            -Omega[1],  Omega[0],         0;
  uddotx <<         0, -uddot[2],  uddot[1],
             uddot[2],         0, -uddot[0],
            -uddot[1],  uddot[0],         0;
  Px <<     0, -P[2],  P[1],
         P[2],     0, -P[0],
        -P[1],  P[0],     0;
  Qx <<     0, -Q[2],  Q[1],
         Q[2],     0, -Q[0],
        -Q[1],  Q[0],     0;

  Eigen::Array<double,18,1> sconst2 = sconst.segment<18>(25);
  Jacobian<double,IncrementalRotationVector> dPsidq(sconst2, Eigen::Array<int,0,1>::Zero());
  Eigen::Matrix3d D = dPsidq(q.tail<3>(), t);

  Eigen::Matrix<double,6,6> K;
  if( (q.tail<3>().array() == 0).all() ) {
    K.topLeftCorner<3,3>()     = s2*m*Eigen::Matrix3d::Identity();
    K.topRightCorner<3,3>()    = -m*(-Gx + Rref*(s2*C + s1*(Omegax*C - Qx))*D);
    K.bottomLeftCorner<3,3>()  = s2*m*E;
    K.bottomRightCorner<3,3>() = -0.5*Fx + m*E*uddotx + Rref*(s2*J0 + s1*(Omegax*J0 - Px))*D;
  }
  else {
    std::cerr << "Error: Jacobian<double,Inertial6DofForceFunction>::operator() is not implemented for Δψ != 0\n";
    //Eigen::Matrix3d T, C2;
    //tangential_transf(q, T);
    //directional_deriv2(q, F, C2);
    //return C2 + T.transpose()*Rref*(s2*J + s1*(Vx*J - Px))*dPsidq(q, t);
  }

  return K;
}

} // namespace Simo
#endif
