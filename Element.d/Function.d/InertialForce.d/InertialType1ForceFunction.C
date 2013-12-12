#ifdef USE_EIGEN3
#include <Element.d/Function.d/InertialForce.d/InertialType1ForceFunction.h>
#include <Element.d/Function.d/Rotation.d/IncrementalRotationVector.h>
#include <iostream>

namespace Simo {

template<>
Eigen::Matrix<double,3,3>
Jacobian<double,InertialType1ForceFunction>
::operator() (const Eigen::Matrix<double,3,1>& q, double t)
{
  Eigen::Matrix3d J, R_n, Rref, Fx, Vx, Px;
  Eigen::Vector3d A_n, V_n, A, V, Psi, F, P;

  J = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+0);
  A_n = Eigen::Map<Eigen::Vector3d>(const_cast<double*>(sconst.data())+9);
  V_n = Eigen::Map<Eigen::Vector3d>(const_cast<double*>(sconst.data())+12);
  R_n = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+15);
  Rref = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+24);
  const double &beta   = sconst[33];
  const double &gamma  = sconst[34];
  const double &alphaf = sconst[35];
  const double &alpham = sconst[36];
  const double &dt     = sconst[37];

  double s1 = gamma/(dt*beta);
  double s2 = (1-alpham)/(dt*dt*beta*(1-alphaf));
  if( (q.array() == 0).all() ) {
    mat_to_vec<double>(R_n.transpose()*Rref, Psi);
  }
  else {
    Eigen::Matrix3d R;
    vec_to_mat(q, R);
    mat_to_vec<double>(R_n.transpose()*R*Rref, Psi);
  }
  V = s1*Psi + (1-(1-alphaf)*gamma/beta)*V_n + dt*(1-alphaf)*(2*beta-gamma)/(2*beta)*A_n;
  A = s2*Psi - (1-alpham)/(dt*beta)*V_n + ((alpham-1)/(2*beta)+1)*A_n;
  F = Rref*(J*A + V.cross(J*V));
  P = J*V;
  Fx <<     0, -F[2],  F[1],
         F[2],     0, -F[0],
        -F[1],  F[0],     0;
  Vx <<     0, -V[2],  V[1],
         V[2],     0, -V[0],
        -V[1],  V[0],     0;
  Px <<     0, -P[2],  P[1],
         P[2],     0, -P[0],
        -P[1],  P[0],     0;

  Eigen::Array<double,18,1> sconst2 = sconst.segment<18>(15);
  Jacobian<double,IncrementalRotationVector> dPsidq(sconst2, Eigen::Array<int,0,1>::Zero());
  //std::cerr << "d[mat_to_vec(R_n^T*vec_to_mat(psi)*Rref)]/dq = \n" << dPsidq(q, 0.) << std::endl;

  if( (q.array() == 0).all() ) {
    return -0.5*Fx + Rref*(s2*J + s1*(Vx*J - Px))*dPsidq(q, t);
  }
  else {
    Eigen::Matrix3d T, C2;
    tangential_transf(q, T);
    directional_deriv2(q, F, C2);
    return C2 + T.transpose()*Rref*(s2*J + s1*(Vx*J - Px))*dPsidq(q, t);
  }
}

} // namespace Simo
#endif
