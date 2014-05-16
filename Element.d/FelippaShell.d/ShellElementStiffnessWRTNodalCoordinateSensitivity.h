#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_
#define _SHELLELEMENTSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

template<typename Scalar>
class ShellElementStiffnessWRTNodalCoordinateSensitivity : public MatrixValuedFunction<9,18,18,Scalar,4,0,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Array<Scalar,18,1> globalu; // nodal displacements
    Scalar E, nu, rho, h; // material properties

  public:
    ShellElementStiffnessWRTNodalCoordinateSensitivity(const Eigen::Array<double,4,1>& sconst, const Eigen::Array<int,0,1>& iconst)
    {
      globalu.setZero();
      E = sconst[0];
      nu = sconst[1];
      rho = sconst[2];
      h = sconst[3];
    }

    Eigen::Matrix<Scalar,18,18> operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar)
    {
      // inputs:
      // q = Nodal coordinates [x0, y0, z0, x1, y1, z1, x2, y2, z2]
      Eigen::Matrix<Scalar,3,1> globalx, globaly, globalz; 
      globalx << q[0], q[3], q[6];
      globaly << q[1], q[4], q[7];
      globalz << q[2], q[5], q[8];

      ele.setgpmat(new ShellMaterialType0<Scalar>(E, h, nu, rho));

      // elm      <input>   Finite Element Number                           not actually used
      // nu       <input>   Poisson's Ratio (for an Isotropic Element)
      // globalX  <input>   X- Nodal Coordinates
      // globalY  <input>   Y- Nodal Coordinates
      // globalZ  <input>   Z- Nodal Coordinates
      // estiff   <output>  Element Stiffness 
      // ctyp     <input>   Type of Constitutive Law (0, 1, 2, 3, or 4)     0 for linear elastic
      // flag     <input>   0: return local, 1: return global
      int ctyp = 0;
      int flag = 1;
      Eigen::Array<Scalar,18,18> estiff;
      ele.andesstf(0, estiff.data(), (Scalar*)NULL, nu, globalx.data(), globaly.data(), globalz.data(),
                   globalu.data(), ctyp, flag);

      // return value:
      // element stiffness matrix
      return estiff; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#include <Element.d/Function.d/SpaceDerivatives.h>

inline void ComputeShellElementStiffnessWRTNodalCoordinateSensitivity()
{
  // example function demonstrating how to compute the sensitivities for a matrix-valued function of a vector

  // scalar parameters
  Eigen::Array<double,4,1> dconst;
//  dconst.segment<18>(0) = Eigen::Array<double,18,1>::Random()*1e-6; // some small displacements
  dconst[0] = 200e9; // E
  dconst[1] = 0.3;   // nu
  dconst[2] = 7850;  // rho
  dconst[3] = 0.01;  // h
  // integer parameters
  Eigen::Array<int,0,1> iconst;
  // inputs
  Eigen::Matrix<double,9,1> q;
  q << 0.0, 0.0, 0.0, // x,y,z coordinates of node 0
       0.1, 0.0, 0.0, // x,y,z coordinates of node 1
       0.0, 0.1, 0.0; // x,y,z coordinates of node 2

  // function evaluation K(q)
  ShellElementStiffnessWRTNodalCoordinateSensitivity<double> foo(dconst,iconst);
  Eigen::Matrix<double,18,18> K = foo(q, 0);
  std::cerr << "K = \n" << K << std::endl;

#ifndef DEBUG_SHELL_ELEMENT_TEMPLATE
  // Jacobian evaluation: ∂K(q)/∂q
  Simo::FirstPartialSpaceDerivatives<double,ShellElementStiffnessWRTNodalCoordinateSensitivity> jac(dconst,iconst);
  Eigen::Array<Eigen::Matrix<double,18,18>,1,9> dKdq = jac(q, 0);
  for(int i=0; i<9; ++i) std::cerr << "i = " << i << ", dKdqi = \n" << dKdq[i] << std::endl;
#endif
};

#endif
#endif

