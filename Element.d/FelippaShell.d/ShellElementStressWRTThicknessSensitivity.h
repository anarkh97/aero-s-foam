#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTSTRESSWRTTHICKNESSSENSITIVITY_H_
#define _SHELLELEMENTSTRESSWRTTHICKNESSSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the shell thickness

template<typename Scalar>
class ShellElementStressWRTThicknessSensitivity : public VectorValuedFunction<1,3,Scalar,30,1,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Array<Scalar,3,1> globalx, globaly, globalz; // nodal coordinates
    Eigen::Array<Scalar,18,1> globalu; // displacements & rotations
    Scalar E, nu, rho; // material properties
    int surface; // thru-thickness location at which stresses are to be evaluated

  public:
    ShellElementStressWRTThicknessSensitivity(const Eigen::Array<double,30,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      globalx = sconst.segment<3>(0).cast<Scalar>();
      globaly = sconst.segment<3>(3).cast<Scalar>();
      globalz = sconst.segment<3>(6).cast<Scalar>();
      globalu = sconst.segment<18>(9).cast<Scalar>();
      E = sconst[27];
      nu = sconst[28];
      rho = sconst[29];
      surface = iconst[0];
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,1,1>& q, Scalar)
    {
      // inputs:
      // q[0] = shell thickness

      ele.setnmat(new ShellMaterialType0<Scalar>(E, q[0], nu, rho));

      // elm      <input>   Finite Element Number                           not actually used
      // maxstr   <input>   Maximum Number of Stresses                      7 (6 components of symmetric stress tensor + 1 von mises)
      // nu       <input>   Poisson's Ratio (for an Isotropic Element)
      // globalX  <input>   X- Nodal Coordinates
      // globalY  <input>   Y- Nodal Coordinates
      // globalZ  <input>   Z- Nodal Coordinates
      // globalU  <input>   Global Displacements at the Nodal Joints
      // stress   <output>  Stresses (Von Mises Stress) of the Element
      // ctyp     <input>   Type of Constitutive Law (0, 1, 2, 3, or 4)     0 for linear elastic
      // flag     <input>   0: stress, 1: strain, >= 2: internal variables
      // surface  <input>   1: upper, 2: median, 3: lower
      Eigen::Array<Scalar,7,3> stress;
      ele.andesvms(0, 7, nu, globalx.data(), globaly.data(), globalz.data(), globalu.data(),
                   stress.data(), 0, 0, surface);

      // return value:
      // von mises stresses at nodes
      Eigen::Matrix<Scalar,3,1> v;
      v << stress(6,0), stress(6,1), stress(6,2);

      return v; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/*
#include <Element.d/Function.d/SpaceDerivatives.h>
inline void ComputeShellElementStressWRTThicknessSensitivity()
{
  // example function demonstrating how to compute the sensitivities

  // scalar parameters
  Eigen::Array<double,30,1> dconst;
  dconst[0] = 0; dconst[1] = 0.1; dconst[2] = 0.0; // x coordinates
  dconst[3] = 0; dconst[4] = 0;   dconst[5] = 0.1; // y coordinates
  dconst[6] = 0; dconst[7] = 0;   dconst[8] = 0;   // z coordinates
  dconst.segment<18>(9) = Eigen::Array<double,18,1>::Random()*1e-6; // some small displacements
  dconst[27] = 200e9; // E
  dconst[28] = 0.3;   // nu
  dconst[29] = 7850;  // rho
  // integer parameters
  Eigen::Array<int,1,1> iconst;
  iconst[0] = 1; // upper surface
  // inputs
  Eigen::Matrix<double,1,1> q;
  q[0] = 1e-2; // value of thickness at which jacobian is to be evaluated

  // function evaluation
  ShellElementStressWRTThicknessSensitivity<double> foo(dconst,iconst);
  Eigen::Matrix<double,3,1> S = foo(q, 0);
  std::cerr << "S = " << S.transpose() << std::endl;

  // Jacobian evaluation
  Simo::Jacobian<double,ShellElementStressWRTThicknessSensitivity> dSdh(dconst,iconst);
  Eigen::Matrix<double,3,1> J = dSdh(q, 0);
  std::cerr << "J = " << J.transpose() << std::endl;
};
*/

#endif
#endif
