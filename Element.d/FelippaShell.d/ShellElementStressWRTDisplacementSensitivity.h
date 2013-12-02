#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTSTRESSWRTDISPLACEMENTSENSITIVITY_H_
#define _SHELLELEMENTSTRESSWRTDISPLACEMENTSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the nodal displacements

template<typename Scalar>
class ShellElementStressWRTDisplacementSensitivity : public VectorValuedFunction<18,3,Scalar,13,1,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Array<Scalar,3,1> globalx, globaly, globalz; // nodal coordinates
    Scalar E, nu, rho, h; // material properties
    int surface; // thru-thickness location at which stresses are to be evaluated

  public:
    ShellElementStressWRTDisplacementSensitivity(const Eigen::Array<double,13,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      globalx = sconst.segment<3>(0).cast<Scalar>();
      globaly = sconst.segment<3>(3).cast<Scalar>();
      globalz = sconst.segment<3>(6).cast<Scalar>();
      E = sconst[9];
      nu = sconst[10];
      rho = sconst[11];
      h = sconst[12];
      surface = iconst[0];
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,18,1>& q, Scalar)
    {
      // inputs:
      // q = Global Displacements at the Nodal Joints

      ele.setnmat(new ShellMaterialType0<Scalar>(E, h, nu, rho));
      Eigen::Matrix<Scalar,18,1> globalu = q;

      // elm      <input>   Finite Element Number                           not actually used
      // maxstr   <input>   Maximum Number of Stresses                      7 (6 components of symmetric stress tensor + 1 von mises)
      // nu       <input>   Poisson's Ratio (for an Isotropic Element)
      // globalX  <input>   X- Nodal Coordinates
      // globalY  <input>   Y- Nodal Coordinates
      // globalZ  <input>   Z- Nodal Coordinates
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

#endif
#endif
