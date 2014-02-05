#ifdef USE_EIGEN3
#ifndef _THREENODESHELLSTRESSWRTTHICKNESSSENSITIVITY_H_
#define _THREENODESHELLSTRESSWRTTHICKNESSSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Shell.d/ShellElementSemiTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the element thickness

template<typename Scalar>
class ThreeNodeShellStressWRTThicknessSensitivity : public VectorValuedFunction<1,3,Scalar,29,1,double>
{
  public:
    ShellElementSemiTemplate<Scalar> ele;
    Eigen::Array<Scalar,3,1> globalx, globaly, globalz, ndTemp; // nodal coordinates
    Eigen::Array<Scalar,18,1> globalu; // nodal displacements
    Scalar E, nu, W, h, Ta; // material properties
    int surface; // thru-thickness location at which stresses are to be evaluated
    int isTemp;

  public:
    ThreeNodeShellStressWRTThicknessSensitivity(const Eigen::Array<double,29,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      globalx = sconst.segment<3>(0).cast<Scalar>();
      globaly = sconst.segment<3>(3).cast<Scalar>();
      globalz = sconst.segment<3>(6).cast<Scalar>();
      globalu = sconst.segment<18>(9).cast<Scalar>();
      E = sconst[27];
      nu = sconst[28];
      surface = iconst[0];
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,1,1>& q, Scalar)
    {
      // inputs:
      // q = Global Thicknesss at the Nodal Joints
      
      Eigen::Matrix<Scalar,3,1> thick;
      thick[0] = thick[1] = thick[2] = q[0];
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
      Scalar thermalStrain = 0.0;
      ele.sands8(globalx.data(), globaly.data(), globalz.data(), E, nu, thick.data(), 
                 globalu.data(), stress.data(), 0, surface, thermalStrain);

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

