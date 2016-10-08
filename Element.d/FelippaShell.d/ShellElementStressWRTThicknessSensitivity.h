#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTSTRESSWRTTHICKNESSSENSITIVITY_H_
#define _SHELLELEMENTSTRESSWRTTHICKNESSSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterial.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType1.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the shell thickness

template<typename Scalar>
class ShellElementStressWRTThicknessSensitivity : public VectorValuedFunction<1,3,Scalar,86,3,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Array<Scalar,3,1> globalx, globaly, globalz; // nodal coordinates
    Eigen::Array<Scalar,18,1> globalu; // displacements & rotations
    Scalar E, nu, rho, Ta, W; // material properties
    Eigen::Array<Scalar,9,1> cframe; // composite frame
    Eigen::Array<Scalar,42,1> coefs;
    Eigen::Array<Scalar,3,1> ndtemps;
    int surface; // thru-thickness location at which stresses are to be evaluated
    int type;
    int sflg;

  public:
    ShellElementStressWRTThicknessSensitivity(const Eigen::Array<double,86,1>& sconst, const Eigen::Array<int,3,1>& iconst)
    {
      globalx = sconst.segment<3>(0).cast<Scalar>();
      globaly = sconst.segment<3>(3).cast<Scalar>();
      globalz = sconst.segment<3>(6).cast<Scalar>();
      globalu = sconst.segment<18>(9).cast<Scalar>();
      E = sconst[27];
      nu = sconst[28];
      rho = sconst[29];
      surface = iconst[0];
      type = iconst[1];
      sflg = iconst[2];
      if(type == 1) {
        cframe = sconst.segment<9>(30).cast<Scalar>();
        coefs = sconst.segment<42>(39).cast<Scalar>();
      }
      Ta = sconst[81];
      W = sconst[82];
      ndtemps = sconst.segment<3>(83).cast<Scalar>();
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,1,1>& q, Scalar)
    {
      // inputs:
      // q[0] = shell thickness

      ShellMaterial<Scalar> *nmat;
      switch(type) {
        case 0 :
          nmat = new ShellMaterialType0<Scalar>(E, q[0], nu, rho, Ta, W);
          break;
        case 1 :
          nmat = new ShellMaterialType1<Scalar>(coefs.data(), cframe.data(), rho, q[0], Ta);
          break;
        default :
          std::cerr << " *** ERROR: ShellElementStressWRTThicknessSensitivity is not defined for this case.\n";
          exit(-1);
      }
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
      ele.andesvms(1, 7, nu, globalx.data(), globaly.data(), globalz.data(), globalu.data(),
                   stress.data(), type, nmat, 0, surface, sflg, ndtemps.data());
      delete nmat;

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