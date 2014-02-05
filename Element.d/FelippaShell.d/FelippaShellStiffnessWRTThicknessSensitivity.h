#ifdef USE_EIGEN3
#ifndef _FELIPPASHELLSTIFFNESSWRTTHICKNESSSENSITIVITY_H_
#define _FELIPPASHELLSTIFFNESSWRTTHICKNESSSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the nodal displacements

template<typename Scalar>
class FelippaShellStiffnessWRTThicknessSensitivity : public MatrixValuedFunction<1,18,18,Scalar,30,0,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Array<Scalar,3,1> globalx, globaly, globalz; // nodal coordinates
    Eigen::Array<Scalar,18,1> globalu;                  // nodal displacements
    Scalar E, nu, rho; // material properties

  public:
    FelippaShellStiffnessWRTThicknessSensitivity(const Eigen::Array<double,30,1>& sconst, const Eigen::Array<int,0,1>& iconst)
    {
      globalx = sconst.segment<3>(0).cast<Scalar>();
      globaly = sconst.segment<3>(3).cast<Scalar>();
      globalz = sconst.segment<3>(6).cast<Scalar>();
      globalu = sconst.segment<18>(9).cast<Scalar>(); 
      E = sconst[27];
      nu = sconst[28];
      rho = sconst[29];
    }

    Eigen::Matrix<Scalar,18,18> operator() (const Eigen::Matrix<Scalar,1,1>& h, Scalar)
    {
      // inputs:
      // h = thickness

      ele.setgpmat(new ShellMaterialType0<Scalar>(E, h[0], nu, rho));

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
      ele.andesstf(0, estiff.data(), (Scalar*)NULL, nu, globalx.data(), globaly.data(), globalz.data(), globalu.data(),
                   ctyp, flag);

      // return value:
      // element stiffness matrix

      return estiff; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif

