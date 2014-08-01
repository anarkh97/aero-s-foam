#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTSTRESSWRTNODALCOORDINATESENSITIVITY_H_
#define _SHELLELEMENTSTRESSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterial.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType1.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

template<typename Scalar>
class ShellElementStressWRTNodalCoordinateSensitivity : public VectorValuedFunction<9,3,Scalar,69,2,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Array<Scalar,18,1> globalu; // element displacements
    Scalar E, nu, rho, h, Ta, W; // material properties
    int surface; // thru-thickness location at which stresses are to be evaluated
    Eigen::Array<Scalar,9,1> cframe;   // composite frame
    Eigen::Array<Scalar,36,1> coefs;
    int type;

  public:
    ShellElementStressWRTNodalCoordinateSensitivity(const Eigen::Array<double,69,1>& sconst, const Eigen::Array<int,2,1>& iconst)
    {
      globalu = sconst.segment<18>(0).cast<Scalar>();
      E = sconst[18];
      nu = sconst[19];
      rho = sconst[20];
      h = sconst[21];
      surface = iconst[0];
      type = iconst[1];
      if(type == 1) {
        cframe = sconst.segment<9>(22).cast<Scalar>();
        coefs = sconst.segment<36>(31).cast<Scalar>(); 
      }
      Ta = sconst[67];
      W = sconst[68];
    }
    

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar)
    {
      // inputs:
      // q = Global Displacements at the Nodal Joints

      ele.setgpnmat(new ShellMaterialType0<Scalar>(E, h, nu, rho, Ta, W));
      if(type == 1) { 
        ele.setgpnmat(new ShellMaterialType1<Scalar>(coefs.data(), cframe.data(), rho, h, Ta)); 
      }
      else if(type == 2 || type == 3) { std::cerr << " ... Error: ShellElementStiffnessWRTNodalCoordinateSensitivity is not defined for this case\n"; exit(-1); }
      else if(type > 4)  { std::cerr << " ... Error: wrong material type\n"; exit(-1); }
      Eigen::Array<Scalar,3,1> globalx; 
      Eigen::Array<Scalar,3,1> globaly; 
      Eigen::Array<Scalar,3,1> globalz; 
      globalx << q[0], q[3], q[6];
      globaly << q[1], q[4], q[7];
      globalz << q[2], q[5], q[8];

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
      std::cerr << stress(6,0) << " " << stress(6,1) << " " << stress(6,2) << std::endl;

      return v; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
