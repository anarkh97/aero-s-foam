#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_
#define _SHELLELEMENTSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterial.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType1.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

template<typename Scalar>
class ShellElementStiffnessWRTNodalCoordinateSensitivity : public MatrixValuedFunction<9,18,18,Scalar,49,1,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Array<Scalar,18,1> globalu; // nodal displacements
    Scalar E, nu, rho, h;              // material properties
    Eigen::Array<Scalar,9,1> cframe;   // composite frame
    Eigen::Array<Scalar,36,1> coefs;
    int type;

  public:
    ShellElementStiffnessWRTNodalCoordinateSensitivity(const Eigen::Array<double,49,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      globalu.setZero();
      E = sconst[0];
      nu = sconst[1];
      rho = sconst[2];
      h = sconst[3];
      type = iconst[0];
      cframe.setZero();    coefs.setZero();
      if(type == 1) {
        cframe = sconst.segment<9>(4).cast<Scalar>();
        coefs = sconst.segment<36>(13).cast<Scalar>(); 
      } 
    }

    Eigen::Matrix<Scalar,18,18> operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar)
    {
      // inputs:
      // q = Nodal coordinates [x0, y0, z0, x1, y1, z1, x2, y2, z2]
      Eigen::Matrix<Scalar,3,1> globalx, globaly, globalz; 
      globalx << q[0], q[3], q[6];
      globaly << q[1], q[4], q[7];
      globalz << q[2], q[5], q[8];

      ele.setgpnmat(new ShellMaterialType0<Scalar>(E, h, nu, rho)); 
      if(type == 1) { 
        ele.setgpnmat(new ShellMaterialType1<Scalar>(coefs.data(), cframe.data(), rho, h)); 
      }
      else if(type == 2 || type == 3) { std::cerr << " ... Error: ShellElementStiffnessWRTNodalCoordinateSensitivity is not defined for this case\n"; exit(-1); }
      else if(type > 4)  { std::cerr << " ... Error: wrong material type\n"; exit(-1); }
//      if(type == 2) ele.setgpnmat(new ShellMaterialTypes2And3<Scalar>(nlays, lData, false, cframe)); 
//      if(type == 3) ele.setgpnmat(new ShellMaterialTypes2And3<Scalar>(nlays, lData, true, cframe)); 

      // elm      <input>   Finite Element Number                           not actually used
      // nu       <input>   Poisson's Ratio (for an Isotropic Element)
      // globalX  <input>   X- Nodal Coordinates
      // globalY  <input>   Y- Nodal Coordinates
      // globalZ  <input>   Z- Nodal Coordinates
      // estiff   <output>  Element Stiffness 
      // type     <input>   Type of Constitutive Law (0, 1, 2, 3, or 4)     0 for linear elastic
      // flag     <input>   0: return local, 1: return global
      int flag = 1;
      Eigen::Array<Scalar,18,18> estiff;
      ele.andesstf(0, estiff.data(), (Scalar*)NULL, nu, globalx.data(), globaly.data(), globalz.data(),
                   globalu.data(), type, flag);

      // return value:
      // element stiffness matrix
      return estiff; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
