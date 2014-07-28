#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTMASSWRTNODALCOORDINATESENSITIVITY_H_
#define _SHELLELEMENTMASSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
//#include <Element.d/FelippaShell.d/ShellMaterial.hpp>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
//#include <Element.d/FelippaShell.d/ShellMaterialType1.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

template<typename Scalar>
class ShellElementMassWRTNodalCoordinateSensitivity : public VectorValuedFunction<9,1,Scalar,49,1,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Scalar E, nu, rho, h; // material properties
    Eigen::Array<Scalar,9,1> cframe;   // composite frame
    Eigen::Array<Scalar,36,1> coefs;
    int type;

  public:
    ShellElementMassWRTNodalCoordinateSensitivity(const Eigen::Array<double,49,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      E = sconst[0];
      nu = sconst[1];
      rho = sconst[2];
      h = sconst[3];
      type = iconst[0];
      cframe.setZero();   coefs.setZero();
      if(type == 1) {
        cframe = sconst.segment<9>(4).cast<Scalar>();
        coefs = sconst.segment<36>(13).cast<Scalar>();
      }
    }

    Eigen::Matrix<Scalar,1,1> operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar)
    {
      // inputs:
      // q = Global Displacements at the Nodal Joints

      ele.setgpmat(new ShellMaterialType0<Scalar>(E, h, nu, rho));
//      if(type == 1) { 
//        ele.setgpnmat(new ShellMaterialType1<Scalar>(coefs.data(), cframe.data(), rho, h)); 
//      }
//      else if(type == 2 || type == 3) { std::cerr << " ... Error: ShellElementStiffnessWRTNodalCoordinateSensitivity is not defined for this case\n"; exit(-1); }
//      else if(type > 4)  { std::cerr << " ... Error: wrong material type\n"; exit(-1); }
      Eigen::Array<Scalar,3,1> globalx; 
      Eigen::Array<Scalar,3,1> globaly; 
      Eigen::Array<Scalar,3,1> globalz; 
      globalx << q[0], q[3], q[6];
      globaly << q[1], q[4], q[7];
      globalz << q[2], q[5], q[8];

      Scalar *grvfor = NULL;
      bool grvflg = false, massflg = true;
      Eigen::Array<Scalar,3,1> gravityAcceleration;
      Eigen::Matrix<Scalar,18,18> ElementMassMatrix;

      Scalar totmas = 0.0;
      ele.andesms(1, globalx.data(), globaly.data(), globalz.data(), ElementMassMatrix.data(),
                   gravityAcceleration.data(), grvfor, grvflg, totmas, massflg);
 
      Eigen::Matrix<Scalar,1,1> mass;
      mass[0] = totmas;

      return mass; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
