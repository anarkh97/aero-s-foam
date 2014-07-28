#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTGRAVITYFORCEWRTNODALCOORDINATESENSITIVITY_H_
#define _SHELLELEMENTGRAVITYFORCEWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterial.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/ShellMaterialType1.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the shell thickness

template<typename Scalar>
class ShellElementGravityForceWRTNodalCoordinateSensitivity : public VectorValuedFunction<9,18,Scalar,52,2,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Matrix<Scalar,3,1> gravityAcceleration;
    Scalar E, nu, rho, h; // material properties
    int gravflg;
    Eigen::Array<Scalar,9,1> cframe;   // composite frame
    Eigen::Array<Scalar,36,1> coefs;
    int type;

  public:
    ShellElementGravityForceWRTNodalCoordinateSensitivity(const Eigen::Array<double,52,1>& sconst, const Eigen::Array<int,2,1>& iconst)
    {
      gravityAcceleration = sconst.segment<3>(0).cast<Scalar>();
      E = sconst[3];
      nu = sconst[4];
      rho = sconst[5];
      h = sconst[6]; 
      gravflg = iconst[0];
      type = iconst[1];
      cframe.setZero();   coefs.setZero();
      if(type == 1) {
        cframe = sconst.segment<9>(7).cast<Scalar>();
        coefs = sconst.segment<36>(16).cast<Scalar>();
      } 
    }

    Eigen::Matrix<Scalar,18,1> operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar)
    {
      // inputs:
      // q = nodal coordinates
      
      Scalar x[3] = { q[0], q[3], q[6] };
      Scalar y[3] = { q[1], q[4], q[7] };
      Scalar z[3] = { q[2], q[5], q[8] };
      ele.setgpnmat(new ShellMaterialType0<Scalar>(E, h, nu, rho));
      if(type == 1) { 
        ele.setgpnmat(new ShellMaterialType1<Scalar>(coefs.data(), cframe.data(), rho, h)); 
      }
      else if(type == 2 || type == 3) { std::cerr << " ... Error: ShellElementStiffnessWRTNodalCoordinateSensitivity is not defined for this case\n"; exit(-1); }
      else if(type > 4)  { std::cerr << " ... Error: wrong material type\n"; exit(-1); }

      Scalar ElementMassMatrix[18][18];
      Eigen::Matrix<Scalar,18,1> gravityForce;
      Scalar grvfor[3];
      bool grvflg = true, masflg = false;
      Scalar totmas = 0;
      ele.andesgf(1, x, y, z, gravityForce.data(),
                  (Scalar *)ElementMassMatrix, gravityAcceleration.data(), grvfor, grvflg, totmas, masflg, gravflg);

      // return value:
      return gravityForce; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
