#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTGRAVITYFORCEWRTTHICKNESSSENSITIVITY_H_
#define _SHELLELEMENTGRAVITYFORCEWRTTHICKNESSSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the shell thickness

template<typename Scalar>
class ShellElementGravityForceWRTThicknessSensitivity : public VectorValuedFunction<1,18,Scalar,15,1,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Eigen::Matrix<Scalar,3,1> gravityAcceleration, globalx, globaly, globalz;
    Scalar E, nu, rho; // material properties
    int gravflg;

  public:
    ShellElementGravityForceWRTThicknessSensitivity(const Eigen::Array<double,15,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      gravityAcceleration = sconst.segment<3>(0).cast<Scalar>();
      globalx = sconst.segment<3>(3).cast<Scalar>();
      globaly = sconst.segment<3>(6).cast<Scalar>();
      globalz = sconst.segment<3>(9).cast<Scalar>();
      E = sconst[12];
      nu = sconst[13];
      rho = sconst[14];
      gravflg = iconst[0];
    }

    Eigen::Matrix<Scalar,18,1> operator() (const Eigen::Matrix<Scalar,1,1>& q, Scalar)
    {
      // inputs:
      // q[0] = shell thickness

      ele.setgpmat(new ShellMaterialType0<Scalar>(E, q[0], nu, rho));

      Eigen::Matrix<Scalar,18,18> ElementMassMatrix;
      Eigen::Matrix<Scalar,18,1> gravityForce;
      Eigen::Matrix<Scalar,3,1> grvfor;
      bool grvflg = true, masflg = false;
      Scalar totmas = 0;
      ele.andesgf(1, globalx.data(), globaly.data(), globalz.data(), gravityForce.data(),
                  ElementMassMatrix.data(), gravityAcceleration.data(), grvfor.data(), grvflg, totmas, masflg, gravflg);

      // return value:
      return gravityForce; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
