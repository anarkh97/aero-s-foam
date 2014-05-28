#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTMASSWRTNODALCOORDINATESENSITIVITY_H_
#define _SHELLELEMENTMASSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/FelippaShell.d/ShellMaterialType0.cpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>
#include <Element.d/FelippaShell.d/ShellElementTemplate.cpp>

template<typename Scalar>
class ShellElementMassWRTNodalCoordinateSensitivity : public VectorValuedFunction<9,1,Scalar,4,0,double>
{
  public:
    ShellElementTemplate<Scalar,EffMembraneTriangle,AndesBendingTriangle> ele;
    Scalar E, nu, rho, h; // material properties

  public:
    ShellElementMassWRTNodalCoordinateSensitivity(const Eigen::Array<double,4,1>& sconst, const Eigen::Array<int,0,1>& iconst)
    {
      E = sconst[0];
      nu = sconst[1];
      rho = sconst[2];
      h = sconst[3];
    }

    Eigen::Matrix<Scalar,1,1> operator() (const Eigen::Matrix<Scalar,9,1>& q, Scalar)
    {
      // inputs:
      // q = Global Displacements at the Nodal Joints

      ele.setgpmat(new ShellMaterialType0<Scalar>(E, h, nu, rho));
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
