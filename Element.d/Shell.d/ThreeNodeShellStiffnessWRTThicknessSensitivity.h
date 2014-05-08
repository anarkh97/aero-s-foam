#ifdef USE_EIGEN3
#ifndef _THREENODESHELLSTIFFNESSWRTTHICKNESSSENSITIVITY_H_
#define _THREENODESHELLSTIFFNESSWRTTHICKNESSSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Shell.d/ShellElementSemiTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the nodal displacements

template<typename Scalar>
class ThreeNodeShellStiffnessWRTThicknessSensitivity : public MatrixValuedFunction<1,18,18,Scalar,11,1,double>
{
  public:
    ShellElementSemiTemplate<Scalar> ele;
    Eigen::Array<Scalar,3,1> globalx, globaly, globalz; // nodal coordinates
    Scalar E, nu;  // material properties
    int flg;

  public:
    ThreeNodeShellStiffnessWRTThicknessSensitivity(const Eigen::Array<double,11,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      globalx = sconst.segment<3>(0).cast<Scalar>();
      globaly = sconst.segment<3>(3).cast<Scalar>();
      globalz = sconst.segment<3>(6).cast<Scalar>();
      E = sconst[9];
      nu = sconst[10];
      flg = iconst[0];
    }

    Eigen::Matrix<Scalar,18,18> operator() (const Eigen::Matrix<Scalar,1,1>& h, Scalar)
    {
      // inputs:
      // h = thickness

      Scalar thickness[3];
      thickness[0] = thickness[1] = thickness[2] = h[0];
      Eigen::Array<Scalar,18,18> estiff;
      ele.tria3d(flg, globalx.data(), globaly.data(), globalz.data(), E, nu, thickness, estiff.data()); 

      return estiff; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif


