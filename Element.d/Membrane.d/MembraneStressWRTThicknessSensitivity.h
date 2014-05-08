#ifdef USE_EIGEN3
#ifndef _MEMBRANESTRESSWRTTHICKNESSSENSITIVITY_H_
#define _MEMBRANESTRESSWRTTHICKNESSSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Membrane.d/MembraneElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the nodal von mises stress w.r.t the element thickness

template<typename Scalar>
class MembraneStressWRTThicknessSensitivity : public VectorValuedFunction<1,3,Scalar,29,1,double>
{
  public:
    MembraneElementTemplate<Scalar> ele;
    Eigen::Array<Scalar,3,1> globalx, globaly, globalz, ndTemp; // nodal coordinates
    Eigen::Array<Scalar,18,1> globalu; // nodal displacements
    Scalar E, nu, W, h, Ta; // material properties

  public:
    MembraneStressWRTThicknessSensitivity(const Eigen::Array<double,29,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      globalx = sconst.segment<3>(0).cast<Scalar>();
      globaly = sconst.segment<3>(3).cast<Scalar>();
      globalz = sconst.segment<3>(6).cast<Scalar>();
      globalu = sconst.segment<18>(9).cast<Scalar>();
      E = sconst[27];
      nu = sconst[28];
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,1,1>& h, Scalar)
    {
      // inputs:
      // q = Global Thicknesss at the Nodal Joints
      
      Eigen::Matrix<Scalar,3,1> thick;
      thick[0] = thick[1] = thick[2] = h[0];
      Eigen::Array<Scalar,7,3> stress;
      ele.sands19(globalx.data(), globaly.data(), globalz.data(), E, nu, thick.data(), 
                  globalu.data(), stress.data(), 0);

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


