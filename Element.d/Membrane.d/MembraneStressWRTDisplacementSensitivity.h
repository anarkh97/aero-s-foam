#ifdef USE_EIGEN3
#ifndef _MEMBRANESTRESSWRTDISPLACEMENTSENSITIVITY_H_
#define _MEMBRANESTRESSWRTDISPLACEMENTSENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Membrane.d/MembraneElementTemplate.cpp>

template<typename Scalar>
class MembraneStressWRTDisplacementSensitivity : public VectorValuedFunction<18,3,Scalar,12,1,double>
{
  public:
    MembraneElementTemplate<Scalar> ele;
    Eigen::Array<Scalar,3,1> globalx, globaly, globalz, ndTemp; // nodal coordinates
    Scalar E, nu, W, h, Ta; // material properties
    int surface; // thru-thickness location at which stresses are to be evaluated
    int isTemp;

  public:
    MembraneStressWRTDisplacementSensitivity(const Eigen::Array<double,12,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      globalx = sconst.segment<3>(0).cast<Scalar>();
      globaly = sconst.segment<3>(3).cast<Scalar>();
      globalz = sconst.segment<3>(6).cast<Scalar>();
      E = sconst[9];
      nu = sconst[10];
      h = sconst[11];
      surface = iconst[0];
    }

    Eigen::Matrix<Scalar,3,1> operator() (const Eigen::Matrix<Scalar,18,1>& q, Scalar)
    {
      // inputs:
      // q = Global Displacements at the Nodal Joints
      
      Eigen::Matrix<Scalar,18,1> globalu = q;
      Eigen::Matrix<Scalar,3,1> thick;
      thick[0] = thick[1] = thick[2] = h;
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

