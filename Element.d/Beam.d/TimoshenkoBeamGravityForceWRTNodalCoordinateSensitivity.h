#ifdef USE_EIGEN3
#ifndef _TIMOSHENKOBEAMGRAVITYFORCEWRTNODALCOORDINATESENSITIVITY_H_
#define _TIMOSHENKOBEAMGRAVITYFORCEWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Beam.d/BeamElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of gravity force w.r.t the nodal displacements

template<typename Scalar>
class TimoshenkoBeamGravityForceWRTNodalCoordinateSensitivity : public VectorValuedFunction<6,12,Scalar,13,1,double>
{
  public:
    BeamElementTemplate<Scalar> ele;
    Scalar massPerNode; // material properties
    Eigen::Array<Scalar,9,1> elemframe;
    Eigen::Array<Scalar,3,1> gravityAcceleration;
    int gravflg;

  public:
    TimoshenkoBeamGravityForceWRTNodalCoordinateSensitivity(const Eigen::Array<double,13,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      for(int i=0; i<3; ++i) gravityAcceleration[i] = sconst[i];
      for(int i=0; i<9; ++i) elemframe[i] = sconst[3+i];
      massPerNode = sconst[12];

      gravflg = iconst[0];
    }

    Eigen::Matrix<Scalar,12,1> operator() (const Eigen::Matrix<Scalar,6,1>& q, Scalar)
    {
      // inputs:
      // q = Global Displacements at the Nodal Joints

      Eigen::Matrix<Scalar,2,1> globalx, globaly, globalz;
      globalx << q[0], q[3];
      globaly << q[1], q[4];
      globalz << q[2], q[5];

      Eigen::Matrix<Scalar,12,1> gravityForce;
      ele.gForce(massPerNode, globalx.data(), globaly.data(), globalz.data(),
                 gravityAcceleration.data(), elemframe.data(), gravityForce.data(), gravflg);

      return gravityForce; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
