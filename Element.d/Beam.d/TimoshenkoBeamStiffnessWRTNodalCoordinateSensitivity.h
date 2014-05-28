#ifdef USE_EIGEN3
#ifndef _TIMOSHENKOBEAMSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_
#define _TIMOSHENKOBEAMSTIFFNESSWRTNODALCOORDINATESENSITIVITY_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Beam.d/BeamElementTemplate.cpp>

// class template to facilitate computation of the sensitivities of the stiffness w.r.t the nodal coordinates

template<typename Scalar>
class TimoshenkoBeamStiffnessWRTNodalCoordinateSensitivity : public MatrixValuedFunction<6,12,12,Scalar,18,0,double>
{
  public:
    BeamElementTemplate<Scalar> ele;
    Scalar A, E, Ixx, Iyy, Izz, alphaY, alphaZ, c1, nu; // material properties
    Eigen::Array<Scalar,9,1> elemframe;

  public:
    TimoshenkoBeamStiffnessWRTNodalCoordinateSensitivity(const Eigen::Array<double,18,1>& sconst, const Eigen::Array<int,0,1>& iconst)
    {
      for(int i=0; i<9; ++i) {
          elemframe[i] = sconst[i];
      }
      E = sconst[9];
      A = sconst[10];
      Ixx = sconst[11];
      Iyy = sconst[12];
      Izz = sconst[13];
      alphaY = sconst[14];
      alphaZ = sconst[15];
      c1 = sconst[16];
      nu = sconst[17];
    }

    Eigen::Matrix<Scalar,12,12> operator() (const Eigen::Matrix<Scalar,6,1>& q, Scalar)
    {
      // inputs:
      // q = Nodal Coordinates

      Eigen::Matrix<Scalar,2,1> globalx, globaly, globalz;
      globalx << q[0], q[3];
      globaly << q[1], q[4];
      globalz << q[2], q[5];

      Eigen::Matrix<Scalar,12,12> estiff;
      ele.modmstif7(estiff.data(), A, E, elemframe.data(),
                    Ixx, Iyy, Izz, alphaY, alphaZ, c1,   
                    nu, globalx.data(), globaly.data(), globalz.data(), 1);

      // return value:

      return estiff; 
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
#endif
