#include <Driver.d/EFrameData.h>
#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

void
NFrameData::transformVector(double *data, bool hasRot)
{
  // transform 3x1 or 6x1 node vector from basic to DOF_FRM coordinates
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);

  if(hasRot) {
    Eigen::Map<Eigen::Matrix<double,6,1> > v(data);
    v.head<3>() = (T*v.head<3>()).eval();
    v.tail<3>() = (T*v.tail<3>()).eval();
  }
  else {
    Eigen::Map<Eigen::Vector3d> v(data);
    v = (T*v).eval();
  }
#endif
}

void
NFrameData::transformVectorInv(double *data, bool hasRot)
{
  // transform 3x1 or 6x1 node vector from DOF_FRM to basic coordinates
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);

  if(hasRot) {
    Eigen::Map<Eigen::Matrix<double,6,1> > v(data);
    v.head<3>() = (T.transpose()*v.head<3>()).eval();
    v.tail<3>() = (T.transpose()*v.tail<3>()).eval();
  }
  else {
    Eigen::Map<Eigen::Vector3d> v(data);
    v = (T.transpose()*v).eval();
  }
#endif
}

