#include "VecBasis.h"

namespace Rom{

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::project(GenDistrVector<double> &x, GenDistrVector<double> &_result) {

#ifdef USE_EIGEN3
  Eigen::Matrix<double, Eigen::Dynamic, 1> GenCoordinates;
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > f(x.data(), x.size());
  Eigen::SparseVector<double> sparsef(f.size());
//  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > basis(vectors_[0].data(), vectors_[0].size(), numVectors());

  for (int i = 0; i < f.size(); i++) {
   if(f(i) != 0.) {
    sparsef.insert(i) = f(i);
    }
  }
//  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Unaligned, Eigen::OuterStride<> > basisT(vectors_[0].data(), vectors_[0].size(), vectorCount(), Eigen::OuterStride<>(size()));

  GenCoordinates = basis.transpose()*sparsef;

  if(structCom)
    structCom->globalSum(GenCoordinates.size(), GenCoordinates.data());
  
  result = basis*GenCoordinates;
#else
  cerr << "USE_EIGEN3 is not defined here in GenVecBasis::project\n";
  exit(-1);
#endif
  return _result;
}

}
