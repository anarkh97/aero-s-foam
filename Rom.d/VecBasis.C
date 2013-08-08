#include "VecBasis.h"
#include "Utils.d/DistHelper.h"
#include "Math.d/Vector.h"

namespace Rom {

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::project(GenDistrVector<double> &x, GenDistrVector<double> &_result) {
//Nothing to do here, Fext & Fint are already projected
// code left for reference
/*
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
*/
  return _result;
}

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::projectUp(GenDistrVector<double> &x, GenDistrVector<double> &_result) {
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > GenCoordinates(x.data(), x.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());

  //full coordinates distributed over MPI processes
  if(compressedKey.size() > 0) {
    Eigen::VectorXd resultBuffer(compressedKey.size());
    resultBuffer = compressedBasis*GenCoordinates;
    for(int i = 0; i < compressedKey.size(); i++)
      result(compressedKey[i]) = resultBuffer(i);
  }
  else {
    result = basis*GenCoordinates;
  }
#endif
  return _result;
}

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::projectUp(std::vector<double> &x, GenDistrVector<double> &_result) {
#ifdef USE_EIGEN3
  //this instantiation is for the post processor, need to fix it to use GenDistrVector
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > GenCoordinates(x.data(), x.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());

  //full coordinates distributed over MPI processes
  result = basis*GenCoordinates;
#endif
  return _result;
}

template <>
GenDistrVector<double> &
GenVecBasis<double, GenDistrVector>::projectDown(GenDistrVector<double> &x, GenDistrVector<double> &_result) {
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > FullCoordinates(x.data(), x.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > result(_result.data(), _result.size());
// fix this portion so we can use the compressed Basis, its way faster but 
// gives the wrong result, just use the sparse vec mult for now 
  if(compressedKey.size() > 0) {
    Eigen::VectorXd coordBuffer(compressedKey.size());
    for(int i = 0; i < compressedKey.size(); i++)
       coordBuffer(i) = FullCoordinates(compressedKey[i]); 
    result = compressedBasis.transpose()*coordBuffer;
/*    Eigen::SparseVector<double> sparsef(FullCoordinates.rows());
    for(int i = 0; i < FullCoordinates.rows(); i++){
      if(FullCoordinates(i) != 0){
        sparsef.insert(i) = FullCoordinates(i);}}
    result = basis.transpose()*sparsef;*/
  }
  else {
    result = basis.transpose()*FullCoordinates;
  }
  //each process gets a copy of reduced coordinates
  if(structCom)
    structCom->globalSum(result.size(), result.data());
#endif 
  return _result;
}

template<>
void
GenVecBasis<double, GenDistrVector>::makeSparseBasis(const std::vector<std::vector<int> > & nodeVec, DofSetArray **dsa)
{
#ifdef USE_EIGEN3
  int dof1, numdofs;

  int dof0 = 0;
  for(int n=0; n<nodeVec.size(); n++) {
    for(int i = 0; i < nodeVec[n].size(); i++) {
      dof1 = dsa[n]->firstdof(nodeVec[n][i]);
      numdofs = dsa[n]->weight(nodeVec[n][i]);
      for(int j = 0; j < numdofs; j++) {
        compressedKey.push_back(dof0+dof1+j);
      }
    }
    dof0 += dsa[n]->size();
  }

  new (&compressedBasis) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(compressedKey.size(),vectorCount()); // O Col major, 1 RowMajor

  for(int i = 0; i < compressedKey.size(); i++) {
    for(int j = 0; j < vectorCount(); j++){
      compressedBasis(i,j) = basis(compressedKey[i],j);
    }
  }
#endif
}

template<>
void
GenVecBasis<double, GenVector>::makeSparseBasis(const std::vector<int> & nodeVec, DofSetArray *dsa)
{
#ifdef USE_EIGEN3
  int dof1, numdofs;

  for(int i = 0; i < nodeVec.size(); i++) {
    dof1 = dsa->firstdof(nodeVec[i]);
    numdofs = dsa->weight(nodeVec[i]);
    for(int j = 0; j < numdofs; j++) {
      compressedKey.push_back(dof1+j);
    }
  }

  new (&compressedBasis) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(compressedKey.size(),vectorCount()); // O Col major, 1 RowMajor

  for(int i = 0; i < compressedKey.size(); i++) {
    for(int j = 0; j < vectorCount(); j++) {
      compressedBasis(i,j) = basis(compressedKey[i],j);
    }
  }
#endif
}

}
