#ifndef ROM_EIGALERKINPROJECTIONSOLVER_H
#define ROM_EIGALERKINPROJECTIONSOLVER_H

#ifdef USE_EIGEN3
#include <Solvers.d/Solver.h>
#include <Math.d/EiSparseMatrix.h>

#include "VecBasis.h"
#include "BasisOps.h"

#include <cstddef>
#include <stdexcept>

#ifdef ROM_TIMING
#include <Eigen/Dense>
#endif

namespace Rom {

template <typename Scalar>
class GenEiSparseGalerkinProjectionSolver : public GenPodProjectionSolver<Scalar>, public GenEiSparseMatrix<Scalar> {
public:
  GenEiSparseGalerkinProjectionSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa);

  using GenEiSparseMatrix<Scalar>::neqs;
/*
  // Full-order matrix assembly
  virtual void zeroAll();
  virtual void add(GenFullSquareMatrix<Scalar> &, int *);
  virtual void addDiscreteMass(int, Scalar);
*/
  // Solution
  virtual void factor();
  virtual void reSolve(GenVector<Scalar> &rhs);
  double projectAndComputeNorm(const GenVector<Scalar> &rhs); // next reSolve must use same rhs

  // Reduced basis parameters
  int basisSize() const { return V.cols(); }
  const GenVecBasis<Scalar> &projectionBasis() const { 
    std::cerr << "ERROR: GenEiSparseGalerkinProjectionSolver does not implement projectionBasis\n";
    exit(-1);
  }
  void projectionBasisIs(const GenVecBasis<Scalar> &); // Passed objects must be kept alive by owner
  
  // Data collection
  const GenVector<Scalar> &lastReducedSolution() const { 
    std::cerr << "ERROR: GenEiSparseGalerkinProjectionSolver does not implement lastReducedSolution\n";
    exit(-1);
  }
  const GenVecBasis<Scalar> &lastReducedMatrixAction() const { 
    std::cerr << "ERROR: GenEiSparseGalerkinProjectionSolver does not implement lastReducedMatrixAction\n";
    exit(-1); 
  }

protected:
  Eigen::Matrix<Scalar,Eigen::Dynamic,1> &getReducedSolution() { return reducedSolution; }

private:
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> V;
  Eigen::Matrix<Scalar,Eigen::Dynamic,1> reducedSolution;
  Eigen::LLT<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> > eiSolver;
  bool rhsIsprojected;
  
  void validateRhs(const GenVector<Scalar> &);

  // Disallow copy and assignment
  GenEiSparseGalerkinProjectionSolver(const GenEiSparseGalerkinProjectionSolver<Scalar> &);
  GenEiSparseGalerkinProjectionSolver<Scalar> &operator=(const GenEiSparseGalerkinProjectionSolver<Scalar> &);
};

template <typename Scalar>
GenEiSparseGalerkinProjectionSolver<Scalar>::GenEiSparseGalerkinProjectionSolver(Connectivity *cn,
                                                   DofSetArray *dsa,
                                                   ConstrainedDSA *c_dsa):
  GenEiSparseMatrix<Scalar>(cn, dsa, c_dsa),
  rhsIsprojected(false)
{
}

/*
template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::zeroAll()
{
  GenEiSparseMatrix<Scalar>::zeroAll();
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::add(GenFullSquareMatrix<Scalar> &elMat, int *dofs)
{
  GenEiSparseMatrix<Scalar>::add(elMat, dofs);
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::addDiscreteMass(int dof, Scalar dmass)
{
  GenEiSparseMatrix<Scalar>::addDiscreteMass(dof, dmass);
}
*/

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::projectionBasisIs(const GenVecBasis<Scalar> &reducedBasis)
{
  if (reducedBasis.vectorSize() != neqs()) {
    throw std::domain_error("Vectors of the reduced basis have the wrong size");
  }

  reducedSolution = Eigen::Matrix<Scalar,Eigen::Dynamic,1>::Zero(reducedBasis.vectorCount());
  rhsIsprojected = false;
  V.resize(reducedBasis.vectorSize(), reducedBasis.vectorCount());
  for(int j=0; j<reducedBasis.vectorCount(); ++j)
    for(int i=0; i<reducedBasis.vectorSize(); ++i)
      V(i,j) = reducedBasis[j][i];
}

template <typename Scalar>
inline
void
GenEiSparseGalerkinProjectionSolver<Scalar>::validateRhs(const GenVector<Scalar> &rhs)
{
  if (rhs.size() != neqs()) {
    throw std::domain_error("Rhs has the wrong size");
  }
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::factor()
{
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> eiReducedMatrix(V.cols(),V.cols());
  eiReducedMatrix.template triangularView<Eigen::Lower>()
    = V.transpose()*(this->M.template selfadjointView<Eigen::Upper>()*V);

  eiSolver.compute(eiReducedMatrix);
}
   
template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::reSolve(GenVector<Scalar> &rhs)
{
  validateRhs(rhs);

  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > x(rhs.data(), V.rows());
  if (rhsIsprojected) {
    x = V*eiSolver.solve(reducedSolution);
  }
  else {
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> b = Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> >(rhs.data(), V.rows());
    x = V*eiSolver.solve(V.transpose()*b);
  }
}

template <typename Scalar>
double
GenEiSparseGalerkinProjectionSolver<Scalar>::projectAndComputeNorm(const GenVector<Scalar> &rhs)
{
  validateRhs(rhs);

  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > b(rhs.data(), V.rows());
  reducedSolution = V.transpose()*b;
  return reducedSolution.norm();
}

typedef GenEiSparseGalerkinProjectionSolver<double> EiSparseGalerkinProjectionSolver;

} /* end namespace Rom */
#endif

#endif /* ROM_EIGALERKINPROJECTIONSOLVER_H */
