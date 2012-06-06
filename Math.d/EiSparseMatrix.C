#include <iostream>
#include <Math.d/Vector.h>
#include <Driver.d/Communicator.h>

template<class Scalar>
GenEiSparseMatrix<Scalar>::GenEiSparseMatrix(Connectivity *cn, DofSetArray *dsa, int *rCN, bool _selfadjoint)
: SparseData(dsa,rCN,cn,int(!_selfadjoint)),
  selfadjoint(_selfadjoint),
  nnz(xunonz[numUncon]),
  unonz(new Scalar[nnz]),
  M(numUncon, numUncon, nnz, xunonz, rowu, unonz)
#ifdef EIGEN_SUPERLU_SUPPORT
  ,lu(0)
#endif
{
  for(int i=0; i<numUncon+1; ++i) xunonz[i]--;
  for(int i=0; i<nnz; ++i) rowu[i]--;
  zeroAll();
}

template<class Scalar>
GenEiSparseMatrix<Scalar>::GenEiSparseMatrix(Connectivity *cn, DofSetArray *dsa, DofSetArray *c_dsa, bool _selfadjoint)
: SparseData(dsa,c_dsa,cn,int(!_selfadjoint)),
  selfadjoint(_selfadjoint),
  nnz(xunonz[numUncon]),
  unonz(new Scalar[nnz]),
  M(numUncon, numUncon, nnz, xunonz, rowu, unonz)
#ifdef EIGEN_SUPERLU_SUPPORT
  ,lu(0)
#endif
{
  for(int i=0; i<numUncon+1; ++i) xunonz[i]--;
  for(int i=0; i<nnz; ++i) rowu[i]--;
  zeroAll();
}

template<class Scalar>
GenEiSparseMatrix<Scalar>::GenEiSparseMatrix(Connectivity *cn, EqNumberer *eqNums)
: SparseData(cn,eqNums,1.0E-6,0),
  selfadjoint(true),
  nnz(xunonz[numUncon]),
  unonz(new Scalar[nnz]),
  M(numUncon, numUncon, nnz, xunonz, rowu, unonz)
#ifdef EIGEN_SUPERLU_SUPPORT
  ,lu(0)
#endif
{
  for(int i=0; i<numUncon+1; ++i) xunonz[i]--;
  for(int i=0; i<nnz; ++i) rowu[i]--;
  zeroAll();
}

template<class Scalar>
GenEiSparseMatrix<Scalar>::~GenEiSparseMatrix() 
{
  if(unonz) { delete [] unonz; unonz=0; }
#ifdef EIGEN_SUPERLU_SUPPORT
  if(lu) delete lu;
#endif
}

template<class Scalar> 
void
GenEiSparseMatrix<Scalar>::zeroAll()
{
  for(int i=0; i<nnz; ++i) unonz[i] = 0;
}

template<class Scalar> 
void
GenEiSparseMatrix<Scalar>::print()
{
  std::cerr << M << std::endl;
}

template<class Scalar>
void
GenEiSparseMatrix<Scalar>::add(int idof, int jdof, Scalar s)
{
  // TODO non-selfadjoint case
  if((idof < 0) || (jdof < 0)) return;
  int dofsi = (jdof > idof) ? idof : jdof;
  int dofsj = (jdof > idof) ? jdof : idof;

  int k,l;
  if((k = unconstrNum[dofsi]) == -1 || (l = unconstrNum[dofsj]) == -1) return;
  int mstart = xunonz[l];
  int mstop  = xunonz[l+1];
  for(int m = xunonz[l]; m < xunonz[l+1]; ++m) {
    if(rowu[m] == k) {
      unonz[m] += s;
      break;
    }
  }
}

template<class Scalar>
void
GenEiSparseMatrix<Scalar>::addDiscreteMass(int dof, Scalar dmass)
{
  int cdof = unconstrNum[dof];
  if(cdof < 0) return;
  int diagLocator = xunonz[cdof+1]-1; // This should be the diagonal
  unonz[diagLocator] += dmass;
}

template<class Scalar> 
void
GenEiSparseMatrix<Scalar>::add(FullSquareMatrix &kel, int *dofs)
{
  int k,l;
  for(int i = 0; i < kel.dim(); ++i) {
    if(dofs[i] < 0 || (k = unconstrNum[dofs[i]]) < 0) continue; // Skip undefined/constrained dofs
    for(int j = 0; j < kel.dim(); ++j) {
      if(selfadjoint && dofs[i] > dofs[j]) continue; // Work with upper symmetric half
      if(dofs[j] < 0 || (l = unconstrNum[dofs[j]]) < 0) continue;  // Skip undefined/constrained dofs
      for(int m = xunonz[l]; m < xunonz[l+1]; ++m) {
        if(rowu[m] == k) {
          unonz[m] += kel[i][j];
          break;
        }
      }
    }
  }
}

template<class Scalar> 
double
GenEiSparseMatrix<Scalar>::getMemoryUsed()
{
  return sizeof(Scalar)*nnz/(1024.0*1024.0);
}

template<class Scalar> 
void
GenEiSparseMatrix<Scalar>::mult(const Scalar *_rhs, Scalar *_result)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(const_cast<Scalar*>(_rhs),numUncon,1), result(_result,numUncon,1);
  if(selfadjoint)
    result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
    result = M*rhs;
}

template<class Scalar>
void
GenEiSparseMatrix<Scalar>::mult(const GenVector<Scalar> &_rhs, Scalar *_result) 
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1), result(_result,numUncon,1);
  if(selfadjoint)
    result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
    result = M*rhs;
}

template<class Scalar>
void
GenEiSparseMatrix<Scalar>::transposeMult(const GenVector<Scalar> &_rhs, GenVector<Scalar> &_result)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1), result(_result.data(),numUncon,1);
  if(selfadjoint)
    result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
    result = M.adjoint()*rhs;
}

template<class Scalar>
void
GenEiSparseMatrix<Scalar>::multAdd(const Scalar *_rhs, Scalar *_result)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(const_cast<Scalar*>(_rhs),numUncon,1), result(_result,numUncon,1);
  if(selfadjoint)
    result += M.template selfadjointView<Eigen::Upper>()*rhs;
  else
    result += M*rhs;
}

template<class Scalar> 
void
GenEiSparseMatrix<Scalar>::mult(const GenVector<Scalar> &_rhs, GenVector<Scalar> &_result)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1), result(_result.data(),numUncon,1);
  if(selfadjoint)
    result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
    result = M*rhs;
}

template<class Scalar> 
Scalar
GenEiSparseMatrix<Scalar>::diag(int dof) const
{
  std::cerr << "GenEiSparseMatrix<Scalar>::diag is not implemented\n";
  return Scalar();
}

template<class Scalar> 
Scalar &
GenEiSparseMatrix<Scalar>::diag( int dof )
{
  static Scalar defaultValue;
  std::cerr << "GenEiSparseMatrix<Scalar>::diag is not implemented\n";
  return defaultValue;
}

template<class Scalar> 
long
GenEiSparseMatrix<Scalar>::size()
{
  return (numUncon) ? nnz : 0;
}

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
#endif

template<class Scalar> 
void
GenEiSparseMatrix<Scalar>::unify(FSCommunicator *communicator)
{
#ifdef DISTRIBUTED
  communicator->globalSum(nnz, unonz);
#endif
}

template<class Scalar>
void
GenEiSparseMatrix<Scalar>::factor()
{
  if(selfadjoint)
    llt.compute(M);
  else {
#ifdef EIGEN_SUPERLU_SUPPORT
    if(lu) delete lu;
    lu = new Eigen::SuperLU<Eigen::SparseMatrix<Scalar> >(M);
    if(lu->info() != Eigen::Success) std::cerr << "factor with SuperLU failed\n";
#else
    std::cerr << "You need SuperLU to factor this matrix\n";
    exit(-1);
#endif
  }
}

template<class Scalar>
void 
GenEiSparseMatrix<Scalar>::solve(Scalar *_rhs, Scalar *_solution)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1), solution(_solution,numUncon,1);
  if(selfadjoint)
    solution = llt.solve(rhs);
  else {
#if EIGEN_SUPERLU_SUPPORT
    //Eigen::Matrix<Scalar, Eigen::Dynamic, 1> rhs_copy = rhs;
    //rhs_copy = lu->solve(rhs_copy);
    //solution = rhs_copy;
    solution = lu->solve(rhs);
    if(lu->info() != Eigen::Success) std::cerr << "solve with SuperLU failed\n";
#else
    std::cerr << "You need SuperLU to solve this system\n";
    exit(-1);
#endif
  }
}

template<class Scalar>
void
GenEiSparseMatrix<Scalar>::solve(GenVector<Scalar> &_rhs, GenVector<Scalar> &_solution)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1), solution(_solution.data(),numUncon,1);
  if(selfadjoint)
    solution = llt.solve(rhs);
  else {
#if EIGEN_SUPERLU_SUPPORT
    //Eigen::Matrix<Scalar, Eigen::Dynamic, 1> rhs_copy = rhs;
    //rhs_copy = lu->solve(rhs_copy);
    //solution = rhs_copy;
    solution = lu->solve(rhs);
    if(lu->info() != Eigen::Success) std::cerr << "solve with SuperLU failed\n";
#else
    std::cerr << "You need SuperLU to solve this system\n";
    exit(-1);
#endif
  }
}

template<class Scalar>
void
GenEiSparseMatrix<Scalar>::reSolve(Scalar *_rhs)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  if(selfadjoint)
    rhs = llt.solve(rhs);
  else {
#if EIGEN_SUPERLU_SUPPORT
    //Eigen::Matrix<Scalar, Eigen::Dynamic, 1> rhs_copy = rhs;
    //rhs_copy = lu->solve(rhs_copy);
    //rhs = rhs_copy;
    rhs = lu->solve(rhs);
    if(lu->info() != Eigen::Success) std::cerr << "solve with SuperLU failed\n";
#else
    std::cerr << "You need SuperLU to solve this system\n";
    exit(-1);
#endif
  }
}

template<class Scalar>
void
GenEiSparseMatrix<Scalar>::reSolve(GenVector<Scalar> &_rhs)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1);
  if(selfadjoint)
    rhs = llt.solve(rhs);
  else {
#if EIGEN_SUPERLU_SUPPORT
    //Eigen::Matrix<Scalar, Eigen::Dynamic, 1> rhs_copy = rhs;
    //rhs_copy = lu->solve(rhs_copy);
    //rhs = rhs_copy;
    rhs = lu->solve(rhs);
    if(lu->info() != Eigen::Success) std::cerr << "solve with SuperLU failed\n";
#else
    std::cerr << "You need SuperLU to solve this system\n";
    exit(-1);
#endif
  }
}
