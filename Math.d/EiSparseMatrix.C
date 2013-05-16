#include <iostream>
#include <algorithm>
#include <Math.d/Vector.h>
#include <Driver.d/Communicator.h>

template<typename Scalar, typename SolverClass>
GenEiSparseMatrix<Scalar,SolverClass>::GenEiSparseMatrix(Connectivity *cn, DofSetArray *dsa, int *rCN, bool _selfadjoint)
: SparseData(dsa,rCN,cn,int(!_selfadjoint)),
  selfadjoint(_selfadjoint),
  nnz(xunonz[numUncon]),
  unonz(new Scalar[nnz]),
  M(numUncon, numUncon, nnz, xunonz, rowu, unonz)
{
  for(int k=0; k < numUncon; k++)
    std::sort(rowu + xunonz[k]-1, rowu + xunonz[k+1]-1);
 
  for(int i=0; i<numUncon+1; ++i) xunonz[i]--;
  for(int i=0; i<nnz; ++i) rowu[i]--;
  zeroAll();
}

template<typename Scalar, typename SolverClass>
GenEiSparseMatrix<Scalar,SolverClass>::GenEiSparseMatrix(Connectivity *cn, DofSetArray *dsa, DofSetArray *c_dsa, bool _selfadjoint)
: SparseData(dsa,c_dsa,cn,int(!_selfadjoint)),
  selfadjoint(_selfadjoint),
  nnz(xunonz[numUncon]),
  unonz(new Scalar[nnz]),
  M(numUncon, numUncon, nnz, xunonz, rowu, unonz)
{
  for(int k=0; k < numUncon; k++)
    std::sort(rowu + xunonz[k]-1, rowu + xunonz[k+1]-1);

  for(int i=0; i<numUncon+1; ++i) xunonz[i]--;
  for(int i=0; i<nnz; ++i) rowu[i]--;
  zeroAll();
}

template<typename Scalar, typename SolverClass>
GenEiSparseMatrix<Scalar,SolverClass>::GenEiSparseMatrix(Connectivity *cn, EqNumberer *eqNums)
: SparseData(cn,eqNums,1.0E-6,0),
  selfadjoint(true),
  nnz(xunonz[numUncon]),
  unonz(new Scalar[nnz]),
  M(numUncon, numUncon, nnz, xunonz, rowu, unonz)
{
  for(int k=0; k < numUncon; k++)
    std::sort(rowu + xunonz[k]-1, rowu + xunonz[k+1]-1);

  for(int i=0; i<numUncon+1; ++i) xunonz[i]--;
  for(int i=0; i<nnz; ++i) rowu[i]--;
  zeroAll();
}

template<typename Scalar, typename SolverClass>
GenEiSparseMatrix<Scalar,SolverClass>::~GenEiSparseMatrix() 
{
  if(unonz) { delete [] unonz; unonz=0; }
}

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::zeroAll()
{
  for(int i=0; i<nnz; ++i) unonz[i] = 0;
}

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::print()
{
  std::cerr << M << std::endl;
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::add(int idof, int jdof, Scalar s)
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

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::addDiscreteMass(int dof, Scalar dmass)
{
  GenFullSquareMatrix<Scalar> m(1);
  m[0][0] = dmass;
  add(m, &dof);
}

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::add(FullSquareMatrix &kel, int *dofs)
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

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::addImaginary(FullSquareMatrix &kel, int *dofs)
{
  int k,l;
  for(int i = 0; i < kel.dim(); ++i) {
    if(dofs[i] < 0 || (k = unconstrNum[dofs[i]]) < 0) continue; // Skip undefined/constrained dofs
    for(int j = 0; j < kel.dim(); ++j) {
      if(selfadjoint && dofs[i] > dofs[j]) continue; // Work with upper symmetric half
      if(dofs[j] < 0 || (l = unconstrNum[dofs[j]]) < 0) continue;  // Skip undefined/constrained dofs
      for(int m = xunonz[l]; m < xunonz[l+1]; ++m) {
        if(rowu[m] == k) {
          ScalarTypes::addComplex(unonz[m], std::complex<double>(0,kel[i][j]));
          break;
        }
      }
    }
  }
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::add(FullSquareMatrixC &kel, int *dofs)
{
  int k,l;
  for(int i = 0; i < kel.dim(); ++i) {
    if(dofs[i] < 0 || (k = unconstrNum[dofs[i]]) < 0) continue; // Skip undefined/constrained dofs
    for(int j = 0; j < kel.dim(); ++j) {
      if(selfadjoint && dofs[i] > dofs[j]) continue; // Work with upper symmetric half
      if(dofs[j] < 0 || (l = unconstrNum[dofs[j]]) < 0) continue;  // Skip undefined/constrained dofs
      for(int m = xunonz[l]; m < xunonz[l+1]; ++m) {
        if(rowu[m] == k) {
          ScalarTypes::addComplex(unonz[m], kel[i][j]);
          break;
        }
      }
    }
  }
} 

template<typename Scalar, typename SolverClass> 
double
GenEiSparseMatrix<Scalar,SolverClass>::getMemoryUsed()
{
  return sizeof(Scalar)*nnz/(1024.0*1024.0);
}

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::mult(const Scalar *_rhs, Scalar *_result)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(const_cast<Scalar*>(_rhs),numUncon,1), result(_result,numUncon,1);
  if(selfadjoint)
    result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
    result = M*rhs;
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::mult(const GenVector<Scalar> &_rhs, Scalar *_result) 
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1), result(_result,numUncon,1);
  if(selfadjoint)
    result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
    result = M*rhs;
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::transposeMult(const GenVector<Scalar> &_rhs, GenVector<Scalar> &_result)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1), result(_result.data(),numUncon,1);
  if(selfadjoint)
    result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
    result = M.adjoint()*rhs;
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::multAdd(const Scalar *_rhs, Scalar *_result)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(const_cast<Scalar*>(_rhs),numUncon,1), result(_result,numUncon,1);
  if(selfadjoint)
    result += M.template selfadjointView<Eigen::Upper>()*rhs;
  else
    result += M*rhs;
}

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::mult(const GenVector<Scalar> &_rhs, GenVector<Scalar> &_result)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1), result(_result.data(),numUncon,1);
  if(selfadjoint)
    result = M.template selfadjointView<Eigen::Upper>()*rhs;
  else
    result = M*rhs;
}

template<typename Scalar, typename SolverClass> 
Scalar
GenEiSparseMatrix<Scalar,SolverClass>::diag(int dof) const
{
  std::cerr << "GenEiSparseMatrix<Scalar,SolverClass>::diag is not implemented\n";
  return Scalar();
}

template<typename Scalar, typename SolverClass> 
Scalar &
GenEiSparseMatrix<Scalar,SolverClass>::diag( int dof )
{
  static Scalar defaultValue;
  std::cerr << "GenEiSparseMatrix<Scalar,SolverClass>::diag is not implemented\n";
  return defaultValue;
}

template<typename Scalar, typename SolverClass> 
long
GenEiSparseMatrix<Scalar,SolverClass>::size()
{
  return (numUncon) ? nnz : 0;
}

#ifdef DISTRIBUTED
#include <Comm.d/Communicator.h>
#endif

template<typename Scalar, typename SolverClass> 
void
GenEiSparseMatrix<Scalar,SolverClass>::unify(FSCommunicator *communicator)
{
#ifdef DISTRIBUTED
  communicator->globalSum(nnz, unonz);
#endif
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::factor()
{
  solver.compute(M);
  if(solver.info() != Eigen::Success) std::cerr << "sparse factor failed\n";
}

template<typename Scalar, typename SolverClass>
void 
GenEiSparseMatrix<Scalar,SolverClass>::solve(Scalar *_rhs, Scalar *_solution)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1), solution(_solution,numUncon,1);
  solution = solver.solve(rhs);
  if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::solve(GenVector<Scalar> &_rhs, GenVector<Scalar> &_solution)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1), solution(_solution.data(),numUncon,1);
  solution = solver.solve(rhs);
  if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::reSolve(Scalar *_rhs)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs,numUncon,1);
  rhs = solver.solve(rhs);
  if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
}

template<typename Scalar, typename SolverClass>
void
GenEiSparseMatrix<Scalar,SolverClass>::reSolve(GenVector<Scalar> &_rhs)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > rhs(_rhs.data(),numUncon,1);
  rhs = solver.solve(rhs);
  if(solver.info() != Eigen::Success) std::cerr << "sparse solve failed\n";
}
