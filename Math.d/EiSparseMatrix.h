#ifndef EI_SPARSEMATRIX_H_
#define EI_SPARSEMATRIX_H_

#include <complex>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>

template<class Scalar>
class GenEiSparseMatrix : public SparseData, public GenSparseMatrix<Scalar>, public GenSolver<Scalar>
{
   // this is a symmetric sparse matrix using CSR storage (upper triangluar part only) 
   // and Eigen 3 implementation via MappedSparseMatrix and Simplicial Cholesky
   int nnz;
   Scalar *unonz;
   Eigen::MappedSparseMatrix<Scalar, Eigen::ColMajor, int> M;
   Eigen::SimplicialCholesky<Eigen::SparseMatrix<Scalar>,Eigen::Upper> llt;

 public:
   GenEiSparseMatrix(Connectivity *, DofSetArray *, DofSetArray *);
   GenEiSparseMatrix(Connectivity *, DofSetArray *, int *);
   GenEiSparseMatrix(Connectivity *, EqNumberer *);
   virtual ~GenEiSparseMatrix();

   Eigen::MappedSparseMatrix<Scalar, Eigen::ColMajor, int> getEigenSparse() { return M; }

   // GenSparseMatrix assembly
   void add(FullSquareMatrix &, int *dofs);
   void addDiscreteMass(int dof, Scalar diMass);
   void add(int, int, Scalar);

   // GenSparseMatrix matrix-vector multiplications
   void mult(const GenVector<Scalar> &, GenVector<Scalar> &);
   void mult(const GenVector<Scalar> &rhs, Scalar *result);
   void mult(const Scalar *, Scalar *);
   void multAdd(const Scalar *, Scalar *);
   void transposeMult(const GenVector<Scalar> & rhs, GenVector<Scalar> & result);

   // GenSparseMatrix miscellaneous
   void zeroAll();
   int dim() { return numUncon; }
   double getMemoryUsed();
   int numRow() { return numUncon; }
   int numCol() { return numUncon; }
   Scalar diag(int dof) const;
   Scalar &diag(int dof);

   // GenSolver factor
   void factor();

   // GenSolver solve
   void solve(Scalar *rhs, Scalar *solution);
   void solve(GenVector<Scalar> &rhs, GenVector<Scalar> &solution );
   void reSolve(Scalar *rhs);
   void reSolve(GenVector<Scalar> &rhs);
/* TODO
   void reSolve(int nRHS, Scalar **RHS);
   void reSolve(int nRHS, GenVector<Scalar> * RHS);
   void reSolve(GenFullM<Scalar> *);
*/
   // GenSolver miscellaneous
   int  neqs() { return numUncon; }
   long size();
   void unify(FSCommunicator *communicator);
   void print();
};

typedef GenEiSparseMatrix<double> EiSparseMatrix;
typedef GenEiSparseMatrix<std::complex<double> > EiComplexSparseMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/EiSparseMatrix.C>
#endif

#endif
