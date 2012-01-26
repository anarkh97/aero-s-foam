#ifndef EI_SPARSEMATRIX_H_
#define EI_SPARSEMATRIX_H_

#include <complex>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#ifdef EIGEN_CHOLMOD_SUPPORT
#include <Eigen/CholmodSupport>
#endif
#ifdef EIGEN_SUPERLU_SUPPORT
#include <Eigen/SuperLUSupport>
#endif
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>

template<class Scalar>
class GenEiSparseMatrix : public SparseData, public GenSparseMatrix<Scalar>, public GenSolver<Scalar>
{
   // this is a symmetric sparse matrix using CSR storage (upper triangluar part only in self-adjoint case) 
   // and Eigen 3 implementation via MappedSparseMatrix and Cholmod (self-adjoint case) or SuperLU backend solvers.
   // If EIGEN_CHOLMOD_SUPPORT is not defined then the Eigen implementation of the Simplicial Cholesky will be
   // used in the self-adjoint case
   bool selfadjoint;
   int nnz;
   Scalar *unonz;
   Eigen::MappedSparseMatrix<Scalar, Eigen::ColMajor, int> M;
#ifdef EIGEN_CHOLMOD_SUPPORT
   Eigen::CholmodDecomposition<Eigen::SparseMatrix<Scalar>,Eigen::Upper> llt;
#else
   Eigen::SimplicialCholesky<Eigen::SparseMatrix<Scalar>,Eigen::Upper> llt;
#endif
#ifdef EIGEN_SUPERLU_SUPPORT
   Eigen::SuperLU<Eigen::SparseMatrix<Scalar> > *lu;
#endif

 public:
   GenEiSparseMatrix(Connectivity *, DofSetArray *, DofSetArray *, bool=true);
   GenEiSparseMatrix(Connectivity *, DofSetArray *, int *, bool=true);
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

template<class Scalar>
class WrapEiSparseMat : public GenEiSparseMatrix<Scalar>
{
  public:
    struct CtorData {
      Connectivity *cn;
      DofSetArray *dsa;
      DofSetArray *cdsa;
      bool flg;
      CtorData(Connectivity *c, DofSetArray *d, DofSetArray *dc, bool f=true) {
        cn = c;
        dsa = d;
        cdsa = dc;
        flg = f;
      }
    };

    WrapEiSparseMat(CtorData &ctd) : GenEiSparseMatrix<Scalar>(ctd.cn, ctd.dsa, ctd.cdsa, ctd.flg) {}
};


#ifdef _TEMPLATE_FIX_
#include <Math.d/EiSparseMatrix.C>
#endif

#endif
