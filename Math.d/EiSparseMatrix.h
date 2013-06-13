#ifndef EI_SPARSEMATRIX_H_
#define EI_SPARSEMATRIX_H_

#ifdef USE_EIGEN3
#include <complex>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#ifdef EIGEN_SPARSELU_SUPPORT
#include <Eigen/SparseLU>
#endif
#ifdef EIGEN_CHOLMOD_SUPPORT
#include <Eigen/CholmodSupport>
#endif
#ifdef EIGEN_UMFPACK_SUPPORT
#include <Eigen/UmfPackSupport>
#endif
#ifdef EIGEN_SUPERLU_SUPPORT
#include <Eigen/SuperLUSupport>
#endif
#ifdef EIGEN_SRQR_SUPPORT
#include <Eigen/SPQRSupport>
#endif
#ifdef EIGEN_SPARSEQR_SUPPORT
#include <Eigen/SparseQR>
#endif
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>

template<class Scalar, class SolverClass = Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >
class GenEiSparseMatrix : public SparseData, public GenSparseMatrix<Scalar>, public GenSolver<Scalar>
{
 protected:
   // this is a symmetric sparse matrix using CSR storage (upper triangluar part only in self-adjoint case) 
   // and Eigen 3 implementation via MappedSparseMatrix
   bool selfadjoint;
   int nnz;
   Scalar *unonz;
   Eigen::MappedSparseMatrix<Scalar, Eigen::ColMajor, int> M;
   SolverClass solver;

 public:
   GenEiSparseMatrix(Connectivity *, DofSetArray *, DofSetArray *, bool=true);
   GenEiSparseMatrix(Connectivity *, DofSetArray *, int *, bool=true);
   GenEiSparseMatrix(Connectivity *, EqNumberer *);
   virtual ~GenEiSparseMatrix();

   Eigen::MappedSparseMatrix<Scalar, Eigen::ColMajor, int>& getEigenSparse() { return M; }

   // GenSparseMatrix assembly
   void add(FullSquareMatrix &, int *dofs);
   void add(GenAssembledFullM<Scalar> &, int *);
   void addImaginary(FullSquareMatrix &, int *dofs);
   void add(FullSquareMatrixC &, int *dofs); 
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
   void printSparse(const std::string& filename);
};


template<class Scalar, class SolverClass = Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>,Eigen::Upper> >
class WrapEiSparseMat : public GenEiSparseMatrix<Scalar, SolverClass>
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

    WrapEiSparseMat(CtorData &ctd) : GenEiSparseMatrix<Scalar,SolverClass>(ctd.cn, ctd.dsa, ctd.cdsa, ctd.flg) {}
};


#ifdef _TEMPLATE_FIX_
#include <Math.d/EiSparseMatrix.C>
#endif

#endif
#endif
