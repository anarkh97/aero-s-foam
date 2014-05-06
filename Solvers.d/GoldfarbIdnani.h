#ifndef _GOLDFARB_IDNANI_H_
#define _GOLDFARB_IDNANI_H_

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#include <Math.d/BLKSparseMatrix.h>
#include <Solvers.d/eiquadprog.hpp>

//#define CHECK_G

template<class BaseSolver, class Scalar>
class GoldfarbIdnaniQpSolver : public BaseSolver
{
  int n, p, m;
  int *unconstrNum; // mapping from dsa to c_dsa unconstrained numbering (KKT)
  int *doftype;
  int *dofmap;

  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXd;
#if !defined(SPARSE_G) && !defined(CHECK_G)
  VectorXd diagG;
#else
  MatrixXd G;
#endif
  MatrixXd CE;
  MatrixXd CI;

  double tol;
  bool check;
  Scalar dummy;

public:
  template<class BaseArgs>
  GoldfarbIdnaniQpSolver(BaseArgs &ba, ConstrainedDSA *cdsa, double _tol, bool _check)
   : BaseSolver(ba) {
    tol = _tol;
    check = _check;
    unconstrNum = cdsa->getUnconstrNum();
    n = 0; p = 0; m = 0;
    doftype = new int[cdsa->size()];
    dofmap = new int[cdsa->size()];
    for(int i = 0; i < cdsa->numNodes(); ++i) {
      int dof;
      if((dof = cdsa->locate(i, DofSet::LagrangeE)) > -1) {
        doftype[dof] = 1;
        dofmap[dof] = p++;
      }
      else if((dof = cdsa->locate(i, DofSet::LagrangeI)) > -1) {
        doftype[dof] = 2;
        dofmap[dof] = m++;
      }
      else {
        int myFirstDof = cdsa->firstdof(i);
        int myNumDofs  = cdsa->weight(i);
        for(int j = 0; j < myNumDofs; ++j) {
          doftype[myFirstDof+j] = 0;
          dofmap[myFirstDof+j] = n++;
        }
      }
    }
#if !defined(SPARSE_G) && !defined(CHECK_G)
    diagG.resize(n); diagG.setZero(); // in this case G stores only the diagonal
#else
    G.resize(n,n); G.setZero();
#endif
    CE.resize(n,p); CE.setZero();
    CI.resize(n,m); CI.setZero();
  }
  ~GoldfarbIdnaniQpSolver() {
     delete [] doftype;
     delete [] dofmap;
  }
  void add(FullSquareMatrix &kel, int *dofs) {
    int I,J;
    for(int i = 0; i < kel.numRow(); ++i ) {
      if((I = unconstrNum[dofs[i]]) < 0 || doftype[I] != 0) continue;
      for(int j = 0; j < kel.numCol(); ++j) {
        if((J = unconstrNum[dofs[j]]) < 0) continue;
        switch(doftype[J]) {
          case 0: 
#if !defined(SPARSE_G) && !defined(CHECK_G)
            if(I==J) diagG[dofmap[I]] += kel[i][j];
#else
            G(dofmap[I],dofmap[J]) += kel[i][j];
#endif
            break;
          case 1:
            CE(dofmap[I],dofmap[J]) += kel[i][j];
            break;
          case 2:
            CI(dofmap[I],dofmap[J]) -= kel[i][j]; // changing sign because inequality constraints 
                                                  // for this solver need to be CE^T*x + ce0 >= 0
            break;
        }
      }
    }
    BaseSolver::add(kel,dofs);
  }
  void addDiscreteMass(int dof, Scalar s) {
    int I;
    if((I = unconstrNum[dof]) < 0 || doftype[I] != 0) return;
#if !defined(SPARSE_G) && !defined(CHECK_G)
    diagG[dofmap[I]] += s;
#else
    G(dofmap[I],dofmap[I]) += s;
#endif
    BaseSolver::addDiscreteMass(dof,s);
  }
  void zeroAll() { 
#if !defined(SPARSE_G) && !defined(CHECK_G)
    diagG.setZero();
#else
    G.setZero(); 
#endif
    CE.setZero();
    CI.setZero();
    BaseSolver::zeroAll();
  }
  int  dim() { return n+p+m; }
  Scalar diag(int dof) const { return Scalar(0); } // TODO
  Scalar &diag(int dof) { return dummy; } // TODO

  void solve(Scalar* rhs, Scalar* sol); 
  void solve(GenVector<Scalar> &rhs, GenVector<Scalar> &sol) {
    solve(rhs.data(), sol.data());
  }
  int neqs() { return n+p+m; }
  long size() { return 0; } // TODO
  void factor() { 
    // XXXX currently we are factoring inside eiquadprog.hpp
    //BaseSolver::factor();
  }
  void reSolve(Scalar *rhs) {  
    Scalar *rhs_copy = new Scalar[neqs()];
    for(int i=0; i<neqs(); ++i) rhs_copy[i] = rhs[i];
    solve(rhs_copy, rhs);
    delete [] rhs_copy;
  }
  void reSolve(GenVector<Scalar> &rhs) {
    reSolve(rhs.data());
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<>
void
GoldfarbIdnaniQpSolver<WrapSparseMat<double>,double>::solve(double*, double*);

template<>
void
GoldfarbIdnaniQpSolver<WrapSparseMat<complex<double> >,complex<double> >::solve(complex<double>*, complex<double>*);

template<>
void
GoldfarbIdnaniQpSolver<WrapEiSparseMat<double>,double>::solve(double*, double*);

template<>
void
GoldfarbIdnaniQpSolver<WrapEiSparseMat<complex<double> >,complex<double> >::solve(complex<double>*, complex<double>*);
#endif

#endif
