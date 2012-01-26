#ifdef USE_EIGEN3
#include <iostream>
#include <Solvers.d/GoldfarbIdnani.h>
#include <Solvers.d/eiquadprog.hpp>
using std::cerr; using std::endl;

template<>
void
GoldfarbIdnaniQpSolver<WrapSparseMat<double>,double>::solve(double* _rhs, double* _sol)
{
#ifndef SPARSE_G
  Eigen::Map<VectorXd> rhs(_rhs,n+p+m), sol(_sol,n+p+m);

  VectorXd g0(n), ce0(p), ci0(m);
  for(int i = 0; i < neqs(); ++i) {
    switch(doftype[i]) {
      case 0 : g0(dofmap[i])  = -rhs(i); break;
      case 1 : ce0(dofmap[i]) = -rhs(i); break;
      case 2 : ci0(dofmap[i]) =  rhs(i); break;
    }
  }
  VectorXd x(n);
  VectorXd lambda(p), mu(m);
  double f = Eigen::solve_quadprog(G, g0, CE, ce0, CI, ci0, x, &lambda, &mu);
  for(int i = 0; i < neqs(); ++i) {
    switch(doftype[i]) {
      case 0 : sol(i) = x(dofmap[i]); break;
      case 1 : sol(i) = -lambda(dofmap[i]); break;
      case 2 : sol(i) = mu(dofmap[i]); break;
    }
  }
#endif
}

template<>
void
GoldfarbIdnaniQpSolver<WrapSparseMat<std::complex<double> >,std::complex<double> >::solve(std::complex<double>*, std::complex<double>*)
{
  cerr << "GoldfarbIdnaniQpSolver not implemented for complex\n";
}

template<>
void
GoldfarbIdnaniQpSolver<WrapEiSparseMat<double>,double>::solve(double* _rhs, double* _sol)
{
  Eigen::Map<VectorXd> rhs(_rhs,n+p+m), sol(_sol,n+p+m);

  VectorXd g0(n), ce0(p), ci0(m);
  for(int i = 0; i < neqs(); ++i) {
    switch(doftype[i]) {
      case 0 : g0(dofmap[i])  = -rhs(i); break;
      case 1 : ce0(dofmap[i]) = -rhs(i); break;
      case 2 : ci0(dofmap[i]) =  rhs(i); break;
    }
  }
  VectorXd x(n);
  VectorXd lambda(p), mu(m);
  try {
#ifdef SPARSE_G
    Eigen::SparseMatrix<double> SparseG = EiSparseMatrix::getEigenSparse();
    double f = Eigen::solve_quadprog(SparseG, G.trace(), g0, CE, ce0, CI, ci0, x, &lambda, &mu, tol);
#else
    double f = Eigen::solve_quadprog(G, g0, CE, ce0, CI, ci0, x, &lambda, &mu, tol);
#endif
  }
  catch(std::runtime_error& e) {
    cerr << "exception: " << e.what() << endl;
  }
  for(int i = 0; i < neqs(); ++i) {
    switch(doftype[i]) {
      case 0 : sol(i) = x(dofmap[i]); break;
      case 1 : sol(i) = -lambda(dofmap[i]); break;
      case 2 : sol(i) = mu(dofmap[i]); break;
    }
  }
}

template<>
void
GoldfarbIdnaniQpSolver<WrapEiSparseMat<std::complex<double> >,std::complex<double> >::solve(std::complex<double>*, std::complex<double>*)
{
  cerr << "GoldfarbIdnaniQpSolver not implemented for complex\n";
}
#endif
