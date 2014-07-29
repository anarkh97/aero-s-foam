#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <algorithm>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <list>
#include <utility>
#include <vector>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

struct double_int {
  double val;
  int rank;
};

struct long_int {
  long index;
  int sub;
};

// AUTHOR: Todd Chapman, Stanford University 2014
// This function solves the Compressed Sensing problem
//
// min ||x||_0 s.t. A*x = b
//
// which is equivalent to 
//
// min ||s||_1 s.t. [A,-A]*s = b & s >= 0 
//
// and its dual problem
//
// min b^T*c s.t. [A,-A]^T*c <= 1
//
// where s_i   = x_i > 0 for i = 1:N
// and   s_2*i = x_i < 0 for i = 1:N
//
// The solution framework is known as Basis Pursuit and the algorithm is called Gradient Polytope 
// Faces Pursuit. This is a paraticular implementation of PFP which uses the Conjugate Gradient method to 
// solve a system of equations at each iteration who's dimensionality changes from one iteration to the 
// next, thus only matrix vector products are required and the scheme is easily parallelizable.
// This algorithm produces identical residuals to that of the standard Polytope Faces Pursuit
//
// References:
// 1. Plumbley, M & Gretsistas, I. "Gradient Polytope Faces Pursuit for Large Sparse Recovery Problems" ICASSP
// 2. Blumensath, T & Davies, M. "Gradient Pursuits" IEEE
//
// AURGUMENTS:
// A       = where A*x= b
// b       = target vector 
// rnorm   = residual norm 
// maxsze  = maximum allowable sparsity
// reltol  = stopping criteria
// scaling = flag to turn on/off unit normalization of columns of A
//
Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pgpfp(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       double maxsze, double reltol, bool verbose, bool scaling)
{
  // each A[i] is the columnwise block of the global A matrix assigned to a subdomain on this mpi process
  // each x[i] of the return value x is the corresponding row-wise block of the global solution vector
  // n is the number of columns in the global A matrix (note this is only used to define the stopping criteria)
  using namespace Eigen;
  int myrank;
#ifdef USE_MPI
  const MPI_Comm mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank(mpicomm, &myrank);
#else
  myrank = 0;
#endif
  struct double_int s;
  struct long_int p;
 


}                   

#endif
