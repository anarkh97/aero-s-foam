#ifdef USE_SCALAPACK
#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <Timers.d/GetTime.h>
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

#include "Plh.h"

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
splh(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool project, double &dtime) {
  // Setup
  int mypid;
  //MPI_Comm_rank(MPI_COMM_WORLD, &mypid);
  //std::cout << "mypid " << mypid << " Entering splh." << std::endl;

  Plh solver(A);
  solver.setRtol(reltol);
  solver.init(A, b);
  solver.summary();

  // Solve
  int max_iter = 5000;
  int nfree = solver.solve(max_iter);

  // Output
  solver.printTimes(true);
  rnorm = solver.getResidualNorm();
  dtime = solver.getComputeTime() + solver.getDistributeMatrixTime();
  Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1> x(n);
  return solver.getSolution();
}

#endif
#endif
