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
  //int mypid;
  //MPI_Comm_rank(MPI_COMM_WORLD, &mypid);
  //std::cout << "mypid " << mypid << ", maxsze = " << maxsze << std::endl;

  // Instantiate the solver and set some parameters
  Plh solver(A);
  int max_np = 2000;
  solver.setMaxNP(max_np);
  solver.setRtol(reltol);
  solver.setMaxIterRatio(maxite);
  if (scaling) solver.setColumnScaling();

  //int mb=64;
  //int nb=64;
  //solver.setABlockSize(mb,nb);

  // Best to keep the default of mp=NROCS and np=1
  //int mp=12;
  //int np=4;
  //solver.setAProcGrid(mp,np);

  // Loads the matrix and RHS
  solver.init(A, b);
  solver.summary();

  // Solve
  int nfree = solver.solve();

  // Output
  solver.printTimes(true);
  rnorm = solver.getResidualNorm();
  dtime = solver.getDownDateTime();
  Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1> x(n);
  return solver.getSolution();
}

#endif
#endif
