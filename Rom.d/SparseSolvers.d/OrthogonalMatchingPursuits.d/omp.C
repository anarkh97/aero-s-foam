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

Eigen::VectorXd
nncgp(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm,
      long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, double &dtime);

Eigen::VectorXd
omp(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm,
      long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool project, double &dtime)
{
  using namespace Eigen;

  const long int m = b.rows();
  const long int maxvec = std::min(m, (long int)(maxsze*A.cols()));
  const long int maxit = maxite*A.cols();
  double bnorm = b.norm();
  double abstol = reltol*bnorm;

  VectorXd x_(maxvec), y(maxvec), r(A.rows()), g(A.cols()), h(A.cols()), DtGDinv(maxvec), g_(maxvec), a(maxvec), S(A.cols());
  x_.setZero();
  r = b;
  rnorm = bnorm;
  info = 1;
  MatrixXd B(A.rows(),maxvec), D(maxvec,maxvec), GD(maxvec,maxvec);
  Matrix<double,Dynamic,Dynamic,ColMajor> BD(A.rows(),maxvec);
  B.setZero(); D.setZero(); BD.setZero(); GD.setZero();
  long int k = 0;
  std::vector<long int> indices;
  std::vector<long int> nld_indices;

  if(scaling) for(int i=0; i<A.cols(); ++i) { double s = A.col(i).norm(); S[i] = (s != 0) ? 1/s : 0; }
  else S.setOnes();

  int iter   = 0; // number of iterations
  int downIt = 0; // number of downdates
  while(true) {

    if(verbose) {
      std::cout << "Iteration = "       << std::setw(9)  << iter   << "    "
                << "Downdate = "        << std::setw(9)  << downIt << "    "
                << "Active set size = " << std::setw(9)  << k      << "    "
                << "Residual norm = "   << std::setw(13) << std::scientific << std::uppercase << rnorm  << "    "
                << "Target = "          << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if(rnorm <= abstol || k+nld_indices.size() == maxvec) break;
    if(iter >= maxit) { info = 3; break; }

    g = S.asDiagonal()*(A.transpose()*r); // gradient
    long int i;
    h = g; for(long int j=0; j<k; ++j) h[indices[j]] = -std::numeric_limits<double>::max(); // make sure the index has not already been selected
    for(long int j=0; j<nld_indices.size(); ++j) h[nld_indices[j]] = -std::numeric_limits<double>::max(); // also make sure near linear dependent indices are not selected
    double gi = h.maxCoeff(&i);
    if(gi <= 0) break;
    B.col(k) = S[i]*A.col(i);
    indices.push_back(i);
    // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
    GD.row(k).head(k) = B.col(k).transpose()*BD.leftCols(k);

    for(long int j=0; j<k+1; ++j) g_[j] = g[indices[j]];
    Block<MatrixXd,Dynamic,1,true> d = D.col(k), c = BD.col(k);
    d.head(k+1) = g_.head(k+1) - D.topLeftCorner(k+1,k).triangularView<Upper>()*(DtGDinv.head(k).asDiagonal()*(GD.topLeftCorner(k+1,k).transpose()*g_.head(k+1))); // direction

    c = B.leftCols(k+1)*d.head(k+1);
    GD.col(k).head(k+1) = B.leftCols(k+1).transpose()*c;
    DtGDinv[k] = 1/c.squaredNorm();
    a[k] = r.dot(c)*DtGDinv[k]; // step length
    if(a[k] < 0) { nld_indices.push_back(i); indices.pop_back(); continue; } else nld_indices.clear(); // check for near linear dependence
    y.head(k+1) = x_.head(k+1) + a[k]*d.head(k+1); // candidate solution
    r -= a[k]*c; // residual
    k++;
    iter++;
    x_.head(k) = y.head(k);

    rnorm = r.norm();
  }

  if(project){
    std::cout << "*** PROJECTING SOLUTION ONTO SELECTED BASIS ***" << std::endl;
    x_.head(k) = nncgp(B.leftCols(k), b, rnorm, info, maxsze, maxite, reltol, verbose, scaling, dtime);
  }

  if(verbose) std::cout.flush();

  VectorXd x = VectorXd::Zero(A.cols());
  for(long int j=0; j<k; ++j) x[indices[j]] = S[indices[j]]*x_[j];
  return x;
}

#endif
