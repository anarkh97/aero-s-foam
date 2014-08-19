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

bool operator== (const long_int& lhs, const long_int& rhs);

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pgpfp(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       double maxsze, double reltol, bool verbose, bool scaling, bool positive)
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

  const int nsub = A.size(); // number of subdomains assigned to this mpi process
  const long int m = b.rows();
  const long int maxvec = std::min(m, (long int)(maxsze*n));
  const long int maxit = 15*n;
  std::vector<long int> maxlocvec(nsub); for(int i=0; i<nsub; ++i) maxlocvec[i] = std::min(A[i].cols(),maxvec);
  double bnorm = b.norm();
  double abstol = reltol*bnorm;

  Array<VectorXd,Dynamic,1> y(nsub), g(nsub), g1(nsub), g2(nsub), g_(nsub), S(nsub);
  for(int i=0; i<nsub; ++i) {
    y[i] = VectorXd::Zero(maxlocvec[i]);
    g_[i].resize(maxlocvec[i]);
    S[i].resize(A[i].cols());
    if(scaling) for(int j=0; j<A[i].cols(); ++j) S[i][j] = 1/A[i].col(j).norm();
    else S[i].setOnes();
  }
  VectorXd r(m), vertex(m), DtGDinv(maxvec), z(maxvec), a(maxvec), lambda(maxvec);
  r = b;
  vertex.setZero();
  rnorm = bnorm;
  a.setZero();
  lambda.setZero();
  Array<MatrixXd,Dynamic,1> B(nsub), D(nsub), GD(nsub);
  for(int i=0; i<nsub; ++i) {
    B[i]  = MatrixXd::Zero(m,maxlocvec[i]);
    D[i]  = MatrixXd::Zero(maxlocvec[i],maxvec);
    GD[i] = MatrixXd::Zero(maxlocvec[i],maxvec);
  }
  Matrix<double,Dynamic,Dynamic,ColMajor> BD(m,maxvec);
  MatrixXd z_loc(maxvec,nsub), c_loc(m,nsub);
  long int k = 0; // k is the dimension of the (global) set of selected indices
  std::vector<long int> l(nsub); // l[i] is the dimension of the subset of selected indices local to a subdomain.
  std::list<std::pair<int,long_int> > gindices; // global indices
  std::vector<std::vector<long int> > indices(nsub); // local indices
  std::vector<std::vector<int> > setKey(nsub);

  Array<long int,Dynamic,1> jk(nsub);
  ArrayXd lamMax(nsub), setMax(nsub), ymin(nsub);

  long int iter   = 0; // number of iterations
  long int downIt = 0;
  while(true) {

    if(myrank == 0 && verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Downdate = " << std::setw(9) << downIt << "    "
                << "Active set size = " << std::setw(9) << k << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if(rnorm <= abstol || k == maxvec || iter >= maxit) break;

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      g[i]  = S[i].asDiagonal()*(A[i].transpose()*r);
      g2[i] = S[i].asDiagonal()*(A[i].transpose()*vertex);

      g1[i] = g[i];
      // compute maximum step length for this subdomain
      for(long int col = 0; col != A[i].cols(); col++){
        double num = g1[i][col];
        double den = g2[i][col];
        if(num >= 0. && den != 1.0) {
          g1[i][col] = num/(1.0-den);
          if(!positive)
            g2[i][col] = -std::numeric_limits<double>::max();
        } else if (num < 0. && den != -1.0) {
          g1[i][col] = -std::numeric_limits<double>::max();
          if(!positive)
            g2[i][col] = (-1.0)*num/(1.0+den);
        } else {
          g1[i][col] = -std::numeric_limits<double>::max();
          if(!positive)
            g2[i][col] = -std::numeric_limits<double>::max();
        }
      }
      // make sure that element has not already been selected
      for(long int j=0; j<l[i]; ++j){
        g1[i][indices[i][j]] = -std::numeric_limits<double>::max();
        if(!positive)
          g2[i][indices[i][j]] = -std::numeric_limits<double>::max();
      }

      bool whichSet = 1;
      long int position1 = 0;
      long int position2 = 0;

      double lam  = 0.0;
      double lam1 = g1[i].maxCoeff(&position1);
      double lam2 = 0.;
      if(!positive)
        lam2 = g2[i].maxCoeff(&position2);

      if(A[i].cols() > 0){
        if (lam1 > lam2 || positive){
          lam = lam1;
          jk[i] = position1;
        } else {
          lam = lam2;
          jk[i] = position2;
          whichSet = -1;
        }
      } else {
        lam = -std::numeric_limits<double>::max();
      }

      lamMax[i] = lam;
      setMax[i] = whichSet;
    }
    int ik; // subdomain which has the max coeff.
    int Set; // correct sign of selected set
    s.val   = (nsub > 0) ? lamMax.maxCoeff(&ik) : -std::numeric_limits<double>::max();
    s.rank  = myrank;
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpicomm);
#endif
    if(s.rank == myrank) {
      p.index = jk[ik];
      p.sub   = ik;
      Set = setMax[ik];
    }
#ifdef USE_MPI
    MPI_Bcast(&p, 1, MPI_LONG_INT, s.rank, mpicomm);
    MPI_Bcast(&Set, 1, MPI_INT, s.rank, mpicomm);
#endif
  
     
    lambda[k] = 1.0/s.val;

    if(s.rank == myrank) {
      B[ik].col(l[ik]) = double(Set)*S[ik][jk[ik]]*A[ik].col(jk[ik]);
      // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
      GD[ik].row(l[ik]).head(k) = B[ik].col(l[ik]).transpose()*BD.leftCols(k);
      indices[ik].push_back(jk[ik]);
      setKey[ik].push_back(Set);
      l[ik]++;
    }
    gindices.push_back(std::pair<int,long_int>(s.rank,p));

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      for(long int j=0; j<l[i]; ++j) g_[i][j] = double(setKey[i][j])*g[i][indices[i][j]];
      z_loc.col(i).head(k) = GD[i].topLeftCorner(l[i],k).transpose()*g_[i].head(l[i]);
    }
    z.head(k) = z_loc.topRows(k).rowwise().sum();
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, z.data(), k, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif

    z.head(k) = (DtGDinv.head(k).asDiagonal()*z.head(k)).eval();
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
      d.head(l[i]) = g_[i].head(l[i]) - D[i].topLeftCorner(l[i],k)*z.head(k);
      c_loc.col(i) = B[i].leftCols(l[i])*d.head(l[i]);
    }
    Block<MatrixXd,Dynamic,1,true> c = BD.col(k) = c_loc.rowwise().sum();
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, c.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif
    DtGDinv[k] = 1/c.squaredNorm();
    a[k] = r.dot(c)*DtGDinv[k];
    r -= a[k]*c;
    vertex += lambda[k]*a[k]*c;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
      y[i].head(l[i]) += a[k]*d.head(l[i]);
      GD[i].col(k).head(l[i]) = B[i].leftCols(l[i]).transpose()*c;
    }
    k++;

    while(true) {
      iter++;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
      for(int i=0; i<nsub; ++i) {
        ymin[i] = (l[i] > 0) ? y[i].head(l[i]).minCoeff(&jk[i]) : std::numeric_limits<double>::max();
      }
      int ik;
      s.val = (nsub > 0) ? ymin.minCoeff(&ik) : std::numeric_limits<double>::max();
      s.rank = myrank;
#ifdef USE_MPI
      MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE_INT, MPI_MINLOC, mpicomm);
#endif

      if(s.val < 0.) {
        downIt++;

        if(s.rank == myrank) {
          p.index = indices[ik][jk[ik]];
          p.sub = ik;
          indices[ik].erase(indices[ik].begin() + jk[ik]); // remove index jk[ik] from the local active set of ik-th subdomain
          setKey[ik].erase(setKey[ik].begin() + jk[ik]);
          for(int j=jk[ik]; j<l[ik]-1; ++j) { y[ik][j] = y[ik][j+1]; } // erase jk[ik]-th element from x_[ik] and y[ik]
          y[ik][l[ik]] = 0;
          l[ik]--;
          //std::cout << "removing index " << p.index << " from subdomain " << ik << " on process with rank " << myrank << std::endl;
        }

#ifdef USE_MPI
        MPI_Bcast(&p, 1, MPI_LONG_INT, s.rank, mpicomm);
#endif
        // remove index from the global active set
        std::list<std::pair<int,long_int> >::iterator pos = std::find(gindices.begin(), gindices.end(), std::pair<int,long_int>(s.rank,p));
        std::list<std::pair<int,long_int> >::iterator fol = gindices.erase(pos);

        // Note: it is necessary to re-G-orthogonalize the basis D now, project the solution y onto the new basis and compute the corresponding residual r.
        // This is done here by starting from the column of D pointed to by fol (because the ones before this are already G-orthogonal), and then
        // following what is the essentially same procedure that is used above to construct the original basis, with a few optimizations when possible.
        k = std::distance(gindices.begin(), fol);
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) {
          y[i].head(l[i]).setZero();
          l[i] = 0; for(std::list<std::pair<int,long_int> >::iterator it = gindices.begin(); it != fol; ++it) { if(it->first == myrank && it->second.sub == i) l[i]++; }
          D[i].bottomRightCorner(maxlocvec[i]-l[i],maxvec-k).setZero(); // XXX this is a larger block than necessary
          y[i].head(l[i]) = D[i].topLeftCorner(l[i],k)*a.head(k);
        }
        r = b - BD.leftCols(k)*a.head(k);
        vertex = BD.leftCols(k)*lambda.head(k).asDiagonal()*a.head(k);
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) {
          g_[i].head(l[i]) = B[i].leftCols(l[i]).transpose()*r;
        }

        for(std::list<std::pair<int,long_int> >::iterator it = fol; it != gindices.end(); ++it) {
          double num = 0.;
          double den = 0.;
          if(it->first == myrank) {
            int i = it->second.sub;
            int slot;
            for(slot = 0; slot != indices[i].size(); slot++){
              if(it->second.index == indices[i][slot])
                break;
            }
            B[i].col(l[i]) = double(setKey[i][slot])*S[i][it->second.index]*A[i].col(it->second.index);
            num = B[i].col(l[i]).transpose()*r;
            den = B[i].col(l[i]).dot(vertex);
            g_[i][l[i]] = num;
            // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
            GD[i].row(l[i]).head(k) = B[i].col(l[i]).transpose()*BD.leftCols(k);
            l[i]++;
          }

#ifdef USE_MPI
          MPI_Bcast(&num, 1, MPI_DOUBLE, it->first, mpicomm);
          MPI_Bcast(&den, 1, MPI_DOUBLE, it->first, mpicomm);
#endif

          lambda[k] = (1.0-den)/num;

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {
            z_loc.col(i).head(k) = GD[i].topLeftCorner(l[i],k).transpose()*g_[i].head(l[i]);
          }
          z.head(k) = z_loc.topRows(k).rowwise().sum();
#ifdef USE_MPI
          MPI_Allreduce(MPI_IN_PLACE, z.data(), k, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif

          z.head(k) = (DtGDinv.head(k).asDiagonal()*z.head(k)).eval();
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {
            Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
            d.head(l[i]) = g_[i].head(l[i]) - D[i].topLeftCorner(l[i],k)*z.head(k);
            c_loc.col(i) = B[i].leftCols(l[i])*d.head(l[i]);
          }
          Block<MatrixXd,Dynamic,1,true> c = BD.col(k) = c_loc.rowwise().sum();
#ifdef USE_MPI
          MPI_Allreduce(MPI_IN_PLACE, c.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif

          DtGDinv[k] = 1/c.squaredNorm();
          a[k] = r.dot(c)*DtGDinv[k];
          r -= a[k]*c;
          vertex += lambda[k]*a[k]*c;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {
            Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
            y[i].head(l[i]) += a[k]*d.head(l[i]);
            GD[i].col(k).head(l[i]) = B[i].leftCols(l[i]).transpose()*c;
            g_[i].head(l[i]) -= a[k]*GD[i].col(k).head(l[i]);
          }
          k++;
        }
      }
      else {
        break;
      }
    }

    rnorm = r.norm();
  }

  if(myrank == 0 && verbose) std::cout.flush();

  Array<VectorXd,Dynamic,1> x(nsub);
  for(int i=0; i<nsub; ++i) {
    x[i] = VectorXd::Zero(A[i].cols());
    for(long int j=0; j<l[i]; ++j) x[i][indices[i][j]] += S[i][indices[i][j]]*double(setKey[i][j])*y[i][j];
  }
  return x;
}

#endif
