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

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
};

Eigen::VectorXd
lars(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm,
      long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool positive)
{
  using namespace Eigen;
  
  const long int m = b.rows();
  const long int maxvec = std::min(m, (long int)(maxsze*A.cols()));
  const long int maxit = maxite*A.cols();

  VectorXd wA(maxvec), rhs(maxvec), vk(maxvec), ylar(maxvec), t(maxvec), sign(maxvec);
  VectorXd g(A.cols()), h(A.cols()), S(A.cols()), Avg(A.cols()), minBuffer(A.cols());
  VectorXd colK(A.rows()), mu(A.rows()), residual(A.rows()), update(A.rows()), standB(A.rows());
  sign.setZero();
  mu.setZero();
  vk.setZero();
  rhs.setOnes();
  colK.setZero(); 
  ylar.setZero();
  wA.setZero();
  t.setZero();
  minBuffer.setZero();
  info = 1;
  MatrixXd B(A.rows(),maxvec), R(maxvec,maxvec);
  B.setZero(); R.setZero();
  long int k = 0;
  std::vector<long int> indices;
  std::vector<long int> nld_indices;

  if(scaling) {
    for(int i=0; i<A.cols(); ++i) { 
      Avg[i] = A.col(i).sum()/A.rows();
      VectorXd dummy = A.col(i).array()-Avg[i];
      double s = dummy.squaredNorm(); 
      S[i] = (s != 0) ? 1/s : 0; 
    }
    standB = b.array() - b.sum()/A.rows();
  } else {
    S.setOnes();
    Avg.setZero();
    standB = b; 
  }
  double  bnorm = standB.norm();
  double abstol = reltol*bnorm;
  rnorm = bnorm;

  int iter   = 0; // number of iterations
  int downIt = 0; // number of downdates
  while(true) {

    if(verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Downdate = " << std::setw(9) << downIt << "    "
                << "Active set size = " << std::setw(9) << k << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if(rnorm <= abstol || k+nld_indices.size() == maxvec) break;
    if(iter >= maxit) { info = 3; break; }

    residual = standB - mu; 

    g = S.asDiagonal()*(A.transpose()*residual - residual.sum()*Avg); // gradient
    long int i;
    h = g; 
    double gi;
    if(positive){
      for(std::vector<long int>::iterator it = indices.begin(); it != indices.end(); ++it) h[*it] = -std::numeric_limits<double>::max(); // make sure the index has not already been selected
      for(std::vector<long int>::iterator it = nld_indices.begin(); it != nld_indices.end(); ++it) h[*it] = -std::numeric_limits<double>::max(); // also make sure near linear dependent indices are not selected
      gi = h.maxCoeff(&i);
      sign[k] = 1;
    } else {
      for(std::vector<long int>::iterator it = indices.begin(); it != indices.end(); ++it) h[*it] = 0.0; // make sure the index has not already been selected
      for(std::vector<long int>::iterator it = nld_indices.begin(); it != nld_indices.end(); ++it) h[*it] = 0.0; 
      gi = h.cwiseAbs().maxCoeff(&i);
      sign[k] = sgn(h[i]);
    }

    if(gi <= 0) break;
    B.col(k) = sign[k]*S[i]*(A.col(i).array()-Avg[i]);
    indices.push_back(i);
   
    //update cholesky factorization R'*R = B'*B where R is upper triangular
    double diagK = B.col(k).squaredNorm();
    colK.head(k) = B.leftCols(k).transpose()*B.col(k);
    vk.head(k)   = R.topLeftCorner(k,k).triangularView<Upper>().transpose().solve(colK.head(k));
    double r     = sqrt(diagK - vk.head(k).squaredNorm());

    Block<MatrixXd,Dynamic,1,true> colR = R.col(k);
    colR.head(k) = vk.head(k);
    colR[k]      = r;
    if(r <= 0) { nld_indices.push_back(i); indices.pop_back(); continue; } else nld_indices.clear(); // check for near linear dependence

    // compute right hand side and solve the symetric system of equations
    wA.head(k+1)  = R.topLeftCorner(k+1,k+1).triangularView<Upper>().solve(R.topLeftCorner(k+1,k+1).triangularView<Upper>().transpose().solve(rhs.head(k+1)));
    double oneNwA = sqrt(rhs.head(k+1).dot(wA.head(k+1)));
    wA.head(k+1)  = wA.head(k+1).array()/oneNwA;
    std::cout << "wA = " << wA.head(k+1).transpose() << std::endl;

    update = B.leftCols(k+1)*wA.head(k+1);
    h = S.asDiagonal()*(A.transpose()*update - update.sum()*Avg); 
    minBuffer = (gi - g.array())/(1.0/oneNwA - h.array());
    for(std::vector<long int>::iterator it = indices.begin(); it != indices.end(); ++it) minBuffer[*it] = std::numeric_limits<double>::max();
    for(long int row = 0; row < minBuffer.rows(); ++row) minBuffer[row] = (minBuffer[row] <= 0) ? std::numeric_limits<double>::max() : minBuffer[row];
    double gamma1 = minBuffer.minCoeff(); 

    if(!positive) {
      minBuffer = (gi + g.array())/(1.0/oneNwA + h.array());
      for(std::vector<long int>::iterator it = indices.begin(); it != indices.end(); ++it) minBuffer[*it] = std::numeric_limits<double>::max();
      for(long int row = 0; row < minBuffer.rows(); ++row) minBuffer[row] = (minBuffer[row] <= 0) ? std::numeric_limits<double>::max() : minBuffer[row];
      double gamma2 = minBuffer.minCoeff();
      gamma1 = std::min(gamma1,gamma2);
    }

    t.head(k) = -1.0*ylar.head(k).array()/(sign.head(k).array()*wA.head(k).array());
    for(long int ele = 0; ele < k; ++ele) t[ele] = (t[ele] <= 0) ? std::numeric_limits<double>::max() : t[ele];
    double gamma_tilde = (k > 0) ? t.head(k).minCoeff(&i) : std::numeric_limits<double>::max();

    bool dropId = false; 
    if(gamma_tilde < gamma1){
      dropId = true;
      gamma1 = gamma_tilde;
    }
  
    //update solution and estimate
    std::cout << "ylar = " << ylar.head(k+1).transpose() << std::endl;
    ylar.head(k+1) = ylar.head(k+1).array() + gamma1*(sign.head(k+1).array()*wA.head(k+1).array()); 
    mu += gamma1*update;

    std::cout << "ylar = " << ylar.head(k+1).transpose() << std::endl;
    std::cout << "sign = " << sign.head(k+1).transpose() << std::endl;

    k++;
    iter++;
    if(dropId) {// if gamma_tilde is selected, remove that index and downdate Cholesky factorization
      downIt++;

      long int leftSide = k;
      k = i;

      //remove selected index
      std::vector<long int>::iterator fol = indices.erase(indices.begin()+k);
      
      //zero out column
      R.col(k).setZero();
      sign[k] = 0; 
      ylar[k] = 0;

      //now zero out diagonal
      for(std::vector<long int>::iterator it = fol; it != indices.end(); ++it){
        //construct Givens rotation
        double l = R.block<2,1>(k,k+1).norm();
        double c = R(k,k+1)/l;
        double s = -1.0*R(k+1,k+1)/l;
        Matrix2d G; G << c, -s, s, c;
        //apply to block of R
        R.block(k,k+1,2,leftSide-k-1) = G*R.block(k,k+1,2,leftSide-k-1);
        R.col(k) = R.col(k+1);
        B.col(k) = B.col(k+1);
        sign[k] = sign[k+1];
        ylar[k] = ylar[k+1]; 
        k++;
      }
      R.col(k).setZero();
      R.row(k).setZero();
      sign[k] = 0; 
      ylar[k] = 0;
    }

    rnorm = residual.norm();
  }

  VectorXd final = B.leftCols(k)*ylar.head(k) - standB; 
  std::cout << "******* final residual = " << final.norm() << " **************"<< std::endl;

  if(verbose) std::cout.flush();

  VectorXd x = VectorXd::Zero(A.cols());
  for(long int j=0; j<k; ++j) x[indices[j]] = S[indices[j]]*ylar[j];
  return x;
}

#endif
