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

  // allocate
  VectorXd wA(maxvec), rhs(maxvec), vk(maxvec), ylar(maxvec), t(maxvec), sign(maxvec);
  VectorXd crlt(A.cols()), h(A.cols()), S(A.cols()), minBuffer(A.cols());
  VectorXd colK(A.rows()), mu(A.rows()), residual(A.rows()), update(A.rows());
 
  // intitialize
  wA.setZero(); rhs.setOnes(); vk.setZero(); ylar.setZero(); t.setZero(); sign.setZero();
  minBuffer.setZero();
  colK.setZero(); mu.setZero(); 

  info = 1;

  MatrixXd B(A.rows(),maxvec), R(maxvec,maxvec);
  B.setZero(); R.setZero();

  std::vector<long int> indices;
  std::vector<long int> nld_indices;

  if(scaling) {
    for(int i=0; i<A.cols(); ++i) { 
      VectorXd dummy = A.col(i).array();
      double s       = dummy.norm(); 
      S[i]           = (s != 0) ? 1/s : 0; 
    }
  } else {
    S.setOnes();
  }
  double  bnorm = b.norm();
  double abstol = reltol*bnorm;
  rnorm = bnorm;
  double C;
  //initialize starting set
  {
    crlt = S.asDiagonal()*(A.transpose()*b);
    long int i;
    if(positive){
      C = crlt.maxCoeff(&i);
    } else {
      C = crlt.cwiseAbs().maxCoeff(&i);
    }
    sign[0] = sgn(crlt[i]);
    indices.push_back(i);
    B.col(0)     = S[indices[0]]*(A.col(indices[0]).array());
    double diagK = B.col(0).squaredNorm();
    double r     = sqrt(diagK);
    R(0,0)       = r;
  }

  bool     dropId  = false;
  long int blockId = 0; 
  long int k       = 0; // current column index <-> Active set size -1
  long int iter    = 0; // number of iterations
  long int downIt  = 0; // number of downdates
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

    residual = b - mu; 

    crlt = S.asDiagonal()*(A.transpose()*residual); // correlations

    long int i = 0; // dummy integer for utilities

    // compute right hand side and solve the symetric system of equations
    wA.head(k+1)  = R.topLeftCorner(k+1,k+1).triangularView<Upper>().transpose().solve(rhs.head(k+1));
    wA.head(k+1)  = R.topLeftCorner(k+1,k+1).triangularView<Upper>().solve(wA.head(k+1));
    double oneNwA = sqrt(rhs.head(k+1).dot(wA.head(k+1)));
    wA.head(k+1) *= 1.0/oneNwA;

    // compute smallest angle at which a new covarient becomes dominant
    update = B.leftCols(k+1)*wA.head(k+1);
    h = S.asDiagonal()*(A.transpose()*update); 
    minBuffer = (C - crlt.array())/(1.0/oneNwA - h.array());
    for(std::vector<long int>::iterator it = indices.begin(); it != indices.end(); ++it) minBuffer[*it] = std::numeric_limits<double>::max();
    for(long int row = 0; row < minBuffer.rows(); ++row) minBuffer[row] = (minBuffer[row] <= 0) ? std::numeric_limits<double>::max() : minBuffer[row];
    if(dropId) minBuffer[blockId] = std::numeric_limits<double>::max();
    double gamma1 = minBuffer.minCoeff(&i); 

    long int j;
    if(!positive) {
      minBuffer = (C + crlt.array())/(1.0/oneNwA + h.array());
      for(std::vector<long int>::iterator it = indices.begin(); it != indices.end(); ++it) minBuffer[*it] = std::numeric_limits<double>::max();
      for(long int row = 0; row < minBuffer.rows(); ++row) minBuffer[row] = (minBuffer[row] <= 0) ? std::numeric_limits<double>::max() : minBuffer[row];
      if(dropId) minBuffer[blockId] = std::numeric_limits<double>::max();
      double gamma2 = minBuffer.minCoeff(&j);
      if(gamma2 < gamma1){
        gamma1 = gamma2;
        i = j;
      }
    }

    // compute smallest angle at which ylars changes sign
    t.head(k+1) = -1.0*ylar.head(k+1).array()/(sign.head(k).array()*wA.head(k+1).array());
    for(long int ele = 0; ele < k+1; ++ele) t[ele] = (t[ele] <= 0) ? std::numeric_limits<double>::max() : t[ele];
    double gamma_tilde = (k > 0) ? t.head(k+1).minCoeff(&j) : std::numeric_limits<double>::max();

    dropId = false; 
    if(gamma_tilde < gamma1){
      dropId = true;
      gamma1 = gamma_tilde;
      i = j; // drop index if gamma_tilde selected
      ylar.head(k+1) = ylar.head(k+1).array() + gamma1*(sign.head(k+1).array()*wA.head(k+1).array());
      mu += gamma1*update;
      blockId = indices[i];
    } else {
      indices.push_back(i); // add index if gamma1 selected
       
      ylar.head(k+1) = ylar.head(k+1).array() + gamma1*(sign.head(k+1).array()*wA.head(k+1).array());
      mu += gamma1*update;  

      B.col(k+1) = S[indices[k+1]]*(A.col(indices[k+1]).array());
      double Correlation = B.col(k+1).transpose()*(b-mu);
      sign[k+1] = sgn(Correlation);
      B.col(k+1) *= sign[k+1];

      //update cholesky factorization R'*R = B'*B where R is upper triangular
      double diagK   = B.col(k+1).squaredNorm();
      colK.head(k+1) = B.leftCols(k+1).transpose()*(B.col(k+1));
      vk.head(k+1)   = R.topLeftCorner(k+1,k+1).triangularView<Upper>().transpose().solve(colK.head(k+1));
      double r       = sqrt(diagK - vk.head(k+1).squaredNorm());

      Block<MatrixXd,Dynamic,1,true> colR = R.col(k+1);
      colR.head(k+1) = vk.head(k+1);
      colR[k+1]      = r;

    }
  
    //update solution and estimate
    C  -= gamma1/oneNwA;

    k++;
    iter++;
    if(dropId) {// if gamma_tilde is selected, remove that index and downdate Cholesky factorization
      downIt++;

      long int leftSide = k;
      k = i;

      //remove selected index
      std::vector<long int>::iterator fol = indices.erase(indices.begin()+k);
       
      //zero out column
      R.col(k).head(k+1).setZero();
      sign[k] = 0;
      ylar[k] = 0;

      //now zero out diagonal
      for(std::vector<long int>::iterator it = fol; it != indices.end(); ++it,k++){
        //construct Givens rotation
        double l = R.block<2,1>(k,k+1).norm();
        double c = R(k,k+1)/l;
        double s = -1.0*R(k+1,k+1)/l;
        Matrix2d G; G << c, -s, s, c;
        //apply to block of R to zero out old diagonal element
        R.block(k,k+1,2,leftSide-k-1) = G*R.block(k,k+1,2,leftSide-k-1);
        R.col(k) = R.col(k+1);
        B.col(k) = B.col(k+1);
        ylar[k] = ylar[k+1]; 
        sign[k] = sign[k+1];
      }
      R.col(k).head(k+1).setZero();
      R.row(k).head(k+1).setZero();
      ylar[k] = 0;
      sign[k] = 0; 
      k--;
    }

    rnorm = residual.norm();
  }

  if(verbose) std::cout.flush();

  VectorXd x = VectorXd::Zero(A.cols());
  for(long int j=0; j<k; ++j) x[indices[j]] = S[indices[j]]*ylar[j];
  return x;
}

#endif
