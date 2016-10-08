#ifndef DIST_PARALLEL_SOLVER
#define DIST_PARALLEL_SOLVER

template<class Scalar> class GenDistrVector;

template<class Scalar> class GenDistrVectorSet;
class Timings;

template <class Scalar> 
class GenParallelSolver {
  public:
    virtual ~GenParallelSolver() = 0;

    virtual void reSolve(GenDistrVector<Scalar> &) = 0;
    virtual double getSolutionTime() = 0;
    virtual void solve(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &) = 0;
    virtual void squareRootMult(GenDistrVector<Scalar> &) {
      std::cerr << "GenParallelSolver::squareRootMult(GenDistrVector<Scalar> &) is not implemented\n"; }
    virtual void inverseSquareRootMult(GenDistrVector<Scalar> &) {
      std::cerr << "GenParallelSolver::inverseSquareRootMult(GenDistrVector<Scalar> &) is not implemented\n"; }
    virtual void getRBMs(double *) {
      std::cerr << "GenParallelSolver::getRBMs(double *) is not implemented\n"; }
    virtual Timings& getTimers() = 0;
    virtual int numRBM() {
      std::cerr << "GenParallelSolver::numRBM() is not implemented\n";
      return 0; } 
    virtual void getRBMs(GenDistrVectorSet<double> &) {
      std::cerr << "GenParallelSolver::getRBMs(DistrVectorSet &) is not implemented\n"; }
    virtual void getRBMs(GenDistrVector<double> *rbms) {
      std::cerr << "GenParallelSolver::getRBMs(DistrVector *) is not implemented\n"; }

    virtual void reconstruct() {
      // Nothing to do
    }
    virtual void refactor() { 
      // Nothing to do
    }
    virtual double getFNormSq(GenDistrVector<Scalar> &f);
};

template <class Scalar>
inline
GenParallelSolver<Scalar>::~GenParallelSolver() {
  // Empty pure virtual destructor must be defined
}

#include <iostream>

template <class Scalar>
inline
double
GenParallelSolver<Scalar>::getFNormSq(GenDistrVector<Scalar> &) {
  std::cerr << "GenParallelSolver::getSqNorm not implemented\n";
  return 0.0;
}

#endif
