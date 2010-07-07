#ifndef DIST_PARALLEL_SOLVER
#define DIST_PARALLEL_SOLVER

template<class Scalar> class GenDistrVector;

template<class Scalar> class GenDistrVectorSet;
class Timings;

template <class Scalar> 
class GenParallelSolver {
  public:
    virtual ~GenParallelSolver() = 0;
  //  virtual void factor()=0;
    virtual void reSolve(GenDistrVector<Scalar> &) = 0;
    virtual double getSolutionTime() = 0;
    virtual void solve(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &) = 0;
    virtual void getRBMs(Scalar *) = 0;
    virtual Timings& getTimers() = 0;
    virtual int numRBM()=0;
    virtual void getRBMs(GenDistrVectorSet<Scalar> &) = 0;
    virtual void reconstruct() {}
    virtual void refactor() {}
    virtual double getFNormSq(GenDistrVector<Scalar> &f) { cerr << "GenParallelSolver::getSqNorm not implemented\n"; }
};

template<class Scalar>
inline GenParallelSolver<Scalar>::~GenParallelSolver() { } // empty by default

#endif
