#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <cstdio>
#include <Utils.d/MyComplex.h>
#include <iostream>
#include <Timers.d/Timing.h>

template <class Scalar> class GenVector;
template <class Scalar> class GenDistrVector;
template <class Scalar> class DistrBlockVector;
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> FullSquareMatrixC;
template <class Scalar> class GenFullM;
class Rbm;
class GeomState;
class Connectivity;
class FSCommunicator;

template<class Scalar> 
class GenSolver {
  protected:
    mutable Timings times;
    mutable double solveTime; // to store solution time
    mutable double memUsed;
    bool print_nullity;
  public:
    GenSolver() { solveTime = 0.0; memUsed = 0.0; print_nullity = true; }
    virtual ~GenSolver() = 0;

    // Rbm functions
    virtual int numRBM();

    int numRBM() const { return const_cast<GenSolver *>(this)->numRBM(); }// TODO Make base const.

    // Functions to return the rigid body modes
    virtual void getRBMs(double *rbms);
    virtual void getRBMs(Vector *rbms);
    virtual void getRBMs(VectorSet &rbms);

    // Get number of equations
    virtual int neqs() const =0;
    virtual int dim() const;

    // Solve functions take two arguments, returning the result in the second
    virtual void solve(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution);
    virtual void solve(const Scalar *rhs, Scalar *solution);
    virtual void solve(const GenDistrVector<Scalar> &rhs, GenDistrVector<Scalar> &solution)
             {std::cerr << "GenSolver::solve(GenDistrVector<Scalar> NOT implemented" << std::endl; } 
    virtual void solve(const DistrBlockVector<Scalar> &rhs, DistrBlockVector<Scalar> &solution)
             {std::cerr << "GenSolver::solve(DistrBlockVector<Scalar> NOT implemented" << std::endl; } 

    // reSolve functions overwrite the rhs vector with the solution
    virtual void reSolve(Scalar *rhs);
    virtual void reSolve(Vector &rhs);
    virtual void reSolve(ComplexVector &rhs);

    // Multiple rhs reSolve functions
    virtual void reSolve(int nRHS, Scalar **rhs);
    virtual void reSolve(int nRHS, Scalar  *rhs);
    virtual void reSolve(int nRHS, GenVector<Scalar> *rhs);
    virtual void reSolve(GenFullM<Scalar> *);

    // Forward substitution for 1 rhs
    virtual void forward(GenVector<Scalar> &rhs);
    virtual void forward(Scalar *rhs);
    // Backward substitution for 1 rhs
    virtual void backward(GenVector<Scalar> &rhs);
    virtual void backward(Scalar *rhs);

    // Multiplication by upper triangular factor of an LU or Cholesky factorization
    virtual void upperMult(Scalar *rhs);

    // Multiplication by lower triangular factor of an LU or Cholesky factorization
    virtual void lowerMult(Scalar *rhs);

    // reBuild functions zero the appropriate structures
    // and reassemble the stiffness matrix and factor
    virtual void reBuild(FullSquareMatrix *kel, int iter=0, int step = 1);
    virtual void reBuild(FullSquareMatrix *kel,FullSquareMatrix *mel,
                         Scalar delta);
    virtual void reBuildGeometricRbms(GeomState *gs);
    virtual void factor();
    virtual void parallelFactor();

    // function to return solution time
    virtual double getSolutionTime() { return solveTime; }
    Timings& getTimers() { return times; }

    // function to return memory used
    virtual double getMemoryUsed() const { return memUsed; }

    virtual long size() const = 0;
    virtual void clean_up();

    // Thuan added functions
    virtual void zeroAll()  { fprintf(stderr, " ... This solver has no zeroAll()\n"); }
    virtual void add(FullSquareMatrix &, int *) 
            { fprintf(stderr, " ... This solver has no add(FullSquareMatrix &, int *)\n"); }
    virtual void addImaginary(FullSquareMatrix &, int *)
            { fprintf(stderr, " ... This solver has no addImaginary(FullSquareMatrix &, int *)\n"); }
    virtual void add(FullSquareMatrixC &kel, int *dofs) // RT added to support PML, DGM
            { fprintf(stderr, " ... This solver has no add(FullSquareMatrixC &kel, int *dofs)\n"); }
      
    virtual Connectivity *getAllDofs()  
            { fprintf(stderr, " ... This solver has no getAllDofs()\n"); return 0; }
    virtual void addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier);
    virtual void addone(Scalar d, int dofi, int dofj);
    virtual Scalar getone(int dofi, int dofj);
    virtual void unify(FSCommunicator *communicator);
    virtual void add(Scalar *d);
    virtual Scalar* getData();
    virtual void print() { };

    virtual void getNullSpace(Scalar *rbm)
            { fprintf(stderr, " ... This solver has no getNullSpace(Scalar *)\n"); }

    void setPrintNullity(bool _p) { print_nullity = _p; }

    virtual void apply(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution) { solve(rhs, solution); }
    virtual double getResidualNorm(const GenVector<Scalar> &v) { return v.norm(); }
};

template<class Scalar>
inline GenSolver<Scalar>::~GenSolver()    { } // empty by default

typedef GenSolver<double> Solver;
typedef GenSolver<DComplex> ComplexSolver;

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/Solver.C>
#endif

#endif
