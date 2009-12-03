#ifndef _NON_LIN_STATIC_H_
#define _NON_LIN_STATIC_H_

#include <Problems.d/StaticDescr.h>

class Domain;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class StaticTimers;
class GeomState;
class Corotator;

class NonLinStatic {
    Domain *domain;
    double *bcx;
    Solver *solver;
    SparseMatrix *spm;
    Solver *prec;
    SparseMatrix *spp;
    FullSquareMatrix *kelArray;
    Corotator **allCorot;
    double firstRes;
    double firstDv;
    double tolerance;
    StaticTimers *times;
 public:
    // Constructor
    NonLinStatic(Domain *d) { domain = d; }

    int  solVecInfo();
    int  sysVecInfo();
    int  elemVecInfo();
    int  getMaxit();
    double getScaleFactor();  // only nlstatic
    double getDeltaLambda0(); // only nlstatic
    double getMaxLambda();    // only maxlambda
    void getRHS(Vector &rhs, GeomState *gs=0); 
    void preProcess();
    Solver *getSolver();
    SingleDomainPostProcessor<double,Vector,Solver> *getPostProcessor();

    int reBuild(int iter, int step, GeomState& geomState);
    GeomState* createGeomState();

    void staticOutput(GeomState *geomState, double lambda, Vector& force, Vector &);
    int checkConvergence(int iter, double normDv, double residualNorm);

    double getStiffAndForce(GeomState& geomState, Vector& residual, 
                          Vector& elementInternalForce, Vector &);

    void updatePrescribedDisplacement(GeomState *geomState, double lambda = 1);

    void printTimers();

    double getTolerance(){return tolerance*firstRes;}
};

#endif
