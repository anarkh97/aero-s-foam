#ifndef _MD_NL_STATIC_H_
#define _MD_NL_STATIC_H_

template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
class Domain;
template <class Scalar> class GenParallelSolver;
typedef GenParallelSolver<double> ParallelSolver;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class Corotator;
class StaticTimers;
class DistrInfo;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
class DistrGeomState;
template <class Scalar> class GenMultiDomainPostProcessor;
typedef GenMultiDomainPostProcessor<double> MultiDomainPostProcessor;

class MDNLStatic 
{
    Domain     *domain;
    DecDomain  *decDomain;
    ParallelSolver *solver;

    FullSquareMatrix **kelArray;
    Corotator ***allCorot;

    double firstDv;
    double firstRes;
    double lastRes;
    double tolerance;
    StaticTimers *times;
    int numSystems;
    double deltaLambda;

    std::map<std::pair<int,int>, double> *mu; // lagrange multipliers for the contact surfaces
    std::vector<double> *lambda; // lagrange multipliers for all the other constraints

 public:
    // Constructor
    MDNLStatic(Domain *d);
    virtual ~MDNLStatic();

    DistrInfo& solVecInfo();
    DistrInfo& sysVecInfo();
    DistrInfo& elemVecInfo();
    int checkConvergence(int iter, double normDv, double normRes);
    int  getMaxit();
    double getScaleFactor();
    double getDeltaLambda0();
    double getMaxLambda();
    void getRHS(DistrVector &);
    ParallelSolver *getSolver();

    void printTimers();

    void staticOutput(DistrGeomState *geomState, double lambda, 
                      DistrVector &force, DistrVector &glRes, DistrGeomState *refState);

    MultiDomainPostProcessor *getPostProcessor();

    void preProcess();

    int reBuild(int iter, int step, DistrGeomState& geomState);

    DistrGeomState* createGeomState();

    void updatePrescribedDisplacement(DistrGeomState *geomState, double l=1.0);

    double getStiffAndForce(DistrGeomState& geomState, DistrVector& residual, 
                            DistrVector& elementInternalForce, DistrVector& gRes,
                            double lambda = 1.0, DistrGeomState *refState = NULL);

    double getTolerance() { return tolerance*firstRes; }

    bool linesearch();
    double getEnergy(double lambda, DistrVector& force, DistrGeomState* geomState)
      { cerr << "MDNLStatic::getEnergy is not implemented\n"; return 0; }

    double getResidualNorm(DistrVector &vec);

  private:
    void getSubStiffAndForce(int isub, DistrGeomState &geomState,
                             DistrVector &res, DistrVector &elemIntForce, double lambda,
                             DistrGeomState *refState);

    void makeSubCorotators(int isub);
    void makeSubKelArrays(int isub);
    void makeSubDofs(int isub);
    void updatePrescribedDisp(int isub, DistrGeomState& geomState);
    void subGetRHS(int isub, DistrVector& rhs);
    void addConstraintForces(int isub, DistrVector &vec);
    void getConstraintMultipliers(int isub);
    void updateConstraintTerms(DistrGeomState* geomState);
};

#endif
