#ifndef _MD_NL_STATIC_H_
#define _MD_NL_STATIC_H_

template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
class Domain;
template <class Scalar> class GenFetiSolver;
typedef GenFetiSolver<double> FetiSolver;
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
    FetiSolver *solver;

    FullSquareMatrix **kelArray;
    Corotator ***allCorot;

    double firstDv;
    double firstRes;
    double lastRes;
    double tolerance;
    StaticTimers *times;
    int numSystems;
    double deltaLambda;
 public:
    // Constructor
    MDNLStatic(Domain *d);

    DistrInfo& solVecInfo();
    DistrInfo& sysVecInfo();
    DistrInfo& elemVecInfo();
    int checkConvergence(int iter, double normDv, double normRes);
    int  getMaxit();
    double getScaleFactor();
    double getDeltaLambda0();
    double getMaxLambda();
    void getRHS(DistrVector &, DistrGeomState *gs=0);
    void addExternalForce(DistrVector &, DistrGeomState *gs, double lambda) { }
    FetiSolver *getSolver();

    void printTimers();

    void staticOutput(DistrGeomState *geomState, double lambda, 
                      DistrVector &force, DistrVector &glRes);

    MultiDomainPostProcessor *getPostProcessor();

    void preProcess();

    int reBuild(int iter, int step, DistrGeomState& geomState);

    DistrGeomState* createGeomState();

    void updatePrescribedDisplacement(DistrGeomState *geomState, double l=1.0);

    double getStiffAndForce(DistrGeomState& geomState, DistrVector& residual, 
                            DistrVector& elementInternalForce, DistrVector&gRes);

    void getSubStiffAndForce(int isub, DistrGeomState &geomState, 
                             DistrVector &res, DistrVector &elemIntForce);

    void makeSubCorotators(int isub);
    void makeSubKelArrays(int isub);
    void makeSubDofs(int isub);
    void updatePrescribedDisp(int isub, DistrGeomState& geomState);

    double getTolerance() {return tolerance*firstRes;}

    double getEnergy(double lambda, DistrVector& force, DistrGeomState* geomState)
      { cerr << "MDNLStatic::getEnergy is not implemented\n"; }

};

#endif
