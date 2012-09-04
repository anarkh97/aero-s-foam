#ifndef _FETIHAXISOLVER_H_
#define _FETIHAXISOLVER_H_

class MDAxiData;
class FetiInfo;
class Connectivity;
class DofSetArray;
class compStruct;
class DistrInfo;
class DistrComplexVector;
template<class Scalar> class GenVector;
typedef GenVector<DComplex> ComplexVector;
template<class Scalar> class GenSolver;
typedef GenSolver<DComplex> ComplexSolver;
class GCRC;
template <class Scalar> class GenFullM;
typedef GenFullM<DComplex> FullMC;
class Timings;


class FetiHAxiSolver {

    MDAxiData **mdAxi;
    int nsub;
    int numModes;

    FetiInfo *fetiInfo;

    Connectivity *subToSub;
    Connectivity *edgeToSub, *subToEdge, *edgeToEdge;
    DofSetArray *coarseEqs; 
    compStruct renumber;

    int halfSize;
    DistrInfo internal, interface;

    double epsilon2;
    int maxiter;

    ComplexSolver **CoarseSolver;
    ComplexVector **CoarseVector;

    GCRC *oSet;

// The following might change into a FetiWorkSpace
    DistrComplexVector *deltaU;
    DistrComplexVector *deltaF;

    int MPCSize;
    FullMC **CoarseMPC;
    SkyMatrixC *MPCSolver;

    Timings times;

  public : 

    FetiHAxiSolver(int, MDAxiData **, Connectivity *, int, FetiInfo *);
    void computeLengths(int &, DistrInfo &, DistrInfo &);

    void defineMatrices(int, int **, int **);
    void constructAndAssembleMat(int, int *, int **, int **); 
    void deleteGlMap(int, int **, int **);
    void factorMatrices(int iP);

    void allocateBuffer(int nSub);
    void countEdges(int iSub, int *edges);
    void numberEdges(int iSub, int *eP, int* edges);
    void receiveNeighbEdgeNums(int iSub, int *eP, int* edges);
    void makeEdgeConnectivity();
    void makeCoarse();
    void allocateCoarseSolver(int, Connectivity *);
    void allocateCoarseVectors(int);
    void makeCoarseGridDofs(int);
    void assembleCoarse(int);
    void factorCoarseMatrices(int);
    void makeStaticLoad(DistrComplexVector &f);
    void makeLocalStaticLoad(int iSub, DistrComplexVector &f);
    DistrInfo &localInfo()  { return internal; }


    void solve(DistrComplexVector &, DistrComplexVector &);


    void computePtd(DistrComplexVector &, DistrComplexVector &);
    void zeroLocalCVec(int);
    void multQtBK(int, DistrComplexVector &);
    void assembleCoarseVector(int);
    void solveCoarsePb(int);
    void subtractQw(int, DistrComplexVector &);
    void multBK(int, DistrComplexVector &, DistrComplexVector &);


    void sendInterf(int, DistrComplexVector &, int &);
    void InterfDiff(int, DistrComplexVector &, int &);


    double preCondition(DistrComplexVector &, DistrComplexVector &);
    void multKbb(int, DistrComplexVector &, DistrComplexVector &,
                 DistrComplexVector &, DistrComplexVector &);
    void normDeltaF(int, DComplex *, DistrComplexVector *);


    void applyFP(DistrComplexVector &, DistrComplexVector &);
    void multTransposeBKQ(int, DistrComplexVector &);
    void finishFPz(int, DistrComplexVector &, DistrComplexVector &);


    void orthogonalize(DistrComplexVector &, DistrComplexVector &,
                       DistrComplexVector &, DistrComplexVector &);
    void gatherHalfInterface(int, DistrComplexVector *,
                      DistrComplexVector *, DComplex *, DComplex *); 
    void scatterHalfInterface(int, DComplex *, DistrComplexVector *,
                              int *);
    void rebuildInterface(int, DistrComplexVector &, int &);
    void orthoAdd(DistrComplexVector &, DistrComplexVector &, DComplex);


    void recoverU(DistrComplexVector &, DistrComplexVector &, 
                  DistrComplexVector &);
    void finishBtPz(int, DistrComplexVector &,DistrComplexVector &);
    void subtractBt(int, DistrComplexVector &, DistrComplexVector &);
    void localSolution(int, DistrComplexVector &,DistrComplexVector &);

    
    Timings& getTimers(); 
    double getSolutionTime();

    void solveMPC(DistrComplexVector &, ComplexVector &, 
                  DistrComplexVector &);
    void buildMPCSolver();
    void startSchur(int Fourier);
    void finishSchur(int Fourier);  
    void computePtd(DistrComplexVector &rl, DistrComplexVector &f, 
                    ComplexVector &g, ComplexVector &muVector, 
                    ComplexVector **lambdaVec);
    void assembleMPCVector(int, DComplex *, double *);
    void copyMu(int, DComplex *, DComplex *);
    void buildRHSSchur(int, DComplex *);
    void recoverSchur(int, ComplexVector **, ComplexVector *);
    void subtractQwMPC(int iSub, DistrComplexVector *f,
                    ComplexVector **lambdaVec, ComplexVector *mu);    
    void applyFPMPC(DistrComplexVector &, DistrComplexVector &, 
                    ComplexVector **lambdaVec);
    void addCKQ(int, ComplexVector **, ComplexVector *);
    void finishFzl(int, DistrComplexVector &, ComplexVector &,
                   ComplexVector **&, DistrComplexVector &);
    void recoverUMPC(DistrComplexVector &, DistrComplexVector &,
                  ComplexVector &, DistrComplexVector &, 
                  ComplexVector &muVector, ComplexVector **lambdaVec);
 
};

#endif
