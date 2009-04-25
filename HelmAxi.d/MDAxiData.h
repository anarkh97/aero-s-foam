#ifndef _MDAXIDATA_H_
#define _MDAXIDATA_H_
#include <Math.d/SymFullMatrix.h>

template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
class FourierHelmBCs;
template<class Scalar> class GenVector;
typedef GenVector<DComplex> ComplexVector;
template<class Scalar> class GenCuCSparse;
typedef GenCuCSparse<DComplex> CuCComplexSparse;
template<class Scalar> class GenDBSparseMatrix;
typedef GenDBSparseMatrix<DComplex> DBComplexSparseMatrix;
template<class Scalar> class GenSolver;
typedef GenSolver<DComplex> ComplexSolver;
class Polygon;
class LineAxiSommer;
template <class Scalar> class GenSparseSet;
typedef GenSparseSet<DComplex> ComplexSparseSet;
template <class Scalar> class GenSymFullMatrix;
typedef GenSymFullMatrix<DComplex> SymFullMatrixC;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> GenFullSquareMatrixC;
class DofSetArray;
class Connectivity;
class MPCData;
template <class Scalar> class GenFullM;
typedef GenFullM<DComplex> FullMC;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
typedef GenSparseMatrix<DComplex> ComplexSparseMatrix;

struct GenHelmParam {
   double kappaHelm;
};

class MDAxiData : public SubDomain {

    GenHelmParam *hParams;

    ComplexSolver **KC;
    CuCComplexSparse **KucC;
    DBComplexSparseMatrix **KbbC;   // For preconditioning
    CuCComplexSparse **KibC;
    ComplexSolver **KiiC;
    DBComplexSparseMatrix **KccC;   // To get normal derivative on body

    int numInterface;
    SommerElement **InterfElem;
    int* InterfaceIndex;    

    int numInQSet;
    ComplexSparseSet **QSet;
    DComplex **BKQ;

    ComplexVector **localCVec;

    public :

    MPCData *locMPCs;
    DComplex **CKQ;
    SymFullMatrixC *CKCt;
    int *localToGlobalMPC;
    int *globalToLocalMPC;

    int halfOffset;

    FourierHelmBCs *locBCs;
    int *locToglDir;

    int InterfSign;
    DComplex *interfBuff;

    int *coarseEqNums;
    int numCGridDofs;

    MDAxiData(Domain *d, int n, Connectivity *cn, Connectivity *sToN);

    void extractBCs(FourierHelmBCs *glBCs, Connectivity *nodeToSub);
    void extractMPCs(MPCData *glMPCs, Connectivity *mpcToSub,
                     Connectivity *nodeToSub);

    int* getInterfaceIndex();

    void renumberData();

    void getInterface(PolygonSet ***allPs); 
    void finishInterface(PolygonSet ***allPs);

    void sendDOFList();
    void gatherDOFList(DofSet ***allBoundary);

    void defineKs(int **, int **);

    void makeKs(int *, int *, int m1=0, int m2=-1);
    void Assemble(int m1=0, int m2=-1);
    void addInterface(int *,int m1=0, int m2=-1);

    void deleteGlMap(int **, int **);

    void factorKC(int m1=0, int m2=-1);
    void factorKiiC(int m1=0, int m2=-1);

    void prepareCoarseData(DofSet ***, int *, int **, int **, int **, 
                           DComplex ***);
    void makeCoarseData(int *, int **, int **, int **, DComplex ***, int f1=0,
                        int f2=-1);
    void deleteCoarseInfo(int **, int **, int **, DComplex ***);
    DComplex CoarseGridValue(int, double, double, DComplex);
    void setNumCoarseDofs();
    int numCoarseDofs();
    void makeCKCt();
    void prepareMPCSet(int f1=0, int f2=-1);

    void assembleCoarse(int, ComplexSparseMatrix *, FullMC **); 
    void computeQtKQ(int, int, int, DComplex **, DComplex *);
    void assembleQtKQ(int, int, int, DComplex *, ComplexSparseMatrix *);
    void assembleCKQSet(int, int, int, DComplex *, FullMC **); 
    void deleteConn();
    int* getCoarseGridDofs(DofSetArray *, Connectivity &);
    void makeCoarseGridDofs(DofSetArray *, Connectivity *);

    void buildRHS(DComplex *, int f1=0, int f2=-1);

    void assembleDUDN(DComplex *u, DComplex **dudn, int f1=0, int f2=-1);
    void localDUDN(DComplex **u, DComplex **dudn, int *, int *, int f1=0,
                   int f2=-1);

    DComplex *getLocalCVec(int);
    void zeroLocalCVec(int f1=0, int f2=-1);

    void multLocalQtBK(DComplex *, int f1=0, int f2=-1);
    void multKf(DComplex *, DComplex *, int shiftu=0, int shiftf=0, int f1=0,
                int f2=-1);
    void computeBw(DComplex *, DComplex *, int f1=0, int f2=-1);
    void sendInterf(DComplex *);
    void interfaceJump(DComplex *);
    void multKbbC(DComplex *, DComplex *, DComplex *, DComplex *, int f1=0,
                  int f2=-1);
    void sendDeltaF(DComplex *);
    DComplex collectAndDotDeltaF(DComplex *);
    void multTransposeBKQ(DComplex *, int f1=0, int f2=-1);
    void finishFPz(ComplexVector **, DComplex *, DComplex *, int f1=0, 
                   int f2=-1);
    void getHalfInterf(DComplex *, DComplex *); 
    void getHalfInterf(DComplex *, DComplex *, DComplex *, DComplex *);
    void scatterHalfInterf(DComplex *, DComplex *);
    void rebuildInterf(DComplex *);
    void subtractBt(DComplex *, DComplex *, int f1=0, int f2=-1);
    void finishBtPz(ComplexVector **, DComplex *,DComplex *);

    double getWaveNumber() { return hParams->kappaHelm; }
    int getMPCSize();

    void mergeSolution(DComplex **, DComplex *);

    void subtractQw(DComplex *f, ComplexVector **w, ComplexVector *mu=0,
                       int f1=0, int f2=-1);
    void addCKQw(ComplexVector **, ComplexVector *, int f1=0, int f2=-1);
    void finishFzl(DComplex *, DComplex *, ComplexVector **, DComplex *, 
                   int f1=0, int f2=-1);
    
};

#endif
