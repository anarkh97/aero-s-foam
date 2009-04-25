#ifndef _FOURIER_DESCR_H_
#define _FOURIER_DESCR_H_


class Domain;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<DComplex> ComplexSparseMatrix;
template <class Scalar> class GenVector;
typedef GenVector<DComplex> ComplexVector;
template <class Scalar> class GenSolver;
typedef GenSolver<DComplex> ComplexSolver;
template <class Scalar> class GenCuCSparse;
typedef GenCuCSparse<DComplex> CuCComplexSparse;
template<class Scalar> class GenBLKSparseMatrix;
typedef GenBLKSparseMatrix<DComplex> BLKSparseMatrixC;
template<class Scalar> class GenDBSparseMatrix;
typedef GenDBSparseMatrix<DComplex> DBComplexSparseMatrix;
template <class Scalar> class GenSkyMatrix;
typedef GenSkyMatrix<DComplex> SkyMatrixC;
class FourierHelmBCs;


struct GenHelmParam {
   double kappaHelm;
};


class FourierStatic {
    Domain *domain;
    ComplexSolver *solver;
    CuCComplexSparse *kuc;
    FourierHelmBCs *bcs;
    ComplexVector bcxC;
    GenHelmParam *hParams;

 public:

    FourierStatic(Domain *d, FourierHelmBCs *);
    void make_bc(int *bc);
    ComplexVector make_bcxC(int mode);
    void preProcess();
    int solVecInfo();
    ComplexSolver *getSolver();
    int getModes();
    int getSlices();
    void constructK();
    void BuildKs(int mode);
    SkyMatrixC *constructSkyMatrixC();
    BLKSparseMatrixC *constructBLKSparseMatrixC();
    void makeComplexK(int mode, ComplexSparseMatrix *K);
    void buildRHS(ComplexVector &force, int mode);
    ComplexVector MergeSolBCs(ComplexVector sol);
    void postProcessing(); 
    void Reconstruction3D(ComplexVector *globalsol);
    Domain *getDomain() {return domain;};
    double getWave() {return (*hParams).kappaHelm;};

};

#endif
