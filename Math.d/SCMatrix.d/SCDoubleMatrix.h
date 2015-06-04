#ifndef SCDOUBLEMATRIX_H_
#define SCDOUBLEMATRIX_H_

#include <mpi.h>
#include <vector>
#include <Eigen/Core>

#ifdef NNLS_DEV
#include "SCBaseMatrix.h"
#include "SCIntMatrix.h"
#else
#include "Math.d/SCMatrix.d/SCBaseMatrix.h"
#include "Math.d/SCMatrix.d/SCIntMatrix.h"
#endif

#define ZERO_TOL 0.0

#define NEGATIVE_INF -1.0e307
#define POSITIVE_INF  1.0e307

#define SCDBL_TIME_GETMAXLOC 0
#define SCDBL_TIME_GETMAX 1
#define SCDBL_TIME_GETMINLOC 2
#define SCDBL_TIME_GETMIN 3
#define SCDBL_N_TIMES 4

typedef struct {
    double x;
    int i;
} DoubleInt;


class SCDoubleMatrix : public SCBaseMatrix {

    public:
        SCDoubleMatrix(int context, int m, int n, int mb, int nb);
        SCDoubleMatrix(const SCDoubleMatrix& matrix);
        SCDoubleMatrix(const SCDoubleMatrix& matrix, int ncols);
        SCDoubleMatrix(std::string filename, int context, int mb, int nb);
        ~SCDoubleMatrix();

        void setA(int i, int j, double val);
        void setColA(int j, double *val);
        void write(std::string fname, bool compact=false, int m=0, int n=0);
        int pivot(int *ip, int *desc_ip);
        int setMatrixRow(int i, double *row);
        int setMatrixColumn(int j, double *col);
        int getMatrixRow(int i, double *row, char scope);
        int getMatrixColumn(int j, double *col, char scope);
        double norm2(int n=0);
        int multiply(SCDoubleMatrix &x, SCDoubleMatrix &y, char trans, int m, int n, double alpha, double beta);
        int multiply(char trans, int m, int n, double alpha, int ia, int ja,
                     SCDoubleMatrix &x, int ix, int jx, int incx, double beta,
                     SCDoubleMatrix &y, int iy, int jy, int incy);
        int multiply(SCDoubleMatrix &B, SCDoubleMatrix &C, char transA='N', char transB='N', double alpha=1.0, double beta=0.0, int m=0, int n=0, int k=0);
        int hadamardProduct(SCDoubleMatrix &x);
        int zero();
        void zero(int ix, int jx, int ni, int nj);
        int set(double val);
        void set(double val, int ix, int jx, int ni, int nj);
        int permuteOld(char direc, SCIntMatrix &ip, int m=0, int n=0);
        int permute(char direc, char rowcol, SCIntMatrix &ip, int m=0, int n=0);
        double * getMatrix() {return _matrix;};
        double * getTau() {return _tau->getMatrix();};
        double getElement(int i, int j);
        double getElement(int i, int j, int rsrc, int csrc);
        void setElement(int i, int j, double value);
        void setElementsLocal(const SCDoubleMatrix& matrix);
        void project(double value=0.0);
        int loadIdentityMatrix(double value=1.0);
        int copy(SCDoubleMatrix& matrix, int n);
        int copy(SCDoubleMatrix& dest, int n, SCIntMatrix& order);
        int add(SCDoubleMatrix& matrix, char trans, int n, double a, double b);
        int add(SCDoubleMatrix& matrix, char trans, int m, int n, double a, double b);
        double dot(SCDoubleMatrix& matrix, int n);
        int reorder(SCIntMatrix& order, int npts=0);
        double getLocalValue(int i, int j=0) {return _matrix[i+_nlocal*j];};
        int distributeVector();
        int qr(int n=-1);
        SCDoubleMatrix * getQ();
        DoubleInt getMinLoc(int begglo=-1, int endglo=-1);
        DoubleInt getMin(int begglo=-1, int endglo=-1);
        DoubleInt getMaxLoc(int begglo=-1, int endglo=-1);
        DoubleInt getMax(int begglo=-1, int endglo=-1);
        void writeLocal(std::string filename);
        int house(int j);
        int pdcopy(int n, int ix, int jx, int incx,
            SCDoubleMatrix & y, int iy, int jy, int incy);
        int pdlacp2(char uplo, int m, int n, int ia, int ja, SCDoubleMatrix& b, int ib, int jb);
        int initqr();
        void initRandom(int seed=1, double lo=-1.0, double hi=1.0);
        int copyRedist( int m, int n, int ia, int ja, SCDoubleMatrix& B, int ib, int jb, int ctxt);
        int initMatrix(double *A);
        void swap(int i, int j);
        double froNorm();
        double amaxElement();
        void normalize(double fac);
        void scalarMultiply( double s);
        int loadMatrix(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A);
        int loadMatrix(Eigen::Ref<Eigen::VectorXd> &b);
        int loadRhs(const Eigen::Ref<const Eigen::VectorXd> &b);
        void startTime(int i) {_wallclock[i] = -SCBaseMatrix::getWallTime();}
        void stopTime(int i) {_wallclock[i] += SCBaseMatrix::getWallTime(); _wallclock_total[i] += _wallclock[i];}
        double getTime(int i) { return _wallclock_total[i]; }
        double getMaxTime(int i);
        int norm2Colunns(SCDoubleMatrix& colnorms);
        void columnScaling(SCDoubleMatrix& colScale);
        void elementWiseInverse();
        void scaleColumnsByL2Norm(SCDoubleMatrix& colScale);
        bool isFeasible();
        void AtA(SCDoubleMatrix& A, int n);
        int solve(SCDoubleMatrix& x, SCDoubleMatrix& b, int n);
        int choldecomp(int n);
        int cholsolve(SCDoubleMatrix& x, SCDoubleMatrix& b, int n);
        double minSingularValue();
        double maxSingularValue();
        double conditionNumber();
        void writeSingularValues(std::string filename = std::string("svd.txt"));

    private:
        double * _matrix;       // Local _mlocal X _nlocal matrix
        SCDoubleMatrix * _tau;  // Needed only if this matrix holds a QR decomposition compactly.
        double * _sing;         // Holds the singular values if computed
        bool _isQR;             // Flag to say if it is an QR decomposition

        double _wallclock[SCDBL_N_TIMES];
        double _wallclock_total[SCDBL_N_TIMES];

        int init();
        int readMatrix(std::string filename);
        double Norm(char normDesignator);
        int singularValues();
};

#endif // SCDOUBLEMATRIX_H_
