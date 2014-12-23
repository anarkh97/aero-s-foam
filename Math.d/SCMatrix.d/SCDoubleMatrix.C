#ifdef NNLS_DEV
#include "SCDoubleMatrix.h"
#else
#include "Math.d/SCMatrix.d/SCDoubleMatrix.h"
#endif

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <Eigen/Core>
#include <ctime>


SCDoubleMatrix::SCDoubleMatrix(std::string filename, int context, int mb, int nb) :
    SCBaseMatrix(filename, context, mb, nb) {
    init();
    readMatrix(filename);
}


SCDoubleMatrix::SCDoubleMatrix(int context, int m, int n, int mb, int nb) :
    SCBaseMatrix(context, m, n, mb, nb) {
    init();
}


SCDoubleMatrix::SCDoubleMatrix(const SCDoubleMatrix& matrix) : 
    SCBaseMatrix(matrix._context, matrix._m, matrix._n, matrix._mb, matrix._nb) {
    SCBaseMatrix::init();
    SCDoubleMatrix::init();
    for (int i=0; i<_sizelocal; i++) {
        _matrix[i] = matrix._matrix[i];
    }
}


int
SCDoubleMatrix::readMatrix(std::string filename) {
    int m, n, count;
    FILE * fptr = fopen(filename.c_str(), "rb");
    count = fread(&m, sizeof(int), 1, fptr);
    count = fread(&n, sizeof(int), 1, fptr);
    if (m != _m || n != _n) {
        std::cout << "Problem in SCDoubleMatrix::readMatrix." << std::endl;
        exit(1);
    }
    double * col = new double[_m];
    count = 0;
    for (int j=1; j<=_n; j++) {
        count += fread(col, sizeof(double), _m, fptr);
        setMatrixColumn(j, col);
    }
    fclose(fptr);
    delete col;
    return count;
}


int
SCDoubleMatrix::init() {
    _matrix = new double[_sizelocal];
    _isQR = false;
    _tau = NULL;
    for (int i=0; i < SCDBL_N_TIMES; i++) {
        _wallclock[i] = 0.0;
        _wallclock_total[i] = 0.0;
    }   
    return 0;
}

int
SCDoubleMatrix::initqr() {
    int d = std::min(_m, _n);
    _tau = new SCDoubleMatrix(_context,  1,   _m, _mb, _nb);
    _isQR = true;
    return 0;
}


SCDoubleMatrix::~SCDoubleMatrix() {
    delete[] _matrix;
    if (_tau != NULL) {
        delete _tau;
    }
}


void
SCDoubleMatrix::write(std::string filename, bool compact) {
    char scope = 'A';
    char blank = ' ';
    double alpha;
    FILE *f;
    if (_mypid == 0) {
        f = fopen(filename.c_str(), "w");
    }
    for (int i=1; i<=_m; i++) {
        for (int j=1; j<=_n; j++) {
            _FORTRAN(pdelget)(&scope, &blank, &alpha, _matrix, &i, &j, _desc);
            if (_mypid == 0) {
                if (compact) {
                    fprintf(f, "%12.4e ", alpha);
                } else {
                    fprintf(f, "%24.16e ", alpha);
                }
                if (_m == 1 || _n == 1) {
                    fprintf(f, "\n");
                }
            }
        }
        if (_mypid == 0) {
            if ( !(_m == 1 || _n == 1) ) {
                fprintf(f, "\n");
            }
        }
    }
    if (_mypid == 0) {
        fprintf(f, "\n");
        fclose(f);
    }
    if (_isQR) {
        std::string taufile = "tau_" + filename;
        _tau->write(taufile);
    }
}


void
SCDoubleMatrix::writeLocal(std::string filename) {
    float alpha;
    FILE *f;
    f = fopen(filename.c_str(), "w");
    int m = std::max(1, _mlocal);
    int n = std::max(1, _nlocal);
    int ig, jg;
    int zero = 0;
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            ig = _FORTRAN(indxl2g)(&i, &_mb, &_myrow, &zero, &_mprow);
            jg = _FORTRAN(indxl2g)(&j, &_nb, &_mycol, &zero, &_npcol);
            alpha = _matrix[i + j%m];
            if (_mypid == 0) {
                fprintf(f, "%12.4e (%d,%d)\n", alpha, ig, jg);
            }
        }
    }
    fclose(f);
}


// j is an element of [1,_n]; Starts at 1 for Fortran
int
SCDoubleMatrix::setMatrixColumn(int j, double *col) {
    if (j > 0 && j <= _n) {
        for (int i=1; i<=_m; i++) {
            _FORTRAN(pdelset)(_matrix, &i, &j, _desc, &(col[i-1]));
        }   
    } else {
        std::cerr << "Problem in SCDoubleMatrix::setMatrixColumn. j = " << j << " is not a valid column index." << std::endl;
    }
    return 0;
}


int
SCDoubleMatrix::setMatrixRow(int i, double *row) {
    if (i > 0 && i <= _m) {
        for (int j=1; j<=_n; j++) {
            _FORTRAN(pdelset)(_matrix, &i, &j, _desc, &(row[j-1]));
        }   
    } else {
        std::cerr << "Problem in SCDoubleMatrix::setMatrixRow. i = " << i << " is not a valid row index." << std::endl;
    }
    return 0;
}



// Norm returned the process row or column 2 norm associated with the vector
double
SCDoubleMatrix::norm2(int n) {
    double norm2 = 0.0;
    if (_m != 1 || _n != 1) { // mx1 or 1xn only
        int one = 1;
        int m = n;
        if (m == 0) {
            if (_m == 1 ) {
                m = _n;
            } else {
                m = _m;
            }
         }
        _FORTRAN(pdnrm2)(&m, &norm2, _matrix, &one, &one, _desc, &one);
        SCBaseMatrix::distributeVector(&norm2, 1);
    }
    return norm2;
}


int
SCDoubleMatrix::norm2Colunns(SCDoubleMatrix& colnorms) {
    if (colnorms._n != _n) {
        return 1;
    }
    int zero=0, one=1, jloc, dummy, p;
    double norm;
    for (int j=1; j<=_n; j++) {
        _FORTRAN(pdnrm2)(&_m, &norm, _matrix, &one, &j, _desc, &one);
        p = _FORTRAN(indxg2p)(&j, &(colnorms._nb), &dummy, &zero, &(colnorms._npcol));
        if (p == _mycol) {
            jloc = _FORTRAN(indxg2l)(&j, &(colnorms._nb), &dummy, &dummy, &(colnorms._npcol));
            colnorms._matrix[jloc-1] = norm; // jloc is a Fortran index
        }
    }
    return 0;
}


// Note, only for vectors. Not matrices!
int
SCDoubleMatrix::distributeVector() {
    int zero = 0;
    int one = 1;
    if (_m == 1) {
        char scope = 'C';
        if (_myrow == 0) {
            _FORTRAN(dgebs2d)(&_context, &scope, &_top, &one, &_nlocal, _matrix, &one);
        } else {
            _FORTRAN(dgebr2d)(&_context, &scope, &_top, &one, &_nlocal, _matrix, &one, &zero, &_mycol);
        }
    } else if (_n == 1) {
        char scope = 'R';
        if (_mycol == 0) {
            _FORTRAN(dgebs2d)(&_context, &scope, &_top, &_mlocal, &one, _matrix, &_mlocal);
        } else {
            _FORTRAN(dgebr2d)(&_context, &scope, &_top, &_mlocal, &one, _matrix, &_mlocal, &_myrow, &zero);
        }
    } else {
        return 1;
    }
    return 0;
}


int
SCDoubleMatrix::multiply(SCDoubleMatrix &x, SCDoubleMatrix &y, char trans, int m, int n, double alpha, double beta) {
    int one = 1;
    _FORTRAN(pdgemv)( &trans,  &m,   &n, &alpha,
            _matrix, &one, &one, _desc, 
            x._matrix, &one, &one, x._desc, &one, &beta,
            y._matrix, &one, &one, y._desc, &one);
    return 0;
}


int
SCDoubleMatrix::multiply(char trans, int m, int n, double alpha, int ia, int ja,
                         SCDoubleMatrix &x, int ix, int jx, int incx, double beta,
                         SCDoubleMatrix &y, int iy, int jy, int incy) {
    _FORTRAN(pdgemv)(&trans,  &m,   &n, &alpha,
            _matrix, &ia, &ja, _desc, 
            x._matrix, &ix, &jx, x._desc, &incx, &beta,
            y._matrix, &iy, &jy, y._desc, &incy);
    return 0;
}


int
SCDoubleMatrix::hadamardProduct(SCDoubleMatrix &x) {
    if (_sizelocal != x._sizelocal) {
        return 1;
    }
    for (int i=0; i<_sizelocal; i++) {
        _matrix[i] *= x._matrix[i];
    }
    return 0;
}


int
SCDoubleMatrix::zero() {
    std::memset(_matrix, 0, _sizelocal*sizeof(double));
    return 0;
}


int
SCDoubleMatrix::zero(int ix, int jx, int n) {
    double dzero = 0.0;
    for (int i=ix; i<=n; i++) {
        _FORTRAN(pdelset)(_matrix, &i, &jx, _desc, &dzero);
    }
    return 0;
}


int
SCDoubleMatrix::set(double val) {
    for (int i=0; i<_sizelocal; i++) {
        _matrix[i] = val;
    }
    return 0;
}


// Old permute
int
SCDoubleMatrix::permuteOld(char direc, SCIntMatrix &ip, int m, int n) {
    int zero = 0, one = 1;

    if (m == 0) m = _m;
    if (n == 0) n = _n;

    char rowcol, pivroc;
    if (ip.getNumberOfRows() == 1) {
        rowcol = 'C';
        pivroc = 'R';
    } else {
        rowcol = 'R';
        pivroc = 'C';
    }
    int *iwork = NULL;  // iwork only needed if ipiv is needs to be transposed.
    _FORTRAN(fpdlapiv)(&direc, &rowcol, &pivroc, &m, &n, _matrix, &one, &one, _desc, ip.getMatrix(), &one, &one, ip.getDesc(), iwork);
    return 0;
}


// New permute
int
SCDoubleMatrix::permute(char direc, char rowcol, SCIntMatrix &ip, int m, int n) {
    int zero = 0, one = 1;

    if (m == 0) m = _m;
    if (n == 0) n = _n;

    int *iwork = NULL;  // iwork only needed if ipiv is needs to be transposed.
    int lwork = ip.getLworkPdlapiv(rowcol);
    if (lwork > 0) {
        iwork = (int *) malloc(lwork * sizeof(int));
    }
    char pivroc = ip.getPivroc();
    _FORTRAN(fpdlapiv)(&direc, &rowcol, &pivroc, &m, &n, _matrix, &one, &one, _desc,
              ip.getMatrix(), &one, &one, ip.getDesc(), iwork);
    if (iwork != NULL) {
        free(iwork);
    }
    return 0;
}


double
SCDoubleMatrix::getElement(int i, int j) {
    double value;
    _FORTRAN(pdelget)(&_scope, &_top, &value, _matrix, &i, &j, _desc);
    return value;
}


void
SCDoubleMatrix::setElement(int i, int j, double value) {
    _FORTRAN(pdelset)(_matrix, &i, &j, _desc, &value);
}


double
SCDoubleMatrix::getElement(int i, int j, int rsrc, int csrc) {
    int desc[DLEN_];
    for (int k=0; k<DLEN_; k++) {
        desc[k] = _desc[k];
    }
    desc[RSRC_] = rsrc;
    desc[CSRC_] = csrc;
    double value;
    _FORTRAN(pdelget)(&_scope, &_top, &value, _matrix, &i, &j, desc);
    return value;
}


void
SCDoubleMatrix::setElementsLocal(const SCDoubleMatrix& matrix) {
    for (int i=0; i<_sizelocal; i++) {
        _matrix[i] = matrix._matrix[i];
    }
}


// Only for vectors
void
SCDoubleMatrix::project() {
    if (_m == 1) {
        if (_myrow == 0) {
            for (int i=0; i<_nlocal; i++) {
                if (_matrix[i] < 0.0) {
                    _matrix[i] = 0.0;
                }
            }
        }
    } else if (_n == 1) {
        if (_mycol == 0) {
            for (int i=0; i<_mlocal; i++) {
                if (_matrix[i] < 0.0) {
                    _matrix[i] = 0.0;
                }
            }
        }
    }
}


int
SCDoubleMatrix::loadIdentityMatrix(double value) {
    // Initialize the matrix to the identity
    int m = std::max(_m, _n);
    SCDoubleMatrix::zero();
    for (int i=1; i<=m; i++) {
        _FORTRAN(pdelset)(_matrix, &i, &i, _desc, &value);
    }

    return 0;
}


int
SCDoubleMatrix::pdcopy(int n, int ix, int jx, int incx,
    SCDoubleMatrix & y, int iy, int jy, int incy) {
    _FORTRAN(pdcopy)(&n, _matrix, &ix, &jx, _desc, &incx,
        y.getMatrix(), &iy, &jy, y.getDesc(), &incy);
    return 0;
}


int
SCDoubleMatrix::pdlacp2(char uplo, int m, int n, int ia, int ja, SCDoubleMatrix& b, int ib, int jb) {
    _FORTRAN(pdlacp2)(&uplo, &m, &n, _matrix, &ia, &ja, _desc, b._matrix, &ib, &jb, b._desc);
    return 0;
}


// Vectors only
int
SCDoubleMatrix::copy(SCDoubleMatrix& matrix, int n) {
    int one = 1;
    if ( (_m == 1 && matrix._m == 1) || (_n == 1 && matrix._n == 1) ) {
        _FORTRAN(pdcopy)(&n,
            _matrix,        &one, &one, _desc,        &one,
            matrix._matrix, &one, &one, matrix._desc, &one);
    } else if ( (_m == 1 && matrix._n == 1) || (_n == 1 && matrix._m == 1) ) {
        char scope;
        double value;
        if (_m == 1) {
            scope = 'A';
            for (int j=1; j <=n; j++) {
                _FORTRAN(pdelget)(&scope, &_top, &value, _matrix, &one, &j, _desc);
                matrix.setElement(j, one, value);
            }
        } else if (_n == 1) {
            scope = 'A';
            for (int i=1; i <=n; i++) {
                _FORTRAN(pdelget)(&scope, &_top, &value, _matrix, &i, &one, _desc);
                matrix.setElement(one, i, value);
            }
        } else {
            std::cerr << "Error in SCDoubleMatrix::copy" << std::endl; 
            exit(1);
        }
    } else {
        std::cerr << "SCDoubleMatrix::copy is for vectors only. Must have _m == 1  or _n == 1" << std::endl; 
        exit(1);
    }
    return 0;
}


int
SCDoubleMatrix::copy(SCDoubleMatrix& dest, int n, SCIntMatrix& order) {
    this->setScope();
    dest.setScope();
    order.setScope();
    if (_m == 1) {
        if (_myrow == 0) {
            double value;
            int k;
            if (n == -1) n = _n;
            for (int j=1; j<=n; j++) {
                k = order.getElement(1,j);
                if (k > 0) {
                    value = this->getElement(1,k);
                    dest.setElement(1,j,value);
                }
            }
        }
    } else if (_n == 1) {
        if (_mycol == 0) {
            double value;
            int k;
            if (n == -1) n = _m;
            for (int i=1; i<=n; i++) {
                k = order.getElement(i,1);
                if (k > 0) {
                    value = this->getElement(k,1);
                    dest.setElement(i,1,value);
                }
            }
        }
    }
    return 0;
}


int
SCDoubleMatrix::add(SCDoubleMatrix& matrix, char trans, int m, int n, double a, double b) {
    int one = 1;
    _FORTRAN(pdgeadd)(&trans, &m, &n, &a,
            _matrix,        &one, &one, _desc, &b,
            matrix._matrix, &one, &one, matrix._desc);
    return 0;
}


int
SCDoubleMatrix::add(SCDoubleMatrix& matrix, char trans, int n, double a, double b) {
    int one = 1;
    if (_m == matrix._m && _n == matrix._n) {
         _FORTRAN(pdgeadd)(&trans, &n, &one, &a,
            _matrix,        &one, &one, _desc, &b,
            matrix._matrix, &one, &one, matrix._desc);
    } else {
        std::cerr << "Problem in SCDoubleMatrix::add" << std::endl;
    }
    return 0;
}


// Only for vectors
double
SCDoubleMatrix::dot(SCDoubleMatrix& matrix, int n) {
    int one = 1;
    double dot;
    _FORTRAN(pddot)(&n, &dot,
           _matrix,        &one, &one, _desc,        &one,
           matrix._matrix, &one, &one, matrix._desc, &one);
    return dot;
}


// Only for vectors.
int
SCDoubleMatrix::reorder(SCIntMatrix& order, int npts) {
    SCDoubleMatrix *mtmp = new SCDoubleMatrix(*this);
    if (_m == 1) {
        int n = _n;
        if (npts > 0) {
            n = npts;
        }
        if (_myrow == 0) {
            double value;
            int k;
            mtmp->setScope();
            this->zero();
            for (int j=1; j<=n; j++) {
                k = order.getElement(1,j);
                value = mtmp->getElement(1,j);
                this->setElement(1,k,value);
            }
        }
    } else if (_n == 1) {
        int m = _m;
        if (npts > 0) {
            m = npts;
        }
        if (_mycol == 0) {
            double value;
            int k;
            mtmp->setScope();
            this->zero();
            for (int i=1; i<=m; i++) {
                k = order.getElement(i,1);
                value = mtmp->getElement(i,1);
                this->setElement(k,1,value);
            }
        }
    }
    delete mtmp;
    return 0;
}


int
SCDoubleMatrix::qr(int n) {
    if (!_isQR) {
        initqr();
    }
    int ncols;
    if (n == -1) {
        ncols = _n;
    } else {
        ncols = n;
    }
    int info;
    int zero=0, one=1;
    int lwork = -1;
    double dwork;
    _FORTRAN(pdgeqrf)(&_m, &ncols, _matrix, &one, &one, _desc, _tau->getMatrix(), &dwork, &lwork, &info);

    lwork = int(dwork);
    double * work = new double[lwork];
    _FORTRAN(pdgeqrf)(&_m, &ncols, _matrix, &one, &one, _desc, _tau->getMatrix(), work, &lwork, &info);

    if (info != 0) {
        std::cout << "qr info = " << info << std::endl;
    }
    _tau->write("tau-qr.txt");
    delete[] work;
    return info;
}


DoubleInt
SCDoubleMatrix::getMaxLoc(int begglo, int endglo) {
    this->startTime(SCDBL_TIME_GETMAXLOC);
    DoubleInt maxval;
    int mypc1, mypc2, np, nb;
    if (_m == 1) {
        np    = _npcol;
        nb    = _nb;
        mypc1 = _mycol;
        mypc2 = _myrow;
    } else if (_n == 1) {
        np    = _mprow;
        nb    = _mb;
        mypc1 = _myrow;
        mypc2 = _mycol;
    } else {
        std::cerr << "SCDoubleMatrix::getMaxLoc is for vectors only. Requires _m == 1 or _n == 1" << std::endl;
        MPI_Finalize();
        exit (-1);
    }
    maxval.x = NEGATIVE_INF;
    maxval.i = -1;
    if (mypc2 == 0) {
        int zero = 0;
        if (begglo < 0) { // Check all
            int imaxval = 0;
            maxval.x  = _matrix[imaxval];
            for (int i=1; i<_sizelocal; i++) {
                if (_matrix[i] > maxval.x) {
                    maxval.x  = _matrix[i];
                    imaxval = i;
                }
            }
            imaxval++; // Fortran index
            maxval.i = _FORTRAN(indxl2g)(&imaxval, &nb, &mypc1, &zero, &np);
        } else {
            int j, p, iloc, dummy=0;
            for (int i=begglo; i<=endglo; i++) {
                p = _FORTRAN(indxg2p)(&i, &nb, &dummy, &zero, &np);
                if (p == mypc1) {
                    iloc = _FORTRAN(indxg2l)(&i, &nb, &dummy, &dummy, &np);
                    j = iloc - 1; // Fortran to C++ index
                    if (_matrix[j] > maxval.x) {
                        maxval.x = _matrix[j];
                        maxval.i = i;
                    }
                }
            }
        }
    }
    this->stopTime(SCDBL_TIME_GETMAXLOC);
    return maxval;
}


DoubleInt
SCDoubleMatrix::getMax(int begglo, int endglo) {
    this->startTime(SCDBL_TIME_GETMAX);
    if (! _row_col_comm_set) {
        SCBaseMatrix::setRowColComms();
    }
    MPI_Comm comm, commd;
    int mypc2;
    if (_m == 1) {
        mypc2 = _myrow;
        comm = _row_comm;
        commd = _col_comm;
    } else if (_n == 1) {
        mypc2 = _mycol;
        comm = _col_comm;
        commd = _row_comm;
    } else {
        std::cerr << "SCDoubleMatrix::getMax is for vectors only. Requires _m == 1 or _n == 1" << std::endl;
        MPI_Finalize();
        exit (-1);
    }
    DoubleInt sendbuf, maxval;
    if (mypc2 == 0) {
        sendbuf = getMaxLoc(begglo, endglo);
        MPI_Allreduce(&sendbuf, &maxval, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
    }
    MPI_Barrier(commd);
    MPI_Bcast(&maxval, 1, MPI_DOUBLE_INT, 0, commd);
    this->stopTime(SCDBL_TIME_GETMAX);
    return maxval;
}



DoubleInt
SCDoubleMatrix::getMinLoc(int begglo, int endglo) {
    this->startTime(SCDBL_TIME_GETMINLOC);
    DoubleInt minval;
    int mypc1, mypc2, np, nb;
    if (_m == 1) {
        np    = _npcol;
        nb    = _nb;
        mypc1 = _mycol;
        mypc2 = _myrow;
    } else if (_n == 1) {
        np    = _mprow;
        nb    = _mb;
        mypc1 = _myrow;
        mypc2 = _mycol;
    } else {
        std::cerr << "SCDoubleMatrix::getMaxLoc is for vectors only. Requires _m == 1 or _n == 1" << std::endl;
        MPI_Finalize();
        exit (-1);
    }
    minval.x = POSITIVE_INF;
    minval.i = -1;
    if (mypc2 == 0) {
        int zero = 0;
        if (begglo < 0) { // Check all
            int iminval = 0;
            minval.x  = _matrix[iminval];
            for (int i=1; i<_sizelocal; i++) {
                if (_matrix[i] < minval.x) {
                    minval.x  = _matrix[i];
                    iminval = i;
                }
            }
            iminval++; // Fortran index
            minval.i = _FORTRAN(indxl2g)(&iminval, &nb, &mypc1, &zero, &np);
        } else {
            int j, p, iloc, dummy=0;
            for (int i=begglo; i<=endglo; i++) {
                p = _FORTRAN(indxg2p)(&i, &nb, &dummy, &zero, &np);
                if (p == mypc1) {
                    iloc = _FORTRAN(indxg2l)(&i, &nb, &dummy, &dummy, &np);
                    j = iloc - 1; // Fortran to C++ index
                    if (_matrix[j] < minval.x) {
                        minval.x = _matrix[j];
                        minval.i = i;
                    }
                }
            }
        }
    }
    this->stopTime(SCDBL_TIME_GETMINLOC);
    return minval;
}


DoubleInt
SCDoubleMatrix::getMin(int begglo, int endglo) {
    this->startTime(SCDBL_TIME_GETMIN);
    if (! _row_col_comm_set) {
        SCBaseMatrix::setRowColComms();
    }
    int mypc2;
    MPI_Comm comm, commd;
    if (_m == 1) {
        mypc2 = _myrow;
        comm = _row_comm;
        commd = _col_comm;
    } else if (_n == 1) {
        mypc2 = _mycol;
        comm = _col_comm;
        commd = _row_comm;
    } else {
        std::cerr << "SCDoubleMatrix::getMin is for vectors only. Requires _m == 1 or _n == 1" << std::endl;
        MPI_Finalize();
        exit (-1);
    }
    DoubleInt sendbuf, minval;
    if (mypc2 == 0) {
        sendbuf = getMinLoc(begglo, endglo);
        MPI_Allreduce(&sendbuf, &minval, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
    }
    MPI_Barrier(commd);
    MPI_Bcast(&minval, 1, MPI_DOUBLE_INT, 0, commd);
    this->stopTime(SCDBL_TIME_GETMIN);
    return minval;
}


int
SCDoubleMatrix::house(int j) {
    int retval = -1;
    int m = _m - j + 1;
    int zero=0, one=1, dummy;
    double mu;
    _FORTRAN(pdnrm2)(&m, &mu, _matrix, &j, &j, _desc, &one);
    int ip = _FORTRAN(indxg2p)(&j, &_nb, &dummy, &zero, &_npcol);
    char scope = 'R';
    if (_mycol == ip) {
        _FORTRAN(dgebs2d)(&_context, &scope, &_top, &one, &one, &mu, &one);
    } else {
        _FORTRAN(dgebr2d)(&_context, &scope, &_top, &one, &one, &mu, &one, &_myrow, &ip);
    }
    if (mu != 0.0) {
        double beta, v1;
        scope = 'C';
        _FORTRAN(pdelget)(&scope, &_top, &v1, _matrix, &j, &j, _desc);
        int jp = _FORTRAN(indxg2p)(&j, &_nb, &_mycol, &zero, &_npcol);
        scope = 'R';
        if (_mycol == jp) {
            if (v1 > 0.0) {
                beta = v1 + mu;
            } else {
                beta = v1 - mu;
            }
            _FORTRAN(dgebs2d)(&_context, &scope, &_top, &one, &one, &beta, &one);
        } else {
            _FORTRAN(dgebr2d)(&_context, &scope, &_top, &one, &one, &beta, &one, &_myrow, &jp);
        }
        double v;
        double a = 1.0/beta, b = 0.0;
        char trans = 'N';
        int m = _m-j;
        int ia = j+1;
        _FORTRAN(pdgeadd)(&trans, &m, &one, &a,
            _matrix, &ia, &j, _desc, &b,
            _matrix, &ia, &j, _desc);

        double sigma = fabs(beta/mu);
        if (v1 > 0.0) {
            mu = -mu;
        }
        _FORTRAN(pdelset)(_tau->getMatrix(), &one, &j, _tau->getDesc(), &sigma);
        _FORTRAN(pdelset)(_matrix, &j, &j, _desc, &mu);
        _tau->distributeVector();  // Can do better
        retval = 0;
    }
    return retval;
}


void
SCDoubleMatrix::initRandom(int m, int n) {
    for (int i=0; i<_sizelocal; i++) {
        _matrix[i] = (rand()%1000000)/1000000.0;
    }
    std::cout << "In initRandom() n = " << n << ", m = " << m << std::endl;
    if (n > m) {
        for (int i=0; i<_sizelocal; i++) {
            if (rand()%n > m) {
                _matrix[i] = 0.0;
            }
        }
    }
}


int
SCDoubleMatrix::copyRedist( int m, int n, int ia, int ja, SCDoubleMatrix& B, int ib, int jb, int ctxt) {
    _FORTRAN(pdgemr2d)(&m, &n, _matrix, &ia, &ja, _desc, B.getMatrix(), &ib, &jb, B.getDesc(), &ctxt); 
    return 0;
}


int
SCDoubleMatrix::initMatrix(double *matrix) {
    int iglo, jglo, kglo, k;
    int ii, jj;
    int zero=0;
    int m = std::max(1, _mlocal);
    int n = std::max(1, _nlocal);
    for (int j=0; j<n; j++) {
        for (int i=0; i<m; i++) {
            ii = i+1; // Fortran index
            jj = j+1; // Fortran index
            iglo = _FORTRAN(indxl2g)(&ii, &_mb, &_myrow, &zero, &_mprow)-1; // Back to C index
            jglo = _FORTRAN(indxl2g)(&jj, &_nb, &_mycol, &zero, &_npcol)-1; // Back to C index
            k    = j*m + i;
            kglo = jglo*_m + iglo;
            _matrix[k] = matrix[kglo];
        }
    }
    return 0;
}


void
SCDoubleMatrix::swap(int i, int j) {
    double tmp;
    if (_m == 1) {
        if (_myrow == 0) {
            tmp = this->getElement(1,i);
            this->setElement(1, i, this->getElement(1,j));
            this->setElement(1, j, tmp);
        }
    } else if (_n == 1) {
        if (_mycol == 0) {
            tmp = this->getElement(i,1);
            this->setElement(i, 1, this->getElement(j,1));
            this->setElement(j, 1, tmp);
        }
    } else {
        std::cout << "SCIntMatrix::swap(i,j) only handles vectors" << std::endl;
    }
}


double
SCDoubleMatrix::froNorm() {
    return Norm('F');
}


double
SCDoubleMatrix::amaxElement() {
    return Norm('M');
}


double
SCDoubleMatrix::Norm(char normDesignator) {
    int one = 1;
    double dnorm;
    double *work = NULL;
    if (normDesignator == '1' || normDesignator == 'O' || normDesignator == 'o') {
        work = new double[_nlocal];
    } else if (normDesignator == 'I' || normDesignator == 'i') {
        work = new double[_mlocal];
    }
    dnorm = _FORTRAN(pdlange)(&normDesignator, &_m, &_n, _matrix, &one, &one, _desc, work);
    if (work != NULL) {
        delete[] work;
    }
    return dnorm;
}


void
SCDoubleMatrix::normalize(double fac) {
    if (fac == 0.0) return;
    double dfac = 1.0/fac;
    for (int i=0; i<_sizelocal; i++) {
        _matrix[i] *= dfac;
    }
}


void
SCDoubleMatrix::scalarMultiply( double s) {
    for (int i=0; i<_sizelocal; i++) {
        _matrix[i] *= s;
    }
}


int
SCDoubleMatrix::loadMatrix(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A) {
    int * sizes = new int[_nprocs+1];
    int ldim = A.size() + 1;
    int * sizeslocal = new int[ldim];
    int ncol = 0;
    sizeslocal[0] = 0;
    for (int i=0; i<A.size(); i++) {
        ncol += A[i].cols();
        sizeslocal[i+1] = ncol;
    }
    // std::cout << "Number of columns = " << A[0].cols() << " on processor " << _mypid << std::endl;
    MPI_Allgather(&ncol, 1, MPI_INT, &(sizes[1]), 1, MPI_INT, MPI_COMM_WORLD);

    sizes[0] = 0;
    for (int i=1; i<_nprocs+1; i++) {
        sizes[i] += sizes[i-1];
    }

    double * col = new double[_m];
    int pc, pr, ig, jl, jg, dummy, ioff, zero=0;
    int mlocal, isd, jsd;
    MPI_Status status;
    MPI_Request request;
    for (int p=0; p<_nprocs; p++) {
        if (p == _mypid) {
            isd = 0; // Local sub domain index
            for (int j=sizes[p]; j<sizes[p+1]; j++) {
                jg = j+1;   // Global Eigen column fortran index - actual Matrix column
                jsd = j-sizes[p]-sizeslocal[isd]; // Local Eigen column C index
                if (jsd >= A[isd].cols()) {
                    isd++;
                    jsd = j-sizes[p]-sizeslocal[isd];
                }
                pc  = _FORTRAN(indxg2p)(&jg, &_nb, &dummy, &zero, &_npcol); // Processor col coord of col jg
                for (int k=0; k<_mprow; k++) {
                    pr = _FORTRAN(blacs_pnum)(&_context, &k, &pc);
                    if (pr != p) {
                        mlocal = _FORTRAN(numroc)(&_m, &_mb, &k, &zero, &_mprow);
                        for (int i=1; i<=mlocal; i++) {
                            ig = _FORTRAN(indxl2g)(&i, &_mb, &k, &zero, &_mprow);
                            col[i-1] = A[isd].col(jsd)[ig-1];
                        }   
                        MPI_Send(col, mlocal, MPI_DOUBLE, pr, k, MPI_COMM_WORLD);
                        //MPI_Isend(col, mlocal, MPI_DOUBLE, pr, k, MPI_COMM_WORLD, &request);
                    } else {
                        jl = _FORTRAN(indxg2l)(&jg, &_nb, &zero, &zero, &_npcol);
                        ioff = (jl-1)*_mlocal;
                        for (int i=1; i<=_mlocal; i++) {
                            ig = _FORTRAN(indxl2g)(&i, &_mb, &_myrow, &zero, &_mprow);
                            _matrix[ioff+i-1] = A[isd].col(jsd)[ig-1];
                        }   
                    }
                }
            }
        } else {
            for (int j=sizes[p]; j<sizes[p+1]; j++) {
                jg = j+1;   // Global Eigen column fortran index - actual Matrix column
                pc  = _FORTRAN(indxg2p)(&jg, &_nb, &dummy, &zero, &_npcol); // Processor col-coord of col jg
                if  (_mycol == pc) {
                    MPI_Recv(col, _mlocal, MPI_DOUBLE, p, _myrow, MPI_COMM_WORLD, &status);
                    jl = _FORTRAN(indxg2l)(&jg, &_nb, &zero, &zero, &_npcol);
                    ioff = (jl-1)*_mlocal;
                    for (int i=0; i<_mlocal; i++) {
                        _matrix[ioff+i] = col[i];
                    }   
                }
            }
        }
    }
    delete[] col;
    delete[] sizes;
    delete[] sizeslocal;
    return 0;
}


/* Loads the transpose of an Eigen matrix
int
SCDoubleMatrix::loadTransposeMatrix(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A) {
    int * sizes = new int[_nprocs+1];
    int ldim = A.size() + 1;
    int * sizeslocal = new int[ldim];
    int ncol = 0;
    int one = 1;
    sizeslocal[0] = 0;
    for (int i=0; i<A.size(); i++) {
        ncol += A[i].cols();
        sizeslocal[i+1] = ncol;
    }
    // std::cout << "Number of columns = " << A[0].cols() << " on processor " << _mypid << std::endl;
    MPI_Allgather(&ncol, 1, MPI_INT, &(sizes[1]), 1, MPI_INT, MPI_COMM_WORLD);

    sizes[0] = 0;
    for (int i=1; i<_nprocs+1; i++) {
        sizes[i] += sizes[i-1];
    }

    double * row = new double[_n]; // can do much much better
    int pc, pr, ig, il, jg, dummy, ioff, zero=0;
    int nlocal, isd, jsd;
    MPI_Status status;
    MPI_Request request;
    for (int p=0; p<_nprocs; p++) {
        if (p == _mypid) {
            isd = 0; // Local sub domain index
            for (int j=sizes[p]; j<sizes[p+1]; j++) {
                jg = j+1;   // Global Eigen column fortran index
                jsd = j-sizes[p]-sizeslocal[isd]; // Local Eigen column C index
                if (jsd >= A[isd].cols()) {
                    isd++;
                    jsd = j-sizes[p]-sizeslocal[isd];
                }
                pr = _FORTRAN(indxg2p)(&jg, &_mb, &dummy, &zero, &_mprow); // Processor row coord of Eigen col jg
                for (int k=0; k<_npcol; k++) {
                    pc = _FORTRAN(blacs_pnum)(&_context, &pr, &k);
                    if (pc != p) {
                        nlocal = _FORTRAN(numroc)(&_n, &_nb, &k, &zero, &_npcol);
                        for (int j=1; j<=nlocal; j++) {
                            ig = _FORTRAN(indxl2g)(&j, &_nb, &k, &zero, &_npcol) - 1;
                            row[j-1] = A[isd].col(jsd)[ig];
                        }   
                        MPI_Send(row, nlocal, MPI_DOUBLE, pr, k, MPI_COMM_WORLD);
                        //MPI_Isend(row, nlocal, MPI_DOUBLE, pr, k, MPI_COMM_WORLD, &request);
                    } else {
                        il = _FORTRAN(indxg2l)(&jg, &_mb, &zero, &zero, &_mprow) - 1;
                        for (int i=1; i<=_nlocal; i++) {
                            ig = _FORTRAN(indxl2g)(&i, &_nb, &_mycol, &zero, &_npcol) -1;
                            _matrix[il+(i-1)*_mlocal] = A[isd].col(jsd)[ig];
                        }   
                    }
                }
            }
        } else {
            for (int j=sizes[p]; j<sizes[p+1]; j++) {
                jg = j+1;   // Global Eigen column fortran index - actual Matrix row index
                pr = _FORTRAN(indxg2p)(&jg, &_mb, &dummy, &zero, &_mprow); // Processor row coord of Eigen col jg
                if  (_myrow == pr) {
                    MPI_Recv(row, _nlocal, MPI_DOUBLE, p, _mycol, MPI_COMM_WORLD, &status);
                    il = _FORTRAN(indxg2l)(&jg, &_mb, &zero, &zero, &_mprow) - 1;
                    for (int i=0; i<_nlocal; i++) {
                        _matrix[il+i*_mlocal] = row[i];
                    }   
                }
            }
        }
    }
    delete[] row;
    delete[] sizes;
    delete[] sizeslocal;
    return 0;
}
*/


int
SCDoubleMatrix::loadMatrix(Eigen::Ref<Eigen::VectorXd> &b) {
    Eigen::Map<Eigen::MatrixXd> em = Eigen::Map<Eigen::MatrixXd>(b.data(), b.rows(), 1);
    std::vector< Eigen::Map<Eigen::MatrixXd> > vec;
    vec.push_back(em);
    this->loadMatrix( vec );
    return 0;
}


int
SCDoubleMatrix::loadRhs(const Eigen::Ref<const Eigen::VectorXd> &b) {
    int zero=0, one=1, ig;
    if (_mycol == 0) {
        for (int i=1; i<= _sizelocal; i++) {
            ig = _FORTRAN(indxl2g)(&i, &_mb, &_myrow, &zero, &_mprow);
            _matrix[i-1] = b(ig-1);
        }
    }
    return 0;
}


// Max returned on processor zero
double
SCDoubleMatrix::getMaxTime(int i) {
    double wall_time = this->getTime(i);
    double max_time;
    MPI_Reduce(&wall_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return max_time;
}


void
SCDoubleMatrix::elementWiseInverse() {
    for (int i=0; i<_sizelocal; i++) {
        if (_matrix[i] != 0.0) {
            _matrix[i] = 1.0/_matrix[i];
        }
    }
}


char
SCDoubleMatrix::columnScaling(SCDoubleMatrix& colScale) {
    DoubleInt cmax = colScale.getMax();
    DoubleInt cmin = colScale.getMin();
    if (cmin.x <= 0.0) {
        if (_mypid == 0) {
            std::cout << "Bad scaling matrix. No scaling will be done." << std::endl;
        }
        return 'N';
    }
    char equed;
    SCDoubleMatrix *rowScale = new SCDoubleMatrix(_context, _m, 1, _mb, _nb);
    rowScale->set(1.0);
    rowScale->distributeVector();
    //colScale.distributeVector();  // Assume this is already distributed
    double amax = this->amaxElement();
    double rowcnd = 1.0;
    double colcnd = cmin.x / cmax.x;
    //std::cout << "_mypid = " << _mypid << ", amax = " << amax << std::endl;
    //std::cout << "_mypid = " << _mypid << ", rowcnd = " << rowcnd << std::endl;
    //std::cout << "_mypid = " << _mypid << ", colcnd = " << colcnd << std::endl;
    int one = 1;
    _FORTRAN(pdlaqge)(&_m, &_n, _matrix, &one, &one, _desc, rowScale->_matrix,
                     colScale._matrix, &rowcnd, &colcnd, &amax, &equed);
    delete rowScale;
    return equed;
}
