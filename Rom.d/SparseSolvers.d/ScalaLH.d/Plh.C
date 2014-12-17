#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <sys/time.h>
#include <string>
#include <vector>
#include <iterator>

#include "Plh.h"


Plh::Plh() {
    _m = 0;
    _n = 0;
    initWithSize();
}


Plh::Plh(int m, int n) : _m(m), _n(n) {
    initWithSize();
}


Plh::Plh(const std::vector< Eigen::Map<Eigen::MatrixXd> >& A) {
    initWithEigen(A);
}


void
Plh::initDefaults() {
    _FORTRAN(blacs_pinfo)(&_mypid, &_nprocs);

    // Solver
    _rtol    = 0.0;
    _maxNP   = 0;

    // A Context
    defaultProcGrid();
    _mb    = MBA_DEFAULT;
    _nb    = NBA_DEFAULT;

    // QR Context
    _mprowQR = _nprocs;
    _npcolQR = 1;
    _mbq     = MBQ_DEFAULT;
    _nbq     = NBQ_DEFAULT;

    _residualIncr = -1;
    _residualFilePtr = NULL;

    // Downdate masking. A test...
    _cnmask = false;
    _ddmask = false;
    _wmask = NULL;
    _colnorms = NULL;
}


void
Plh::initWithSize() {
    initDefaults();
    _initializedWithEigen = false;
    _matrixInitialized = false;
    _max_iter = _n*MAXITE_DEFAULT;
}


// Determines m & n and:
//     _nEigenSubDomains
//     _eigenSubDomainSize
//     _eigenColsPerProc
//     _eigenSubdomainStart
void
Plh::initWithEigen(const std::vector< Eigen::Map<Eigen::MatrixXd> >& A) {
    initDefaults();
    if (_mypid == 0) {
        std::cout << "Initializing with a distributed Eigen matrix." << std::endl;
    }

    if (A.size() == 0) {
        std::cout << "No Eigen matrix or processor " << _mypid << ". Exiting..." << std::endl;
        MPI_Finalize();
        exit(1);
    }

    int ncols = 0;
    _nEigenSubDomains = A.size();
    _eigenSubDomainSize = new int[_nEigenSubDomains];
    for (int i=0; i<_nEigenSubDomains; i++) {
        _eigenSubDomainSize[i] = A[i].cols();
        ncols += A[i].cols();
    }

    // Eigen cols per processor
    _eigenColsPerProc = new int[_nprocs];
    MPI_Allgather(&ncols, 1, MPI_INT, _eigenColsPerProc, 1, MPI_INT, MPI_COMM_WORLD);

    // Starting point for Eigen global column indicies 
    _eigenSubdomainStart = new int[_nEigenSubDomains];
    int nstart = 0;
    for (int p=1; p<=_mypid; p++) {
        nstart += _eigenColsPerProc[p-1];
    }
    _eigenSubdomainStart[0] = nstart;
    for (int isd=1; isd<_nEigenSubDomains; isd++) {
        _eigenSubdomainStart[isd] = _eigenSubdomainStart[isd-1] + _eigenSubDomainSize[isd-1];
    }

    // Set matrix size
    _m = A[0].rows();
    _n = 0;
    for (int p=0; p<_nprocs; p++) {
        _n += _eigenColsPerProc[p];
    }

    defaultProcGrid();
    _initializedWithEigen = true;
    _matrixInitialized = false;
    _max_iter = _n*MAXITE_DEFAULT;
}


void
Plh::setMatrixSize(int m, int n) {
    _m = m;
    _n = n;
    _max_iter = _n*MAXITE_DEFAULT;
    defaultProcGrid();
}


void
Plh::setResidualFileName(std::string filename, int incr) {
    _residualIncr = incr;
    _residualFileName = filename;
    _residualFilePtr = NULL;
    if (_mypid == 0) {
        _residualFilePtr = fopen(_residualFileName.c_str(), "w");
    }
}


void
Plh::singleContext() {
    // QR Context
    _mprowQR = _mprow;
    _npcolQR = _npcol;
}


void
Plh::init(const std::vector< Eigen::Map<Eigen::MatrixXd> >& eigenMatrix, const Eigen::Ref<const Eigen::VectorXd>& eigenRhs) {
    init();
    loadMatrix(eigenMatrix);
    loadRhs(eigenRhs);
}


Plh::~Plh() {
    delete _A;
    delete _x;
    delete _b;
    delete _w;
    delete _wmask;
    delete _workm;

    delete _Q;
    delete _Qtb;
    delete _rhs;
    delete _zQR;
    delete _xQR;
    delete _wQR;
    delete _QtoA;

    if (_work_qr) delete[] _work_qr;

    if (_mypid == 0) {
        if (_residualFilePtr != NULL) {
            fclose(_residualFilePtr);
        }
    }

    if (_initializedWithEigen) {
        delete[] _eigenSubDomainSize;
        delete[] _eigenColsPerProc;
        delete[] _eigenSubdomainStart;
    }
}


void
Plh::init() {
    // std::cout << "Begin init()" << std::endl;
    _FORTRAN(blacs_pinfo)(&_mypid, &_nprocs);

    int ic = -1;
    int zero = 0, one = 1;
    char order = 'R';

    _FORTRAN(blacs_get)(&ic, &zero, &_context);
    _FORTRAN(blacs_gridinit)(&_context, &order,  &_mprow, &_npcol );
    _FORTRAN(blacs_gridinfo)(&_context, &_mprow, &_npcol, &_myrow, &_mycol);

    // Comm group for rows and columns
    MPI_Comm_split(MPI_COMM_WORLD, _myrow, _mycol, &_row_comm); 
    MPI_Comm_split(MPI_COMM_WORLD, _mycol, _myrow, &_col_comm); 

    // Get context for QR decomposition. 
    if (_mprow == _mprowQR && _npcol == _npcolQR) {
        _contextQR = _context;
        _myrowQR = _myrow;
        _mycolQR = _mycol;
        _row_commQR = _row_comm;
        _col_commQR = _col_comm;
    } else {
        _FORTRAN(blacs_get)(&ic, &zero, &_contextQR);
        _FORTRAN(blacs_gridinit)( &_contextQR, &order, &_nprocs, &one);
        _FORTRAN(blacs_gridinfo)(&_contextQR, &_mprowQR, &_npcolQR, &_myrowQR, &_mycolQR);
        MPI_Comm_split(MPI_COMM_WORLD, _myrowQR, _mycolQR, &_row_commQR); 
        MPI_Comm_split(MPI_COMM_WORLD, _mycolQR, _myrowQR, &_col_commQR); 
    }

    int dmax = std::max(_m, _n);
    int dmin = std::min(_m, _n);
    if (_maxNP > 0) {
        dmin = std::min(dmin, _maxNP);
    } else {
        _maxNP = dmin;
    }

    // Initialize vectors and matrices
    // Context c_A
    _A     = new SCDoubleMatrix(_context,   _m,   _n, _mb, _nb);
    _x     = new SCDoubleMatrix(_context,    1,   _n, _mb, _nb);
    _b     = new SCDoubleMatrix(_context,   _m,    1, _mb, _nb);
    _w     = new SCDoubleMatrix(_context,    1,   _n, _mb, _nb);
    _workm = new SCDoubleMatrix(_context,   _m,    1, _mb, _nb);
    _Atb   = new SCDoubleMatrix(_context,    1,   _n, _mb, _nb);
    _QtoA  = new SCIntMatrix(   _context,    1,   _n, _mb, _nb);

    // Context c_Q
    _Q       = new SCDoubleMatrix(_contextQR,   _m, dmin, _mbq, _nbq);
    _Qtb     = new SCDoubleMatrix(_contextQR,   _m,    1, _mbq, _nbq);
    _rhs     = new SCDoubleMatrix(_contextQR, dmax,    1, _mbq, _nbq);
    _zQR     = new SCDoubleMatrix(_contextQR,    1,   _n, _mbq, _nbq);
    _xQR     = new SCDoubleMatrix(_contextQR,    1,   _n, _mbq, _nbq);
    _wQR     = new SCDoubleMatrix(_contextQR,    1,   _n, _mbq, _nbq);
    _bQR     = new SCDoubleMatrix(_contextQR,   _m,    1, _mbq, _nbq);
    _rQR     = new SCDoubleMatrix(_contextQR,   _m,    1, _mbq, _nbq);
    _workmQR = new SCDoubleMatrix(_contextQR,   _m,    1, _mbq, _nbq); // Change _mb to _mbq and _nb to _nbq


    _Q->initqr();
    _work_qr = NULL;

    for (int i=0; i < N_TIMES; i++) {
        _wallclock[i] = 0;
        _wallclock_total[i] = 0;
    }

    //_w->testCommunicators();

    _matrixInitialized = true;
}


void
Plh::summary() {
    if (_mypid == 0) {
        std::cout << std::endl;
        std::cout << "################### Solver Summary ########################" << std::endl;
        std::cout << " m               = "               << _m << std::endl;
        std::cout << " n               = "               << _n << std::endl;
        std::cout << "A context:"                        << std::endl;
        std::cout << "    mprow        = " << _mprow     << std::endl;
        std::cout << "    npcol        = " << _npcol     << std::endl;
        std::cout << "    mba          = " << _mb        << std::endl;
        std::cout << "    nba          = " << _nb        << std::endl;
        std::cout << "QR context:"                       << std::endl;
        std::cout << "    mprow QR     = " << _mprowQR   << std::endl;
        std::cout << "    npcol QR     = " << _npcolQR   << std::endl;
        std::cout << "    mbq          = " << _mbq       << std::endl;
        std::cout << "    nbq          = " << _nbq       << std::endl;
        std::cout << "maxNP            = " << _maxNP     << std::endl;
        std::cout << "rtol             = " << _rtol      << std::endl;
        std::cout << "max iterations   = " << _max_iter  << std::endl;
        if (_matrixInitialized) {
            if (_context == _contextQR) {
                std::cout << "Single blacs context used." << std::endl;
            }
        } else {
            if (_mprow == _mprowQR && _npcol == _npcolQR) {
                std::cout << "Single blacs context will be generated." << std::endl;
            }
        }
        std::cout << "Matrix load time was " << getDistributeMatrixTime() << " seconds. " << std::endl;
        std::cout << "############################################################" << std::endl << std::endl;
    }
}


int
Plh::setMatrixColumn(int j, double *col) {
    return _A->setMatrixColumn(j, col);
}


int
Plh::setRHS(double *rhs) {
    return _b->setMatrixColumn(1, rhs);
}


int
Plh::writeMatrix(std::string filename) {
    _A->write(filename, true);
    return 0;
}


int
Plh::writeRhs(std::string filename) {
    _b->write(filename, true);
    return 0;
}


int
Plh::writeX(std::string filename) {
    _x->write(filename);
    return 0;
}


int
Plh::close() {
    MPI_Barrier(MPI_COMM_WORLD);
    if (_context == _contextQR) {
        MPI_Comm_free(&_row_comm);
        MPI_Comm_free(&_col_comm);
    } else {
        MPI_Comm_free(&_row_comm);
        MPI_Comm_free(&_col_comm);
        MPI_Comm_free(&_row_commQR);
        MPI_Comm_free(&_col_commQR);
    }
    //_FORTRAN(blacs_gridexit)(&_context);
    int one=1;
    _FORTRAN(blacs_exit)(&one);
    return 0;
}


int
Plh::distributeVector(int context, char scope, char top, double *vec, int n, int proc_coord) {
    int zero = 0;
    int one = 1;
    std::string col = std::string("C");
    std::string s(1,scope);
    int myrow, mycol;
    if (context == _context) {
        myrow = _myrow;
        mycol = _mycol;
    } else {
        myrow = _myrowQR;
        mycol = _mycolQR;
    }
    if (col.compare(s) == 0) {
        if (myrow == proc_coord) {
            _FORTRAN(dgebs2d)(&context, &scope, &top, &one, &n, vec, &one);
        } else {
            _FORTRAN(dgebr2d)(&context, &scope, &top, &one, &n, vec, &one, &proc_coord, &mycol);
        }
    } else {
        if (mycol == proc_coord) {
            _FORTRAN(dgebs2d)(&context, &scope, &top, &n, &one, vec, &one);
        } else {
            _FORTRAN(dgebr2d)(&context, &scope, &top, &n, &one, vec, &n, &myrow, &proc_coord);
        }
    }
    return 0;
}


int
Plh::distributeVector(int context, char scope, char top, int *ivec, int n, int proc_coord) {
    int zero = 0;
    int one = 1;
    std::string col = std::string("C");
    std::string s(1,scope);
    int myrow, mycol;
    if (context == _context) {
        myrow = _myrow;
        mycol = _mycol;
    } else {
        myrow = _myrowQR;
        mycol = _mycolQR;
    }
    if (col.compare(s) == 0) {
        if (myrow == proc_coord) {
            _FORTRAN(igebs2d)(&context, &scope, &top, &one, &n, ivec, &one);
        } else {
            _FORTRAN(igebr2d)(&context, &scope, &top, &one, &n, ivec, &one, &proc_coord, &mycol);
        }
    } else {
        if (mycol == proc_coord) {
            _FORTRAN(igebs2d)(&context, &scope, &top, &n, &one, ivec, &one);
        } else {
            _FORTRAN(igebr2d)(&context, &scope, &top, &n, &one, ivec, &n, &myrow, &proc_coord);
        }
    }
    return 0;
}


void
Plh::initMatrix(double *A) {
    _A->initMatrix(A);
}


void
Plh::initRhs(double *b) {
    _b->initMatrix(b);
}


double
Plh::getWallTime(){
    struct timeval time;
    gettimeofday(&time,NULL);
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


void
Plh::loadMatrix(const std::vector< Eigen::Map<Eigen::MatrixXd> >& em) {
    MPI_Barrier(MPI_COMM_WORLD); // For timings
    startTime(TIME_LOADMATRIX);
    _A->loadMatrix(em);
    MPI_Barrier(MPI_COMM_WORLD); // For timings
    stopTime(TIME_LOADMATRIX);
}


void
Plh::loadRhs(const Eigen::Ref<const Eigen::VectorXd> &eb){
    startTime(TIME_LOADRHS);
    _b->loadRhs(eb);
    stopTime(TIME_LOADRHS);
}


void
Plh::sub_iteration_output(int iqr) {
    // _rnorm2 = residual2Norm();
    if (_mypid == 0) {
        //if (_iter%HEADER_INCR == 0) {
        //    header();
        //}
        printf("#%5d %8d %12.4e %8d %12.4e %8d %12.4e %12.4e %8d\n",
            _iter, _nP, _wmax.x, _wmax.i, _zmin.x, _zmin.i, _rnorm2, _alpha, iqr);
    }
}


void
Plh::iteration_output() {
    // _rnorm2 = residual2Norm();
    if (_mypid == 0) {
        if (_iter%HEADER_INCR == 0) {
            header();
        }
        printf("%6d %8d %12.4e %8d %12.4e %8d %12.4e\n",
            _iter, _nP, _wmax.x, _wmax.i, _zmin.x, _zmin.i, _rnorm2);
    }
}


void
Plh::header() {
    std::cout << " iter  ";
    std::cout << "   nP    ";
    std::cout << "   wmax(A)   ";
    std::cout << "iwmax(A) ";
    std::cout << "   zmin(Q)   ";
    std::cout << "izmin(Q) ";
    std::cout << "    rnorm    ";
    std::cout << "  alpha(Q)   ";
    std::cout << " iqr(Q)  ";
    std::cout << std::endl;

    std::cout << "------ ";
    std::cout << "-------- ";
    std::cout << "------------ ";
    std::cout << "-------- ";
    std::cout << "------------ ";
    std::cout << "-------- ";
    std::cout << "------------ ";
    std::cout << "------------ ";
    std::cout << "-------- ";
    std::cout << std::endl;;
}


void
Plh::writeSolution(bool compact) {
    _x->write("x.plh", compact);
    writeSet("set.plh");
}


void
Plh::printTimes(bool debug) {
    double times_max[N_TIMES];
    MPI_Allreduce(_wallclock_total, times_max, N_TIMES, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    for (int i=0; i<N_TIMES; i++) {
        _wallclock_total[i] = times_max[i];
    }

    if (debug) {
        double w_maxloc = _w->getMaxTime(SCDBL_TIME_GETMAXLOC);
        double w_max    = _w->getMaxTime(SCDBL_TIME_GETMAX);
        double z_minloc = _zQR->getMaxTime(SCDBL_TIME_GETMINLOC);
        double z_min    = _zQR->getMaxTime(SCDBL_TIME_GETMIN);

        if (_mypid == 0) {
            std::cout << std::endl;
            std::cout << "Wallclock Times (seconds):"                                      << std::endl;
            std::cout << "    gradf            : " << times_max[TIME_GRADF]                << std::endl;
            std::cout << "    mult gradf       : " << times_max[TIME_MULT_GRADF]           << std::endl;
            std::cout << "    pdormq gradf     : " << times_max[TIME_PDORMQR_GRADF]        << std::endl;
            std::cout << "    updateQR         : " << times_max[TIME_UPDATEQR]             << std::endl;
            std::cout << "    qr               : " << times_max[TIME_PDORMQR]              << std::endl;
            std::cout << "    copyRedist       : " << times_max[TIME_COPYREDIST]           << std::endl;
            std::cout << "    copyxQRtox       : " << times_max[TIME_COPYXQRTOX]           << std::endl;
            std::cout << "    mcopyQtoA        : " << times_max[TIME_MCOPYQTOA]            << std::endl;
            std::cout << "    updateQtb        : " << times_max[TIME_UPDATEQTB]            << std::endl;
            std::cout << "    updateX          : " << times_max[TIME_UPDATEX]              << std::endl;
            std::cout << "    SolveR           : " << times_max[TIME_SOLVER]               << std::endl;
            std::cout << "    pdtrsv           : " << times_max[TIME_PDTRSV]               << std::endl;
            std::cout << "    moveFromPToZ     : " << times_max[TIME_MOVEFROMPTOZ]         << std::endl;
            std::cout << "    initnnls         : " << times_max[TIME_INIT_NNLS]            << std::endl;
            std::cout << "    writeResidual    : " << times_max[TIME_WRITE_RESIDUAL]       << std::endl;
            std::cout << "    LD Check         : " << times_max[TIME_LDCHECK]              << std::endl;
            std::cout << "    PZ Check         : " << times_max[TIME_PZCHECK]              << std::endl;
            std::cout << "    nextVector       : " << times_max[TIME_NEXT_VECTOR]          << std::endl;
            std::cout << "    rejectVector     : " << times_max[TIME_REJECT_VECTOR]        << std::endl;
            std::cout << "    getMax           : " << times_max[TIME_GETMAX]               << std::endl;
            std::cout << "    _w MaxLoc        : " << w_maxloc                             << std::endl;
            std::cout << "    _w Max           : " << w_max                                << std::endl;
            std::cout << "    _zQR MinLoc      : " << z_minloc                             << std::endl;
            std::cout << "    _zQR Min         : " << z_min                                << std::endl;
            std::cout << "    loadMatrix       : " << times_max[TIME_LOADMATRIX]           << std::endl;
            std::cout << "    loadRhs          : " << times_max[TIME_LOADRHS]              << std::endl;
            std::cout << "    get x Eigen      : " << times_max[TIME_GET_SOLUTION]         << std::endl;
            std::cout << "    Down Date        : " << times_max[TIME_DOWNDATE]             << std::endl;
            std::cout << "    Loop Total       : " << times_max[TIME_MAIN_LOOP]            << std::endl;
        }
    } else {
        if (_mypid == 0) {
            std::cout << std::endl;
            std::cout << "Wallclock Times (seconds):"                               << std::endl;
            std::cout << "    Solver            : " << times_max[TIME_MAIN_LOOP]    << std::endl;
            std::cout << "        gradf         : " << times_max[TIME_GRADF]        << std::endl;
            std::cout << "        QR            : " << times_max[TIME_UPDATEQR]     << std::endl;
            std::cout << "        Down Date     : " << times_max[TIME_DOWNDATE]     << std::endl;
            std::cout << "    Distribute Matrix : " << times_max[TIME_LOADMATRIX]   << std::endl;
            std::cout << "    Distribute RHS    : " << times_max[TIME_LOADRHS]      << std::endl;
            std::cout << "    get x Eigen       : " << times_max[TIME_GET_SOLUTION] << std::endl;
        }
    }
}


std::string
Plh::intToString(int i) {
    std::ostringstream ss;
    ss << i;
    return ss.str();
}


std::vector<int>
Plh::factor(int n) {    
    std::vector<int> factors;
    for(int i=2; i*i <= n; ++i) {
        while(n%i==0) {
            factors.push_back(i);
            n /= i;
        }
    }
    if(n>1)
        factors.push_back(n);
    return factors;
}



void
Plh::defaultProcGrid() {
    _mprow = _nprocs;
    _npcol = 1;
    /*
    if (_m == 0 || _n ==0 || _nprocs == 1) {
        _mprow = _nprocs;
        _npcol = 1;
    } else {
        std::vector<int> factors = factor(_nprocs);
        if (_mypid == 0) {
            //std::cout << "Factors of nprocs = " << _nprocs << " are: ";
            //std::copy(factors.begin(), factors.end(), std::ostream_iterator<int>(std::cout, ";"));
            std::cout << std::endl;
        }
        int mult, i;
        if (_m >= _n) {
            mult = _m/_n;
            if (mult*_n != _m) mult += 1;
            i = factors.size() - 1;
            _mprow = factors[i];
            _npcol = _nprocs / _mprow;
            while(_mprow < mult*_npcol && i > 0) {
                i--;
                _mprow *= factors[i];
                _npcol = _nprocs / _mprow;
            }
        } else {
            mult = _n/_m;
            if (mult*_m != _n) mult += 1;
            i = factors.size() - 1;
            _npcol = factors[i];
            _mprow = _nprocs / _npcol;
            while(_npcol < mult*_mprow && i > 0) {
                i--;
                _npcol *= factors[i];
                _mprow = _nprocs / _npcol;
            }
            // Error on the side of larger _mprow -> Back out the last factor
            if (_npcol > mult*_mprow && i < factors.size()-1) {
                _npcol = _npcol/factors[i+1];
                _mprow = _nprocs / _npcol;
            }
        }
    }
    */
}


// Move 1 value at a time. Probably can do better, but only called once so why bother.
Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
Plh::getSolution() {
    startTime(TIME_GET_SOLUTION);
    //std::cout << "_mypid = " << _mypid << "Entering getSolution" << std::endl;
    Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1> x(_nEigenSubDomains);
    int jg, ps, pe;
    double val;
    MPI_Status status;
    MPI_Request sendrequest;

    // Send scalapack array. Non-blocking
    if (_myrow == 0) {
        double * matrix = _x->getMatrix();
        for (int j=1; j<=_x->getSizelocal(); j++) {
            jg = _x->getGlobalColIndex(j); // Argument needs Fortran index. Returns Fortran index.
            pe  = getEigenProc(jg-1);      // Argument needs C index. Returns C index.
            //std::cout << "_mypid = " << _mypid << ", Sending jg = " << jg << " to processor " << pe << std::endl;
            MPI_Isend(&(matrix[j-1]), 1, MPI_DOUBLE, pe, jg, MPI_COMM_WORLD, &sendrequest);
        }
    }

    //MPI_Barrier(MPI_COMM_WORLD);
    // Receive and stuff into Eigen array. Blocking.
    for (int isd=0; isd<_nEigenSubDomains; isd++) {
        x[isd] = Eigen::VectorXd::Zero(_eigenSubDomainSize[isd]);
        for (int j=0; j<_eigenSubDomainSize[isd]; j++) {
            jg = j + _eigenSubdomainStart[isd] + 1; // Fortran global index
            ps = _x->getProc(1,jg); // Scalapack processor for global col jg
            //std::cout << "_mypid = " << _mypid << ", Receiving jg = " << jg << " from processor " << ps << std::endl;
            MPI_Recv(&val, 1, MPI_DOUBLE, ps, jg, MPI_COMM_WORLD, &status);
            x[isd][j] = val;
        }
    }


    //std::cout << "_mypid = " << _mypid << "Exiting getSolution" << std::endl;
    stopTime(TIME_GET_SOLUTION);
    return x;
}


void
Plh::write(std::string filename, Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1> & x) {
    //std::cout << "_mypid = " << _mypid << "Writing solution. size = " << x.size() << std::endl;
    std::string fname = filename + "." + intToString(_mypid);
    FILE * f = fopen(fname.c_str(), "w");
    int size = x.size();
    for (int isd=0; isd<x.size(); isd++) {
        for (int j=0; j< x[isd].rows(); j++) {
            fprintf(f, "%24.16e \n", x[isd][j]);
        }
    }
    fclose(f);
}


int
Plh::getEigenProc(int j) {
    if (_nprocs == 1 ) return 0;
    int p = 0;
    int nstart = _eigenColsPerProc[p];
    while (j >= nstart && p < _nprocs-1) {
        p++;
        nstart += _eigenColsPerProc[p];
    }
    return p;
}


void
Plh::testCommunicators() {
    int zero=0;
    int buf = _mypid;
    if (_myrow == 0) {
        buf = 1;
    } else {
        buf = 0;
    }   
    int st = MPI_Bcast(&buf, 1, MPI_INT, 0, _col_comm);
    std::cout << "Col Test: _mypid = " << _mypid << ", _myrow = " << _myrow << ", buf = " << buf << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    if (_mycol == 0) {
        buf = 1;
    } else {
        buf = 0;
    }   
    st = MPI_Bcast(&buf, 1, MPI_INT, 0, _row_comm);
    std::cout << "Row Test: _mypid = " << _mypid << ", _mycol = " << _mycol << ", buf = " << buf << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
}


void
Plh::setDownDateMask() {
    if (_wmask == NULL) {
        _wmask = new SCDoubleMatrix(_context, 1, _n, _mb, _nb);
    }
    _ddmask = true;
}


void
Plh::setColumnNormMask( double cnfac) {
    if (_wmask == NULL) {
        _wmask = new SCDoubleMatrix(_context, 1, _n, _mb, _nb);
    }
    _colnorms = new SCDoubleMatrix(_context, 1, _n, _mb, _nb);
    _cnmask = true;
}
