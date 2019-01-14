//
// Created by Michel Lesoinne on 2018-12-14.
//

#include <complex>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include "CholmodImp.h"

namespace {
#ifdef WITH_CHOLMOD
cholmod_common *getCholmodCommon() {
	static cholmod_common c;
	static bool is_initialized = [] {
		cholmod_start(&c);
		return true;
	} ();
	return &c;
}
#endif // WITH_CHOLMOD
}

CholmodImp::CholmodImp(const SparseData &upperTriangularStructure, bool isComplex)
	: isComplex(isComplex)
{
#ifdef WITH_CHOLMOD
	n = upperTriangularStructure.numCol();
	nnz = upperTriangularStructure.nnz();
	auto c = getCholmodCommon();
	A = cholmod_allocate_sparse(
		n,
		n,
		nnz,
		false,
		true,
		1, // upper triangular
		isComplex ? CHOLMOD_COMPLEX : CHOLMOD_REAL,
		c
	);
	// Create the structure of A
	auto colptr = static_cast<int*>(A->p);
	auto rowind = static_cast<int*>(A->i);
	auto val = static_cast<double*>(A->x);

	auto &colPointers = upperTriangularStructure.colPointers();
	auto &rowIndices = upperTriangularStructure.rowIndices();
	const auto basis = colPointers[0]; // 0 or 1 depending on the SparseData implementation.
	for (int i = 0; i < colPointers.size(); ++i)
		colptr[i] = colPointers[i] - basis;
	for (int i = 0; i < rowIndices.size(); ++i)
	rowind[i] = rowIndices[i] - basis;

	L = cholmod_analyze(A, c);
//	cholmod_factorize(A, L, c);
#endif // WITH_CHOLMOD
}

CholmodImp::~CholmodImp()
{
#ifdef WITH_CHOLMOD
	auto c = getCholmodCommon();
	cholmod_free_factor (&L, c) ; /* free matrices */
	cholmod_free_sparse (&A, c) ;
#endif // WITH_CHOLMOD
}

long CholmodImp::memorySize() const
{
	return 0;
}

void CholmodImp::factorize()
{
	auto c = getCholmodCommon();
	if (!L)
		L = cholmod_analyze(A, c);
	cholmod_factorize(A, L, c);
}

template<typename Scalar>
void CholmodImp::setData(const GenDBSparseMatrix <Scalar> &K)
{
	auto colptr = static_cast<int*>(A->p);
	auto rowind = static_cast<int*>(A->i);
	auto val = static_cast<Scalar*>(A->x);
	auto &colPointers = K.colPointers();
	auto &rowIndices = K.rowIndices();
	const auto basis = colPointers[0]; // 0 or 1 depending on the SparseData implementation.

	const auto &Kvalues = K.values();
	for (size_t i = 0 ; i < Kvalues.size(); ++i)
		val[i] = Kvalues[i];
}

template void CholmodImp::setData<double>(const GenDBSparseMatrix <double> &K);
template void CholmodImp::setData(const GenDBSparseMatrix <std::complex<double>> &K);

void test() {
//	// Step 2: translate sparse matrix into cholmod format
//	//----------------------------------------------------
//
//	cholmod_sparse* A = cholmod_allocate_sparse(
//		n, n, nnz,    // Dimensions and number of non-zeros
//		false,        // Sorted = false
//		true,         // Packed = true
//		1,            // stype (-1 = lower triangular, 1 = upper triangular)
//		CHOLMOD_REAL, // Entries are real numbers
//		&c
//	);
//
//	int* colptr = (int*)A->p;
//	int* rowind = (int*)A->i;
//	double* val = (double*)A->x;
//
//	// Convert Geogram Matrix into CHOLMOD Matrix
//	index_t count = 0 ;
//	for(index_t i=0; i<n; ++i) {
//		colptr[i] = int(count);
//		NLRowColumn& Ri = MM.row[i];
//		for(index_t jj=0; jj<Ri.size; ++jj) {
//			const NLCoeff& C = Ri.coeff[jj];
//			index_t j = C.index;
//			if(j <= i) {
//				val[count] = C.value;
//				rowind[count] = int(j);
//				++count;
//			}
//		}
//	}
//	geo_assert(count == nnz);
//	colptr[n] = int(nnz);
//
//	/*
//	geo_assert(cholmod_check_sparse(A,&c) != 0);
//	if(n < 10) {
//		cholmod_write_sparse(stdout,A,NULL,NULL,&c);
//	}
//	cholmod_print_sparse(A,"A",&c);
//	*/
//
//	// Step 2: construct right-hand side
//	cholmod_dense* b = cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, &c) ;
//	Memory::copy(b->x, b_in, n * sizeof(double)) ;
//
//	// Step 3: factorize and solve
//	cholmod_factor* L = cholmod_analyze(A, &c) ;
//	geo_debug_assert(cholmod_check_factor(L,&c) != 0);
//
//	if(!cholmod_factorize(A, L, &c)) {
//		std::cerr << "COULD NOT FACTORIZE !!!" << std::endl;
//	}
//
//	cholmod_dense* x = cholmod_solve(CHOLMOD_A, L, b, &c) ;
//	Memory::copy(x_out, x->x, n * sizeof(double)) ;
//
//	// Step 4: cleanup
//	cholmod_free_factor(&L, &c) ;
//	cholmod_free_sparse(&A, &c) ;
//	cholmod_free_dense(&x, &c) ;
//	cholmod_free_dense(&b, &c) ;
}