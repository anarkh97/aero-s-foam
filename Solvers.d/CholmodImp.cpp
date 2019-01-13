//
// Created by Michel Lesoinne on 2018-12-14.
//

#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include "CholmodImp.h"

namespace {
#ifdef USE_CHOLMOD
cholmod_common *getCholmodCommon() {
	static cholmod_common c;
	static bool is_initialized = [] {
		cholmod_start(&c);
		return true;
	} ();
	return &c;
}
#endif // USE_CHOLMOD
}

CholmodImp::CholmodImp(const SparseData &fortranBasedUpperMatrix, bool isComplex)
{
#ifdef USE_CHOLMOD
	auto c = getCholmodCommon();
	A = cholmod_allocate_sparse(
		n,
		n,
		nnz,
		false,
		true,
		1, // lower triangular
		isComplex ? CHOLMOD_COMPLEX : CHOLMOD_REAL,
		c
	);
	// Fill A

	L = cholmod_analyze(A, c);
	cholmod_factorize(A, L, c);
#endif // USE_CHOLMOD
}

CholmodImp::~CholmodImp()
{
#ifdef USE_CHOLMOD
	auto c = getCholmodCommon();
	cholmod_free_factor (&L, c) ; /* free matrices */
	cholmod_free_sparse (&A, c) ;
#endif // USE_CHOLMOD
}

long CholmodImp::memorySize() const
{
	return 0;
}

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