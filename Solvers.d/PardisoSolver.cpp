//
// Created by Michel Lesoinne on 2019-01-14.
//

#include "PardisoSolver.h"

#ifdef WITH_PARDISO

#include <mkl_pardiso.h>
#include <mkl_types.h>

class PardisoSolver {

};

template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getPardiso(const Connectivity *cn, const EqNumberer *_dsa) {
	return {nullptr, nullptr};
}

#else
template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getPardiso(const Connectivity *cn, const EqNumberer *_dsa) {
	return {nullptr, nullptr};
}

#endif

#include <complex>

template std::pair<GenSolver<double> *, GenSparseMatrix<double> *>
getPardiso<double>(const Connectivity *cn, const EqNumberer *_dsa);

template std::pair<GenSolver<std::complex<double>> *, GenSparseMatrix<std::complex<double>> *>
getPardiso<std::complex<double>>(const Connectivity *cn, const EqNumberer *_dsa);
