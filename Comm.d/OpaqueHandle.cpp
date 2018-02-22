//
// Created by Michel Lesoinne on 2/5/18.
//

#include <mpi.h>
#include <complex>
#include "OpaqueHandle.h"
#include "MPICompatTraits.h"

OpHandle MaxHandle{MPI_MAX};
OpHandle MinHandle{MPI_MIN};
OpHandle SumHandle{MPI_SUM};
OpHandle ProdHandle{MPI_PROD};
TypeHandle IntHandle{MPI_INT};

template <>
TypeHandle CommTypeTrait<double>::typeHandle() {
	return MPI_DOUBLE;
}

template <>
TypeHandle CommTypeTrait<int>::typeHandle() {
	return MPI_INT;
}

template <>
TypeHandle CommTypeTrait<long>::typeHandle() {
	return MPI_LONG;
}

template <>
TypeHandle CommTypeTrait<std::complex<double>>::typeHandle() {
	return MPI_COMPLEX;
}

RankHandle RankHandle::Any{-1};