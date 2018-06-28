//
// Created by Michel Lesoinne on 2/5/18.
//


#include <complex>
#include "OpaqueHandle.h"
#include "MPICompatTraits.h"

#ifdef USE_MPI
#include <mpi.h>
OpHandle MaxHandle{MPI_MAX};
OpHandle MinHandle{MPI_MIN};
OpHandle SumHandle{MPI_SUM};
OpHandle ProdHandle{MPI_PROD};
TypeHandle IntHandle{MPI_INT};

template <>
TypeHandle CommTypeTrait<double>::typeHandle() {
	return TypeHandle{MPI_DOUBLE};
}

template <>
TypeHandle CommTypeTrait<char>::typeHandle() {
	return TypeHandle{MPI_CHAR};
}

template <>
TypeHandle CommTypeTrait<int>::typeHandle() {
	return TypeHandle{MPI_INT};
}

template <>
TypeHandle CommTypeTrait<long>::typeHandle() {
	return TypeHandle{MPI_LONG};
}

template <>
TypeHandle CommTypeTrait<std::complex<double>>::typeHandle() {
	return TypeHandle{MPI_COMPLEX};
}

#else

OpHandle MaxHandle{1};
OpHandle MinHandle{2};
OpHandle SumHandle{3};
OpHandle ProdHandle{4};

template <>
TypeHandle CommTypeTrait<double>::typeHandle() {
    return TypeHandle{sizeof(double)};
}

template <>
TypeHandle CommTypeTrait<char>::typeHandle() {
    return TypeHandle{sizeof(char)};
}

template <>
TypeHandle CommTypeTrait<int>::typeHandle() {
    return TypeHandle{sizeof(int)};
}

template <>
TypeHandle CommTypeTrait<long>::typeHandle() {
    return TypeHandle{sizeof(long)};
}

template <>
TypeHandle CommTypeTrait<std::complex<double>>::typeHandle() {
    return TypeHandle{sizeof(std::complex<double>)};
}
#endif // USE_MPI

RankHandle RankHandle::Any{-1};