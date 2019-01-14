//
// Created by Michel Lesoinne on 2018-12-14.
//

#ifndef FEM_CHOLMODIMP_H
#define FEM_CHOLMODIMP_H

#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>

#ifdef WITH_CHOLMOD
#include "cholmod.h"
#else
struct cholmod_sparse{};
struct cholmod_factor{};
#endif //USE_CHOLMOD

class Connectivity;
class EqNumberer;

class CholmodImp {
public:
	CholmodImp(const SparseData &upperTriangularStructure, bool isComplex);
	~CholmodImp();

	template <typename Scalar>
	void setData(const GenDBSparseMatrix<Scalar> &K);
	void factorize();

	long memorySize() const;
private:
	cholmod_sparse *A;
	cholmod_factor *L = nullptr;
	size_t n;
	size_t nnz;
	bool isComplex;

};


#endif //FEM_CHOLMODIMP_H
