//
// Created by Michel Lesoinne on 2018-12-14.
//

#ifndef FEM_CHOLMODIMP_H
#define FEM_CHOLMODIMP_H

#ifdef USE_CHOLMOD
#include "cholmod.h"
#else
struct cholmod_sparse{};
struct cholmod_factor{};
#endif //USE_CHOLMOD

class Connectivity;
class EqNumberer;

class CholmodImp {
public:
	CholmodImp(const Connectivity *cn, const EqNumberer *_dsa, bool isComplex);
	~CholmodImp();

	long memorySize() const;
private:
	cholmod_sparse *A;
	cholmod_factor *L = nullptr;
	size_t n;
	size_t nnz;

};


#endif //FEM_CHOLMODIMP_H
