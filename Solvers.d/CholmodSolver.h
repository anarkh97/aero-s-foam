//
// Created by Michel Lesoinne on 2018-12-14.
//

#ifndef FEM_CHOLDMOD_H
#define FEM_CHOLDMOD_H

#include <Utils.d/Connectivity.h>
#include "SolverCntl.h"


template <typename Scalar>
class CholmodSolver {
public:
	CholmodSolver(const Connectivity *cn, const EqNumberer *_dsa, double _tol,
	         const SolverCntl &_scntl, int _ngrbm);
private:

};


#endif //FEM_CHOLDMOD_H
