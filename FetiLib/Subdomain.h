//
// Created by Michel Lesoinne on 2/1/18.
//

#ifndef FEM_SUBDOMAIN_H
#define FEM_SUBDOMAIN_H

#include <Eigen/Sparse>
#include "DOFInfo.h"

namespace FetiLib {

template <typename S>
class Subdomain {
public:
	Subdomain(DOFInfo dofInfo, Eigen::SparseMatrix<S> K);
private:
	DOFInfo dofInfo;
	Eigen::SparseMatrix<S> K;
};

}


#endif //FEM_SUBDOMAIN_H
