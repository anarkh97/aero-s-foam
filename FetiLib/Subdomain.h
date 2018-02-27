//
// Created by Michel Lesoinne on 2/1/18.
//

#ifndef FEM_SUBDOMAIN_H
#define FEM_SUBDOMAIN_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "DOFInfo.h"

namespace FetiLib {

template <typename S>
class Subdomain {
public:
	Subdomain(DOFInfo dofInfo, Eigen::Matrix<double,3,Eigen::Dynamic> X, Eigen::SparseMatrix<S> K);
	const DOFInfo &getDOFInfo() { return dofInfo; }
	const Eigen::Matrix<double, 3, -1> &getX() const {
		return X;
	}

	const Eigen::SparseMatrix<S> &getK() const {
		return K;
	}
private:
	DOFInfo dofInfo;
	Eigen::Matrix<double,3,Eigen::Dynamic> X;
	Eigen::SparseMatrix<S> K;
};

}


#endif //FEM_SUBDOMAIN_H
