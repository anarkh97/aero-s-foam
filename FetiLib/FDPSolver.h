//
// Created by Michel Lesoinne on 10/31/17.
//

#ifndef FEM_FDPSOLVER_H
#define FEM_FDPSOLVER_H

#include <complex>
#include <memory>
#include <experimental/propagate_const>
#ifdef USE_MPI
#include <mpi.h>
namespace FetiLib {
using Com = MPI_Comm;
inline auto default_com() { return MPI_COMM_WORLD; }
}
#else
namespace FetiLib {
using Com = int;
inline auto default_com() { return -1; }
}
#endif

#include <Eigen/Dense>

#include <FetiLib/DOFInfo.h>
#include "Subdomain.h"

namespace FetiLib {

namespace tpl {


class DPSImpl;

template <typename O>
using const_enforcing_unique_ptr = std::experimental::propagate_const<std::unique_ptr<O>>;

template <typename T>
class DPSolver {
	DPSolver(std::vector<Subdomain<T>> subdomains, Com communicator = default_com());
	bool solve(Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> rhs,
	           Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> solution);
private:
	/** \brief Pointer to the opaque implementation of the solver.
	 * \details The pointer is wrapped to enforce const correctness.
	 */
	const_enforcing_unique_ptr<DPSImpl> pImpl;
};

}

using FDPSolver = tpl::DPSolver<double>;
using FDPHSolver = tpl::DPSolver<std::complex<double>>;

}

#endif //FEM_FDPSOLVER_H
