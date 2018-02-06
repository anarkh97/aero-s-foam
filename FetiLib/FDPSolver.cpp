//
// Created by Michel Lesoinne on 10/31/17.
//

#include <Driver.d/SComm.h>
#include <set>
#include "FDPSolver.h"

namespace FetiLib {

namespace tpl {

/** \brief Obtain the set of nodes involved in a set of DOFs. */
std::set<gl_node_t> globalNodeSet(const DOFInfo &dofInfo) {
	// Could be done in one line, using transform iterators.
	std::set<gl_node_t> nodes;
	for(const auto &info : dofInfo)
		nodes.insert(info.first);
	return nodes;
}



//SComm buildSComm() {
//	return SComm{0};
//}

/** \brief The base class of FetiDP(H) implementation */
class DPSImpl {
};

/** \brief Implementation class for the FetiDPSolver. */
template <typename T>
class FetiDP : public DPSImpl {
	FetiDP(std::vector<Subdomain<T>> subdomains);
private:
	std::vector<Subdomain<T>> subdomains;
};

template<typename T>
FetiDP<T>::FetiDP(std::vector<Subdomain<T>> subdomains) : subdomains(std::move(subdomains)){

}

template<typename T>
DPSolver<T>::DPSolver(std::vector<Subdomain<T>> subdomains, Com communicator) {

}

template class DPSolver<double>;
template class DPSolver<std::complex<double>>;

} // namespace tpl

} //namespace FetiLib