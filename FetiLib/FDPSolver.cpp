//
// Created by Michel Lesoinne on 10/31/17.
//

#include "FDPSolver.h"

namespace FetiLib {

namespace tpl {

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