//
// Created by Michel Lesoinne on 10/31/17.
//

#include <mpi.h>
#include <Driver.d/SComm.h>
#include <set>
#include <Comm.d/BaseCommunicator.h>
#include <Comm.d/MPICompatTraits.h>
#include "FDPSolver.h"

namespace FetiLib {

struct SharedDataInfo {
	std::vector<int> cpus;
	Connectivity sharedData;
};

/** \brief A reprentation of uniform segmentation the integer interval [0,n).
 * \details The size of each segment can differ only by one, being larger by 1 when the segment index
 * is less than the remainder of the division.
 */
class Segmenter {
	Segmenter(gl_node_t nNodes, int nProc) : quotient(nNodes/nProc), remainder(nNodes%nProc) {}
	int getProc(gl_node_t node) {
		auto pr = node/quotient;
		auto r = node % quotient;
		if(r < pr && r < remainder)
			return pr-1;
		return pr;
	}
private:
	gl_node_t quotient;
	gl_node_t remainder;
};

namespace tpl {

/** \brief Obtain the set of nodes involved in a set of DOFs.
 *
 * @param dofInfo An array of pair of global node number, DofType for all the DOFs in this CPU.
 * @return The set of all global node numbers that appear in dofInfo.
 */
std::set<gl_node_t> globalNodeSet(const DOFInfo &dofInfo) {
	// Could be done in one line, using transform iterators.
	std::set<gl_node_t> nodes;
	for(const auto &info : dofInfo)
		nodes.insert(info.first);
	return nodes;
}

SharedDataInfo
buildSharedDataInfo(const std::set<gl_node_t> &glNodes, const BaseCommunicator &communicator) {
	auto it = std::max_element(glNodes.begin(), glNodes.end());
	gl_node_t largestNode = (it == glNodes.end()) ? 0 : *it;
	gl_node_t numGlNodes = communicator.globalMax(largestNode+1);

    throw "Unfinished.";
}

//SComm buildSComm() {
//	return SComm{0};
//}

/** \brief The base class of FetiDP(H) implementation */
class DPSImpl {
public:
	DPSImpl(Com communicator) : communicator(CommunicatorHandle(communicator)) {}
	const BaseCommunicator communicator;
};

/** \brief Implementation class for the FetiDPSolver. */
template <typename T>
class FetiDP : public DPSImpl {
	FetiDP(std::vector<Subdomain<T>> subdomains, Com communicator);
private:
	std::vector<Subdomain<T>> subdomains;
};

template<typename T>
FetiDP<T>::FetiDP(std::vector<Subdomain<T>> subdomains, Com communicator) :
		DPSImpl(communicator),
		subdomains(std::move(subdomains)){

}

template<typename T>
DPSolver<T>::DPSolver(std::vector<Subdomain<T>> subdomains, Com communicator) {

}

template class DPSolver<double>;
template class DPSolver<std::complex<double>>;

} // namespace tpl

} //namespace FetiLib