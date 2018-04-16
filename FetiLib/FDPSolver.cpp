//
// Created by Michel Lesoinne on 10/31/17.
//

#include <mpi.h>
#include <Driver.d/SComm.h>
#include <set>
#include <Comm.d/BaseCommunicator.h>
#include <Comm.d/MPICompatTraits.h>
#include "FDPSolver.h"
#include <Feti.d/Feti.h>
#include "ConcreteSub.h"

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
	Segmenter(global_node_index nNodes, int nProc) : quotient(nNodes/nProc), remainder(nNodes%nProc) {}
	int getProc(global_node_index node) {
		auto pr = node/quotient;
		auto r = node % quotient;
		if(r < pr && r < remainder)
			return pr-1;
		return pr;
	}
private:
	global_node_index quotient;
	global_node_index remainder;
};

/*
 getSharedNodes(nodeToSub, subToNode);

 makeCorners();// Corners for FETI-DP

 getSharedDOFs();
 preProcessMPCs();//Multi-Point Constraint

 getSharedFSIs();

 getSharedMPCs();

 paralApply(numSub, subDomain, &BaseSub::mergeInterfaces);
 paralApply(numSub, subDomain, &GenSubDomain<Scalar>::applySplitting);

 //paralApply(numSub, subDomain, &GenSubDomain<Scalar>::initSrc);
 makeInternalInfo();

 makeNodeInfo();
 */

/** \brief Function creating useful connectivities from a matrix structure
 *
 */
template <typename T>
Connectivity getConnectivities(const Subdomain<T> &subdomain) {

}

/**
 * \details makeSubs needs to equip ConcreteSub with:
 * nodeToNode connectivity.
 *
 * @tparam T
 * @param subdomains
 * @return
 */
template <typename T>
std::vector<std::unique_ptr<ConcreteSub<T>>> makeSubs(const std::vector<Subdomain<T>> &subdomains) {

}

namespace tpl {

/** \brief Obtain the set of nodes involved in a set of DOFs.
 *
 * @param dofInfo An array of pair of global node number, DofType for all the DOFs in this CPU.
 * @return The set of all global node numbers that appear in dofInfo.
 */
std::set<global_node_index> globalNodeSet(const DOFInfo &dofInfo) {
	// Could be done in one line, using transform iterators.
	std::set<global_node_index> nodes;
	for(const auto &info : dofInfo)
		nodes.insert(info.first);
	return nodes;
}

SharedDataInfo
buildSharedDataInfo(const std::set<global_node_index> &glNodes, const BaseCommunicator &communicator) {
	auto it = std::max_element(glNodes.begin(), glNodes.end());
	global_node_index largestNode = (it == glNodes.end()) ? 0 : *it;
	global_node_index numGlNodes = communicator.globalMax(largestNode+1);

    throw "Unfinished.";
}

//SComm buildSComm() {
//	return SComm{0};
//}

class OptionSet {
public:
	void addOption(const char *name, double v) {
		doubleOptions.insert({std::string{name}, v});
	}

	double getDouble(const std::string name, double defaultValue) {
		auto it = doubleOptions.find(name);
		return it == doubleOptions.end() ? defaultValue : it->second;
	}

private:
	std::map<std::string, double> doubleOptions;
};

/** \brief The base class of FetiDP(H) implementation */
class DPSImpl {
public:
	DPSImpl(Com communicator) : communicator(CommunicatorHandle(communicator)) {}

	void setOption(const char *optionName, double v) { options.addOption(optionName, v); }

protected:
	const BaseCommunicator communicator;
	OptionSet options;
};

/** \brief Implementation class for the FetiDPSolver. */
template <typename T>
class FetiDP : public DPSImpl {
	FetiDP(std::vector<Subdomain<T>> subdomains, Com communicator);
private:
	std::vector<Subdomain<T>> subdomains;
	std::unique_ptr<GenFetiDPSolver<T>> fetiSolver;
};

template<typename T>
FetiDP<T>::FetiDP(std::vector<Subdomain<T>> subdomains, Com communicator) :
		DPSImpl(communicator),
		subdomains(std::move(subdomains)){
}
namespace {

template <typename T>
std::unique_ptr<FetiSub<T>> makeFetiSub(int localSubIndex, const Subdomain<T> &subData) {
	auto sub = std::make_unique<FetiSub<T>>();
	SComm *sc = nullptr;
//	SComm *sc = new SComm(nConnect[NESubMap[subI]], connectedDomain[NESubMap[subI]],
//	                      remoteID[NESubMap[subI]], interfNode[NESubMap[subI]]);
	sc->locSubNum = localSubIndex;
	sub->setSComm(sc);

	return std::move(sub);
}

}
template<typename T>
DPSolver<T>::DPSolver(std::vector<Subdomain<T>> subdomains, Com communicator) {
	if(subdomains.size() > 1)
		throw "Multiple subs/CPU not yet supported.";
	// Build the FetiSubs
	for(int iSub = 0; iSub < subdomains.size(); ++iSub)
		std::unique_ptr<FetiSub<T>> fetiSub{ makeFetiSub(0, subdomains[0]) };
}

template<typename T>
void DPSolver<T>::setOption(const char *optionName, double v) {
	pImpl->setOption(optionName, v);
}

template class DPSolver<double>;
template class DPSolver<std::complex<double>>;

} // namespace tpl

} //namespace FetiLib