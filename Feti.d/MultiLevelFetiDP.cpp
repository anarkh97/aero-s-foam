//
// Created by Michel Lesoinne on 6/25/18.
//
#include <sstream>
#include <gsl/span>
#include <Paral.d/MDDynam.h>
#include <Element.d/MatrixElement.d/MatrixElement.h>
#include "Feti.h"
#include <Driver.d/DecDomain.h>
#include <FetiLib/Subdomain.h>
#include <Types.h>

extern double t5;

using std::vector;
using std::unique_ptr;

namespace {
template <typename Scalar>
struct SuperElement {
	std::vector<gl_num_t> nodes;
	std::vector<DofSet> dofs;
	GenAssembledFullM<Scalar> *Kel;
};


/// \brief Make a vector out of a node.
Eigen::Map<const Eigen::Vector3d> v(const Node &node) {
	return Eigen::Map<const Eigen::Vector3d>(&node.x, 3, 1);
}

void formSharedLists(const std::vector<gl_sub_idx> &localSubGlobalIndices, const Connectivity &subToNode)
{

}

template <typename Scalar>
vector<unique_ptr<FetiLib::Subdomain<Scalar>>>
formSubdomains(const std::vector<gl_sub_idx> &localSubGlobalIndices, std::vector<SuperElement<Scalar>> &,
               const Connectivity *subToNode) {
	vector<unique_ptr<FetiLib::Subdomain<Scalar>>> subdomains;



	return subdomains;
}

}

/*
 * Solver building does:
 buildSharedNodeComm(nodeToSub, subToNode);

 makeCorners();// Corners for FETI-DP

 getSharedDOFs();

 preProcessMPCs();//Multi-Point Constraint

 getSharedFSIs();

 getSharedMPCs();
  paralApply(subDomain, &BaseSub::mergeInterfaces);
 paralApply(subDomain, &GenSubDomain<Scalar>::applySplitting);

 //paralApply(subDomain, &GenSubDomain<Scalar>::initSrc);
 makeInternalInfo();

 makeNodeInfo();

 */

/** \brief Get a vector with the 'nodes' of the subdomain created super-element
 *
 * @param sub The subdomain for which the nodes are sought.
 * @param corners Corner nodes of the subdomain, not including edge generalized nodes..
 * @param augOffset
 * @param subdomainEdgeIndices
 * @param withEdgeAugmentation
 * @return
 */
std::vector<gl_num_t> subSuperNodes(const FetiBaseSub &sub,
                                    gsl::span<const gl_num_t> corners,
                                    gl_num_t augOffset,
                                    gsl::span<const gl_num_t> subdomainEdgeIndices,
                                    bool withEdgeAugmentation) {
	std::vector<gl_num_t> nodes {corners.begin(), corners.end()};

	if(withEdgeAugmentation) {
		int iEdgeN = 0;
		for (int iNeighb = 0; iNeighb < sub.numNeighbors(); ++iNeighb) {
			if (sub.isEdgeNeighbor(iNeighb)) {
				if (sub.edgeDofs[iNeighb].count() != 0)
					nodes.push_back(augOffset + subdomainEdgeIndices[iEdgeN]);
				iEdgeN++;
			}
		}
	}
	return nodes;
}

template <typename Scalar>
void GenFetiDPSolver<Scalar>::makeMultiLevelDPNew(const Connectivity *subToCorner) {
	int numSub = subToCorner->csize();

	// Get the decomposition of subdomains into super-subdomains.
	Connectivity *decCoarse = this->cpuToSub;
	std::stringstream fn;
	fn << "decomposition." << numSub << std::endl;
	FILE *f;
	if ((f = fopen(fn.str().c_str(),"r")) != NULL) {
		if(verboseFlag)
			filePrint(stderr, " ... Reading Decomposition from file %s ...\n", fn.str().c_str());
		decCoarse = new Connectivity(f,numSub);
		fclose(f);
	}

	Connectivity *elemToSubCoarse = decCoarse->alloc_reverse();
	Connectivity *CPUMapCoarse = this->cpuToSub->transcon(elemToSubCoarse);
	delete elemToSubCoarse;

	// Build the supersubdomains.
//	std::vector<FetiLib::Subdomain<Scalar>> coarseSubdomains;

	std::vector<SuperElement<Scalar>> superElements;
	superElements.reserve(this->subdomains.size());
	std::vector<gl_sub_idx> localSubGlobalIndices;
	localSubGlobalIndices.reserve(this->subdomains.size());
    // Create a super-element for each subdomain and put it in the coarseSubdomains arrray.
	for(int iSub = 0; iSub < this->nsub; ++iSub) {
		// Get the node numbers.
		auto &sub = *(this->subdomains[iSub]);
		auto globalSubIndex = sub.subNum();
		auto edges = (*(this->subToEdge))[globalSubIndex];
		std::vector<gl_num_t> coarsenodes = subSuperNodes(sub, (*subToCorner)[globalSubIndex], augOffset, edges,
			fetiInfo->augmentimpl == FetiInfo::Primal);
		// Get the DOFs
		std::vector<DofSet> dofs;
		dofs.reserve(coarsenodes.size());
		for(auto &ds : this->subdomains[iSub]->cornerDofs)
			dofs.push_back(ds);
		if (fetiInfo->augmentimpl == FetiInfo::Primal) {
			for(int iNeighb = 0; iNeighb < this->subdomains[iSub]->numNeighbors(); ++iNeighb)
				if (this->subdomains[iSub]->isEdgeNeighbor(iNeighb) &&
				    this->subdomains[iSub]->edgeDofs[iNeighb].count())
					dofs.push_back( this->subdomains[iSub]->edgeDofs[iNeighb] );
		}
		// Get the matrix
		auto Kel = this->subdomains[iSub]->Kcc.get();
		if(dofs.size() != coarsenodes.size())
			std::cerr << "BAD SIZE!\n";
		// Build the super-element
		superElements.push_back({std::move(coarsenodes), std::move(dofs), Kel});
		localSubGlobalIndices.push_back( globalSubIndex );
	}

	// Now form the set of nodes.
	std::map<gl_num_t, Eigen::Vector3d> nodes;
	for(int i = 0; i < this->nsub; i++) {
		Eigen::Vector3d xyz;
		int s = this->subdomains[i]->subNum();
		int numCorner = this->subdomains[i]->numCorners();
		const auto &localCornerNodes = this->subdomains[i]->getLocalCornerNodes();
		for(int iCorner = 0; iCorner < numCorner; ++iCorner) {
			auto lcn = localCornerNodes[iCorner];
			const Node *node = this->subdomains[i]->getNodeSet()[localCornerNodes[iCorner]];
			int cornerNum = (*subToCorner)[s][iCorner];
			nodes.insert({cornerNum, v(*node)});
		}
		if (fetiInfo->augmentimpl == FetiInfo::Primal) { // 020314 JAT
			const CoordSet &subnodes = this->subdomains[i]->getNodeSet();
			Connectivity &sharedNodes = *(this->subdomains[i]->getSComm()->sharedNodes);
			int iEdgeN = 0, edgeNum;
			for(int iNeighb = 0; iNeighb < this->subdomains[i]->numNeighbors(); ++iNeighb) {
				if(this->subdomains[i]->isEdgeNeighbor(iNeighb)) {
					edgeNum = augOffset + (*(this->subToEdge))[s][iEdgeN++];
					if (nodes.find(edgeNum) == nodes.end()) {
						xyz.setZero();
						for(auto node : sharedNodes[iNeighb])
							xyz += v(*subnodes[node]);
						xyz /= sharedNodes.num(iNeighb);
						nodes.insert({ edgeNum, xyz });
					}
				}
			}
		}
	}
	std::cout << "Number of coarse nodes: " << nodes.size() << std::endl;//<< " number "



	formSubdomains(localSubGlobalIndices, superElements, subToCorner);

//	std::vector<FetiBaseSub *> baseSubs(decCoarseDomain->getAllSubDomains(),
//	                                    decCoarseDomain->getAllSubDomains()+decCoarseDomain->getNumSub());
//	paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::makeKccDofsExp2,
//	           decCoarseDomain->getNumSub(),
//	           baseSubs.data(),
//	           augOffset, this->subToEdge);
}

template <typename Scalar>
void GenFetiDPSolver<Scalar>::makeMultiLevelDP(const Connectivity *subToCorner) {
	std::cout << "Test new code does not crash.\n";
	makeMultiLevelDPNew(subToCorner);
	std::cerr << "using FETI-DP solver for coarse problem\n";
	Domain *coarseDomain = new Domain();
	coarseDomain->solInfo().solvercntl = fetiInfo->coarse_cntl;
	int numSub = subToCorner->csize();
	coarseDomain->setNumElements(subToCorner->csize());

	// Get the decomposition of subdomains into super-subdomains.
	Connectivity *decCoarse = this->cpuToSub;
	char fn[65];
	FILE *f;
	sprintf(fn,"decomposition.%d",numSub);
	if ((f = fopen(fn,"r")) != NULL) {
		if(verboseFlag) filePrint(stderr, " ... Reading Decomposition from file %s ...\n", fn);
		decCoarse = new Connectivity(f,numSub);
		fclose(f);
	}

	Connectivity *elemToSubCoarse = decCoarse->alloc_reverse();
	Connectivity *CPUMapCoarse = this->cpuToSub->transcon(elemToSubCoarse);
	delete elemToSubCoarse;

	std::vector<int> pointer(this->glNumSub+1, 0);

	// Create a super-element for each subdomain and put it in the coarseDomain.

	// Phase one only works to get node numbers.
	for(int i = 0; i < this->nsub; ++i) {
		int globalSubIndex = this->subdomains[i]->subNum();
		auto &sub = *(this->subdomains[i]);
		auto edges = (*(this->subToEdge))[globalSubIndex];
		std::vector<gl_num_t> coarsenodes = subSuperNodes(sub, (*subToCorner)[globalSubIndex], augOffset, edges,
		                                                  fetiInfo->augmentimpl == FetiInfo::Primal);
		int n = coarsenodes.size();
		// TODO Fix this leak!!!!
		int *elem = new int[n];
		for(int i = 0; i < n; ++i)
			elem[i] = coarsenodes[i];
		coarseDomain->addElem(globalSubIndex,0,n,elem);//.data());
		pointer[globalSubIndex] = n;
	}

	// Set the dofs and stiffness of each super-element.
	Elemset& elems = coarseDomain->getElementSet();
	for(int i = 0; i < this->nsub; ++i) {
		int s = this->subdomains[i]->subNum();
		int n = pointer[s];
		int nc = subToCorner->num(s);
		DofSet *coarseDofs = new DofSet[n];
		for(n = 0; n < nc; n++)
			coarseDofs[n] = this->subdomains[i]->cornerDofs[n];
		if (fetiInfo->augmentimpl == FetiInfo::Primal) {
			for(int iNeighb = 0; iNeighb < this->subdomains[i]->numNeighbors(); ++iNeighb)
				if(this->subdomains[i]->isEdgeNeighbor(iNeighb) &&
				   this->subdomains[i]->edgeDofs[iNeighb].count())
					coarseDofs[n++] = this->subdomains[i]->edgeDofs[iNeighb];
		}
		((MatrixElement*)elems[s])->setDofs(coarseDofs);
		((MatrixElement*)elems[s])->setStiffness(this->subdomains[i]->Kcc.get());
	}

	// Create the nodeset.
	CoordSet& nodes = coarseDomain->getNodes();
	// Add the corner and edge nodes.
	for(int i = 0; i < this->nsub; i++) {
		Eigen::Vector3d xyz;
		int s = this->subdomains[i]->subNum();
		int numCorner = this->subdomains[i]->numCorners();
		const auto &localCornerNodes = this->subdomains[i]->getLocalCornerNodes();
		for(int iCorner = 0; iCorner < numCorner; ++iCorner) {
			auto lcn = localCornerNodes[iCorner];
//			auto gcn = this->subdomains[i]->
			const Node *node = this->subdomains[i]->getNodeSet()[localCornerNodes[iCorner]];
			int cornerNum = (*subToCorner)[s][iCorner];
			if (!nodes[cornerNum]) {
				xyz[0] = node->x; xyz[1] = node->y; xyz[2] = node->z;
				nodes.nodeadd(cornerNum, xyz.data());
			}
		}
		if (fetiInfo->augmentimpl == FetiInfo::Primal) { // 020314 JAT
			const CoordSet &subnodes = this->subdomains[i]->getNodeSet();
			Connectivity &sharedNodes = *(this->subdomains[i]->getSComm()->sharedNodes);
			int iEdgeN = 0, edgeNum;
			for(int iNeighb = 0; iNeighb < this->subdomains[i]->numNeighbors(); ++iNeighb) {
				if(this->subdomains[i]->isEdgeNeighbor(iNeighb)) {
					edgeNum = augOffset + (*(this->subToEdge))[s][iEdgeN++];
					if (!nodes[edgeNum]) {
						xyz.setZero();
						for(auto node : sharedNodes[iNeighb])
							xyz += v(*subnodes[node]);
						xyz /= sharedNodes.num(iNeighb);
						nodes.nodeadd(edgeNum, xyz.data());
					}
				}
			}
		}
	}


	coarseDomain->setNumNodes(cornerToSub->csize());

#ifdef USE_MPI
	Communicator *structCom = new Communicator(CommunicatorHandle{this->fetiCom->getComm()});
#else
	Communicator *structCom = NULL;
#endif
	GenDecDomain<Scalar> *decCoarseDomain = new GenDecDomain<Scalar>(coarseDomain, structCom, false, true);

	const Connectivity *elemToNode; // JAT 220216
	if (fetiInfo->augmentimpl == FetiInfo::Primal) { // JAT 041114
#ifdef DISTRIBUTED
		this->fetiCom->globalSum(this->glNumSub, pointer.data());
#endif
		int total = 0;
		for(int iSub = 0; iSub < this->glNumSub; ++iSub) {
			int tmp = pointer[iSub];
			pointer[iSub] = total;
			total += tmp;
		}
		pointer[this->glNumSub] = total;

		std::vector<int> target(total, 0);
		for(int i = 0; i < this->nsub; i++) {
			int s = this->subdomains[i]->subNum();
			int nc = subToCorner->num(s);
			int n, k = pointer[s];
			for(n = 0; n < nc; n++)
				target[k+n] = (*subToCorner)[s][n];
			int iEdgeN = 0;
			for(int iNeighb = 0; iNeighb < this->subdomains[i]->numNeighbors(); ++iNeighb) {
				if(this->subdomains[i]->isEdgeNeighbor(iNeighb)) {
					if(this->subdomains[i]->edgeDofs[iNeighb].count())
						target[k+n++] = augOffset + (*(this->subToEdge))[s][iEdgeN];
					iEdgeN++;
				}
			}
		}
#ifdef DISTRIBUTED
		this->fetiCom->globalSum(total, target.data());
#endif
		elemToNode = new Connectivity(this->glNumSub, std::move(pointer), std::move(target));
	} else
		elemToNode = subToCorner;

	decCoarseDomain->setElemToNode(elemToNode);
	if (elemToNode == subToCorner)
		elemToNode = 0;
	decCoarseDomain->setSubToElem(decCoarse);
	decCoarseDomain->setCPUMap(CPUMapCoarse);
	decCoarseDomain->preProcess();
	for(int i=0; i<decCoarseDomain->getNumSub(); ++i) decCoarseDomain->getSubDomain(i)->makeAllDOFs();
	GenMDDynamMat<Scalar> ops;
	if(verboseFlag) filePrint(stderr, " ... Factor Kcc solver              ...\n");
	decCoarseDomain->buildOps(ops, 0.0, 0.0, 1.0);
	coarseInfo = &(decCoarseDomain->solVecInfo());
	KccParallelSolver = ops.dynMat;
	std::vector<FetiBaseSub *> baseSubs(decCoarseDomain->getAllSubDomains(),
	                                    decCoarseDomain->getAllSubDomains()+decCoarseDomain->getNumSub());
	paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::makeKccDofsExp2,
	           decCoarseDomain->getNumSub(),
	           baseSubs.data(),
	           augOffset, this->subToEdge); // JAT 101816
}

template
void GenFetiDPSolver<double>::makeMultiLevelDP(const Connectivity *subToCorner);

template
void GenFetiDPSolver<std::complex<double>>::makeMultiLevelDP(const Connectivity *subToCorner);