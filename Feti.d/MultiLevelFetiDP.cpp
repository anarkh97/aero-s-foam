//
// Created by Michel Lesoinne on 6/25/18.
//

#include <Paral.d/MDDynam.h>
#include <Element.d/MatrixElement.d/MatrixElement.h>
#include "Feti.h"
#include <Driver.d/DecDomain.h>

extern double t5;

template <typename Scalar>
void GenFetiDPSolver<Scalar>::makeMultiLevelDP(const Connectivity *subToCorner) {
	std::cerr << "using FETI-DP solver for coarse problem\n";
	Domain *coarseDomain = new Domain();
	coarseDomain->solInfo().solvercntl = fetiInfo->coarse_cntl;
	int numSub = subToCorner->csize();
	coarseDomain->setNumElements(subToCorner->csize());

	Connectivity *decCoarse = this->cpuToSub;
	char fn[65];
	FILE *f;
	sprintf(fn,"decomposition.%d",numSub);
	if ((f = fopen(fn,"r")) != NULL) {
		if(verboseFlag) filePrint(stderr, " ... Reading Decomposition from file %s ...\n", fn);
		decCoarse = new Connectivity(f,numSub);
		fclose(f);
	}

	Connectivity *elemToSubCoarse = decCoarse->reverse();
	Connectivity *CPUMapCoarse = this->cpuToSub->transcon(elemToSubCoarse);
	delete elemToSubCoarse;

	std::vector<int> pointer(this->glNumSub+1, 0);

	for(int i = 0; i < this->nsub; ++i) {
		int s = this->subdomains[i]->subNum();
		int nc = subToCorner->num(s);
		int n = nc;
		if (fetiInfo->augmentimpl == FetiInfo::Primal)
			for(int iNeighb = 0; iNeighb < this->subdomains[i]->numNeighbors(); ++iNeighb)
				if(this->subdomains[i]->isEdgeNeighbor(iNeighb) &&
				   this->subdomains[i]->edgeDofs[iNeighb].count())
					n++;
//		std::vector<int> elem(n);
		// TODO Fix this leak!!!!
		int *elem = new int[n];
		for(n = 0; n < nc; n++)
			elem[n] = (*subToCorner)[s][n];
		if (fetiInfo->augmentimpl == FetiInfo::Primal) {
			int iEdgeN = 0;
			for(int iNeighb = 0; iNeighb < this->subdomains[i]->numNeighbors(); ++iNeighb) {
				if(this->subdomains[i]->isEdgeNeighbor(iNeighb)) {
					if(this->subdomains[i]->edgeDofs[iNeighb].count())
						elem[n++] = augOffset + (*(this->subToEdge))[s][iEdgeN];
					iEdgeN++;
				}
			}
		}
//        coarseDomain->addElem(s, 0, subToCorner->num(s), (*subToCorner)[s]); // 0 is a "matrix" element
		pointer[s] = n;
		coarseDomain->addElem(s,0,n,elem);//.data());
	}

	if(verboseFlag) filePrint(stderr, " ... Assemble Kcc solver            ...\n");
	t5 -= getTime();
	paralApplyToAll(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::multKcc); // create the local Kcc^*
	t5 += getTime();
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

	CoordSet& nodes = coarseDomain->getNodes();
	double xyz[3];
	for(int i = 0; i < this->nsub; i++) {
		int s = this->subdomains[i]->subNum();
		int numCorner = this->subdomains[i]->numCorners();
		const auto &localCornerNodes = this->subdomains[i]->getLocalCornerNodes();
		for(int iCorner = 0; iCorner < numCorner; ++iCorner) {
			const Node *node = this->subdomains[i]->getNodeSet()[localCornerNodes[iCorner]];
			int cornerNum = (*subToCorner)[s][iCorner];
			if (!nodes[cornerNum]) {
				xyz[0] = node->x; xyz[1] = node->y; xyz[2] = node->z;
				nodes.nodeadd(cornerNum, xyz);
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
						xyz[0] = xyz[1] = xyz[2] = 0.0;
						for(int iNode = 0; iNode < sharedNodes.num(iNeighb); iNode++) {
							xyz[0] += subnodes[sharedNodes[iNeighb][iNode]]->x;
							xyz[1] += subnodes[sharedNodes[iNeighb][iNode]]->y;
							xyz[2] += subnodes[sharedNodes[iNeighb][iNode]]->z;
						}
						xyz[0] /= sharedNodes.num(iNeighb);
						xyz[1] /= sharedNodes.num(iNeighb);
						xyz[2] /= sharedNodes.num(iNeighb);
						nodes.nodeadd(edgeNum, xyz);
					}
				}
			}
		}
	}
	/*
	double xyz[3] = { 0.0, 0.0, 0.0 };
	for(int i = 0; i < cornerToSub->csize(); ++i) {
	  nodes.nodeadd(i, xyz);
	}
	*/
	coarseDomain->setNumNodes(cornerToSub->csize());

	// The following loop set the node coordinates... They are needed for the corner maker
	// note that the coordinates of any nodes which are not represented on this mpi process are not correct (zero)
	// however, they shouldn't be needed in any case
	for(int iSub = 0; iSub < this->nsub; ++iSub) {
		int numCorner = this->subdomains[iSub]->numCorners();
		const auto &localCornerNodes = this->subdomains[iSub]->getLocalCornerNodes();
		for(int iCorner = 0; iCorner < numCorner; ++iCorner) {
			const Node *node = this->subdomains[iSub]->getNodeSet()[localCornerNodes[iCorner]];
			int cornerNum = (*subToCorner)[this->subdomains[iSub]->subNum()][iCorner];
			nodes[cornerNum]->x = node->x;
			nodes[cornerNum]->y = node->y;
			nodes[cornerNum]->z = node->z;
		}
	}

#ifdef USE_MPI
	Communicator *structCom = new Communicator(CommunicatorHandle{this->fetiCom->getComm()});
#else
	Communicator *structCom = NULL;
#endif
	GenDecDomain<Scalar> *decCoarseDomain = new GenDecDomain<Scalar>(coarseDomain, structCom, false);

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