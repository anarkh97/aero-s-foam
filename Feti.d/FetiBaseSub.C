//
// Created by Michel Lesoinne on 1/18/18.
//

#include "FetiSub.h"


int
FetiBaseSub::interfLen() const
{
	// Total length for the local interface
	return totalInterfSize;
}

int
FetiBaseSub::halfInterfLen() const
{
	return masterFlagCount;
}


void
FetiBaseSub::setRbmCommSize(int _numRBM, FSCommStructure *pt) const
{
	for(int iSub = 0; iSub < scomm->numT(SComm::std); ++iSub)
		pt->setLen(subNum(), scomm->neighbT(SComm::std,iSub), scomm->lenT(SComm::std,iSub)*_numRBM);
}

void
FetiBaseSub::setCommSize(FSCommStructure *pt, int size) const
{
	for(int iSub = 0; iSub < scomm->numNeighb; ++iSub)
		pt->setLen(subNum(), scomm->subNums[iSub], size);
}

void
FetiBaseSub::setDofCommSize(FSCommStructure *pt) const
{
	for(int iSub = 0; iSub < scomm->numT(SComm::all); ++iSub)
		pt->setLen(subNum(), scomm->neighbT(SComm::all,iSub), scomm->lenT(SComm::all,iSub));
}

void
FetiBaseSub::setMpcNeighbCommSize(FSCommPattern<int> *pt, int size) const
{
	for(int iSub = 0; iSub < scomm->numT(SComm::mpc); ++iSub)
		pt->setLen(subNum(), scomm->neighbT(SComm::mpc,iSub), size);
}

void
FetiBaseSub::computeMasterFlag(const Connectivity &mpcToSub)
{
	// PJSA: 12-13-02  masterFlag to be used in dot product and orthogonalization
	// allows for mpcs or wet interface dofs connecting 1 or > 2 subs
	masterFlag.resize(totalInterfSize);
	int  i, j;

	std::vector<bool> mpcFlag(numMPC, true);

	std::vector<bool> wiFlag(numWIdof, true);

	if(numWIdof && wiMaster.size() == 0) { // build wiMaster
		wiMaster.resize(numWIdof);  // wiMaster[i] is true if this subdomain is the master of the wet interface dof i
		for(i=0; i<numWIdof; ++i) wiMaster[i] = true;
		for(i=0; i < scomm->numT(SComm::wet); ++i) {
			if(scomm->neighbT(SComm::wet, i) < subNum())
				// set wiMaster false if this isn't the lowest numbered subdomain sharing the wet interface dof
				for(j=0; j < scomm->lenT(SComm::wet, i); ++j)
					wiMaster[scomm->wetDofNb(i, j)] = false;
		}
	}

	if(numMPC && !mpcMaster) { // PJSA moved here from SubDomain::scatterHalfInterf
		mpcMaster = new bool[numMPC];  // only allocate & init 1st time, dual mpcs only
		for(i=0; i<numMPC; ++i) mpcMaster[i] = false;
	}
//	std::cout << "Step 4 " << numMPC << " " << numWIdof << " " << masterFlag.size() << " " << boundDofFlag.size() << std::endl;

	int nbdofs = 0;
	masterFlagCount = 0;
	for(int iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
		int rank = (scomm->neighbT(SComm::all,iSub) < subNum()) ? 0 : 1;
		int count = 0;
		for(j=0; j<scomm->lenT(SComm::all,iSub); ++j) {
			// The following should never happen. It is only there to get around a
			// strange bug of clang++ on Mac OS...
			if(boundDofFlag[nbdofs] > 2)
				std::cout << "Strange: " << boundDofFlag[nbdofs] << std::endl;
			switch(boundDofFlag[nbdofs]) {
				case 0: {
					if((count % 2) == rank) {
						masterFlag[nbdofs++] = true;
						masterFlagCount++;
					}
					else masterFlag[nbdofs++] = false;
					count++;
				} break;
				case 1: { // wet interface
					int bdof = scomm->boundDofT(SComm::all,iSub,j);
					int windex = -1 - bdof;
					if(wiMaster[windex]) {
						if(wiFlag[windex]) { // only need one master for each WI dof
							masterFlag[nbdofs++] = true;
							masterFlagCount++;
							wiFlag[windex] = false;
						}
						else masterFlag[nbdofs++] = false;
					}
					else masterFlag[nbdofs++] = false;
				} break;
				case 2: { // mpc
					int bdof = scomm->boundDofT(SComm::all,iSub,j);
					int locMpcNb = -1 - bdof;
					int glMpcNb = localToGlobalMPC[locMpcNb];
					if(subNum() == mpcToSub[glMpcNb][0]) {
						mpcMaster[locMpcNb] = true; // PJSA
						if(mpcFlag[locMpcNb]) { // only need one master for each MPC
							masterFlag[nbdofs++] = true;
							masterFlagCount++;
							mpcFlag[locMpcNb] = false;
						}
						else masterFlag[nbdofs++] = false;
					}
					else masterFlag[nbdofs++] = false;
				} break;
			}
		}
	}
//	std::cout << "Step 5" << std::endl;

}

int
FetiBaseSub::numCoarseDofs()
{
	if(nCDofs == -1) {
		nCDofs = numCornerDofs();
		if(getFetiInfo().augment == FetiInfo::Gs) {
			nCDofs += nGrbm;
			for(int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
				nCDofs += neighbNumGRBMs[iSub];
			}
		}

		if(getFetiInfo().isEdgeAugmentationOn()) {
			for(int iSub = 0; iSub < scomm->numNeighb; ++iSub)
				nCDofs += edgeDofSize[iSub];
		}
		nCDofs += numMPC_primal; // MPC MODIFICATION: add the number of mpc equations
	}
	return nCDofs;
}

void
FetiBaseSub::setMpcCommSize(FSCommStructure *mpcPat) const
{
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc,i);
		int len = (subNum() != neighb) ? scomm->lenT(SComm::mpc,i) : 0;
		mpcPat->setLen(subNum(), neighb, len);
	}
}

void
FetiBaseSub::computeInternalMasterFlag()
{
	const int dofCount = get_c_dsa()->size();
	internalMasterFlag = new bool[dofCount];
	std::fill_n(internalMasterFlag, dofCount, true);

	for(int i = 0; i < scomm->numNeighb; ++i) {
		if(subNum() > scomm->subNums[i]) {
			for(int j = 0; j < scomm->sharedDOFsPlus->num(i); ++j) {
				internalMasterFlag[(*scomm->sharedDOFsPlus)[i][j]] = false;
			}
		}
	}
}

void
FetiBaseSub::sendMatProps(FSCommPattern<double> *matPat)
{
	for(int i = 0; i < scomm->numNeighb; ++i) {
		FSSubRecInfo<double> sInfo = matPat->getSendBuffer(subNum(), scomm->subNums[i]);
		sInfo.data[0] = Ymod;
		sInfo.data[1] = Prat;
		sInfo.data[2] = Dens;
		sInfo.data[3] = Thih;
		sInfo.data[4] = Sspe;
	}
}

void
FetiBaseSub::collectMatProps(FSCommPattern<double> *matPat)
{
	if(!neighbYmod) neighbYmod = new double[scomm->numNeighb];
	if(!neighbPrat) neighbPrat = new double[scomm->numNeighb];
	if(!neighbDens) neighbDens = new double[scomm->numNeighb];
	if(!neighbThih) neighbThih = new double[scomm->numNeighb];
	if(!neighbSspe) neighbSspe = new double[scomm->numNeighb];
	for(int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<double> rInfo = matPat->recData(scomm->subNums[iSub], subNum());
		neighbYmod[iSub] = (Ymod + rInfo.data[0])/2.0;
		neighbPrat[iSub] = (Prat + rInfo.data[1])/2.0;
		neighbDens[iSub] = (Dens + rInfo.data[2])/2.0;
		neighbThih[iSub] = (Thih + rInfo.data[3])/2.0;
		neighbSspe[iSub] = (Sspe + rInfo.data[4])/2.0;
	}
}


void
FetiBaseSub::sendWaveNumbers(FSCommPattern<double> *kPat)
{
  for(int i = 0; i < scomm->numNeighb; ++i) {
    FSSubRecInfo<double> sInfo = kPat->getSendBuffer(subNum(), scomm->subNums[i]);
    sInfo.data[0] = k_p;
    sInfo.data[1] = k_s;
    sInfo.data[2] = k_s2;
    sInfo.data[3] = k_f;
  }
}

void
FetiBaseSub::collectWaveNumbers(FSCommPattern<double> *kPat)
{
  if(!neighbK_p) neighbK_p = new double[scomm->numNeighb];
  if(!neighbK_s) neighbK_s = new double[scomm->numNeighb];
  if(!neighbK_s2) neighbK_s2 = new double[scomm->numNeighb];
  if(!neighbK_f) neighbK_f = new double[scomm->numNeighb];
  for(int i = 0; i < scomm->numNeighb; ++i) {
    FSSubRecInfo<double> rInfo = kPat->recData(scomm->subNums[i], subNum());
    neighbK_p[i] = (k_p + rInfo.data[0])/2.0;
    neighbK_s[i] = (k_s + rInfo.data[1])/2.0;
    neighbK_s2[i] = (k_s2 + rInfo.data[2])/2.0;
    neighbK_f[i] = (k_f + rInfo.data[3])/2.0;
  }
}

void
FetiBaseSub::findEdgeNeighbors()
{
	int count = 0;
	bool *isEdgeNeighb = new bool[scomm->numNeighb];  // deleted in ~SComm()
	for(int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		isEdgeNeighb[iSub] = false;
		for(int j=0; j<scomm->sharedNodes->num(iSub); ++j) {
			if(boundaryDOFs[iSub][j].count() > 0) {
				isEdgeNeighb[iSub] = true;
				count++;
				break;
			}
		}
	}
	scomm->setEdgeNeighb(count, isEdgeNeighb);
}

void
FetiBaseSub::zeroEdgeDofSize()
{
	if(edgeDofSize.size() != 0) {
		for(int i=0; i<scomm->numNeighb; ++i) edgeDofSize[i] = 0;
		nCDofs = -1;
	}
}

std::vector<int>
FetiBaseSub::makeBMaps(const DofSetArray *dof_set_array)
{
// Build the variables boundLen, boundMap[], dualToBoundary[]
// dof_set_array is c_dsa for FETI-DPC or cc_dsa for FETI-DP

	int iDof;

	int lLen = (dof_set_array==0) ? localLen() : dof_set_array->size();
	invBoundMap.resize(lLen);

	boundLen = 0;
	for(iDof = 0; iDof < lLen; ++iDof)
		invBoundMap[iDof] = (dofWeight(iDof) > 1) ? boundLen++ : -1 ;

	boundMap.resize(boundLen);
	boundLen = 0;
	for(iDof = 0; iDof < lLen; ++iDof)
		if(dofWeight(iDof) > 1) boundMap[boundLen++] = iDof;

	int gLen = getDsa()->size();
	std::vector<int> glBoundMap(gLen);
	for(iDof = 0; iDof < gLen; ++iDof) {
		int dofI = dof_set_array->getRCN(iDof);
		glBoundMap[iDof] = (dofI) < 0 ? -1 : invBoundMap[dofI];
	}

	if(dualToBoundary) delete [] dualToBoundary;
	dualToBoundary = new int[totalInterfSize];

	if(numMPC > 0) {
		for(iDof = 0; iDof < totalInterfSize; ++iDof)
			if(allBoundDofs[iDof] < 0)  // no preconditioning on virtual nodes yet
				dualToBoundary[iDof] = -1;
			else {
				dualToBoundary[iDof] = invBoundMap[ccToC[allBoundDofs[iDof]]];
				if(dualToBoundary[iDof] < 0) {
					std::cerr << "Error in makeBMaps " << std::endl;
					std::cerr << "Tried to map " << allBoundDofs[iDof]
					          << " which in cdsa is " << ccToC[allBoundDofs[iDof]] << std::endl;
					exit(-1);
				}
			}
	}
	else {
		for(iDof = 0; iDof < totalInterfSize; ++iDof) {
			if(allBoundDofs[iDof] < 0)  // COUPLED_DPH QUESTION: no preconditioning on wet interface nodes yet
				dualToBoundary[iDof] = -1;
			else{
				dualToBoundary[iDof] = invBoundMap[allBoundDofs[iDof]];
			}
		}
		invBoundMap.clear();
	}

	return glBoundMap;
}


std::vector<int>
FetiBaseSub::makeIMaps(const DofSetArray *dof_set_array)
{
	// Build the variables internalLen, internalMap[]
	int iDof;
	int lLen = (dof_set_array==0) ? localLen() : dof_set_array->size();
	int gLen = getDsa()->size();
	std::vector<int> glInternalMap(gLen);

	internalLen = numWIdof;
	for(iDof = 0; iDof < lLen; ++iDof) if(dofWeight(iDof) == 1) internalLen++; //add ir dofs
	internalMap.resize(internalLen);
	if(numWIdof) wiInternalMap.resize(numWIdof);
	internalLen = 0;
	for(iDof = 0; iDof < gLen; ++iDof) {
		int dofI = dof_set_array->getRCN(iDof);
		if((dofI > -1) && (dofWeight(dofI) == 1)) { // regular ir dof
			internalMap[internalLen] = dofI;
			glInternalMap[iDof] = internalLen++;
		}
		else if((numWIdof) && (wetInterfaceMap[iDof] > -1)) { // wet interface dof
			internalMap[internalLen] = -wetInterfaceMap[iDof];
			wiInternalMap[wetInterfaceMap[iDof]] = internalLen;
			glInternalMap[iDof] = internalLen++;
		}
		else glInternalMap[iDof] = -1;
	}

	return glInternalMap;
}
