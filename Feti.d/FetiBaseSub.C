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
	if(masterFlag) delete [] masterFlag;
	masterFlag = new bool[totalInterfSize];
	int rank, iSub, i, j;

	bool *mpcFlag =  (bool *) dbg_alloca(sizeof(bool)*numMPC);
	for(i=0; i<numMPC; ++i) mpcFlag[i] = true;

	bool *wiFlag = (bool *) dbg_alloca(sizeof(bool)*numWIdof);
	for(i=0; i<numWIdof; ++i) { wiFlag[i] = true; }

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

	int nbdofs = 0;
	masterFlagCount = 0;
	for(iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
		if(scomm->neighbT(SComm::all,iSub) < subNum()) rank = 0; else rank = 1;
		int count = 0;
		for(j=0; j<scomm->lenT(SComm::all,iSub); ++j) {
			int bdof = scomm->boundDofT(SComm::all,iSub,j);
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