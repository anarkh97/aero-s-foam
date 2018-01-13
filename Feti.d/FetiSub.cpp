//
// Created by Michel Lesoinne on 12/6/17.
//
#include <complex>
#include <Eigen/Sparse>
#include <Utils.d/dofset.h>
#include "FetiSub.h"
#include <Driver.d/Mpc.h>
#include <Utils.d/SolverInfo.h>
#include <Utils.d/dbg_alloca.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/SparseSet.h>
#include <Math.d/CuCSparse.h>
#include <Solvers.d/Solver.h>
#include <Math.d/BLAS.h>

void
FetiBaseSub::markCornerDofs(int *glCornerDofs) const
{
	for(int i=0; i<numCRN; ++i)
		glCornerDofs[glCornerNodes[i]] |= cornerDofs[i].list();
}


void
FetiBaseSub::setSComm(SComm *sc)
{
	scomm = sc;
}


void
FetiBaseSub::addMPCsToGlobalZstar(FullM *globalZstar, int startRow, int startCol, int numCol)
{
	FullM *Zmpc = rigidBodyModesG->Zmpc;
	for(int i=0; i<numMPC_primal; ++i)
		for(int j=0; j<numCol; ++j)
			(*globalZstar)[startRow+localToGroupMPC[i]][startCol+j] += (*Zmpc)[i][j];
}

void
FetiBaseSub::addSPCsToGlobalZstar(FullM *globalZstar, int &zRow, int zColOffset)
{
	FullM *Zstar = rigidBodyModesG->Zstar;
	globalZstar->add(*Zstar, zRow, zColOffset);
	zRow += Zstar->numRow();
}


void
FetiBaseSub::setWIoneCommSize(FSCommStructure *pat) const
{
	for(int i = 0; i < scomm->numT(SComm::fsi); ++i)
		if(subNum() != scomm->neighbT(SComm::fsi,i))
			pat->setLen(subNum(), scomm->neighbT(SComm::fsi,i), 1);
}

void
FetiBaseSub::setWICommSize(FSCommStructure *pat) {
	for (int i = 0; i < scomm->numT(SComm::fsi); ++i)
		pat->setLen(subNum(), scomm->neighbT(SComm::fsi, i), numNeighbWIdof[i]);
}

void
FetiBaseSub::setWImapCommSize(FSCommPattern<int> *pat)
{
	for(int i = 0; i < scomm->numT(SComm::fsi); ++i)
		if(subNum() != scomm->neighbT(SComm::fsi,i))
			glToLocalWImap.setCommSize(pat, subNum(), scomm->neighbT(SComm::fsi,i));
}


void FetiBaseSub::setNumGroupRBM(int *ngrbmGr)
{
	groupRBMoffset = 0;
	for(int i=0; i<group; ++i) groupRBMoffset += ngrbmGr[i];
	numGroupRBM = ngrbmGr[group];
}

void FetiBaseSub::getNumGroupRBM(int *ngrbmGr)
{
	//cerr << "in getNumGroupRBM " << subNumber << " " << group << " " << numGroupRBM << std::endl;
	ngrbmGr[group] = numGroupRBM;
}

void
FetiBaseSub::sendNumWIdof(FSCommPattern<int> *pat) const
{
	for(int i = 0; i < scomm->numT(SComm::fsi); ++i) {
		if(subNum() != scomm->neighbT(SComm::fsi,i)) {
			FSSubRecInfo<int> sInfo = pat->getSendBuffer(subNum(), scomm->neighbT(SComm::fsi,i));
			sInfo.data[0] = numWIdof;
		}
	}
}

void
FetiBaseSub::recvNumWIdof(FSCommPattern<int> *pat)
{
	numNeighbWIdof.resize(scomm->numT(SComm::fsi));
	for(int i = 0; i < scomm->numT(SComm::fsi); ++i) {
		if(subNum() != scomm->neighbT(SComm::fsi,i)) {
			FSSubRecInfo<int> rInfo = pat->recData(scomm->neighbT(SComm::fsi,i), subNum());
			numNeighbWIdof[i] = rInfo.data[0];
		}
		else numNeighbWIdof[i] = 0;
	}
}

void
FetiBaseSub::sendWImap(FSCommPattern<int> *pat)
{
	for(int i=0; i< scomm->numT(SComm::fsi); ++i)
		if(subNum() != scomm->neighbT(SComm::fsi,i))
			glToLocalWImap.pack(pat,subNum(), scomm->neighbT(SComm::fsi,i));
}

void
FetiBaseSub::recvWImap(FSCommPattern<int> *pat)
{
	if(neighbGlToLocalWImap) delete [] neighbGlToLocalWImap;
	neighbGlToLocalWImap = new GlobalToLocalMap[scomm->numT(SComm::fsi)];
	for(int i=0; i<scomm->numT(SComm::fsi); ++i)
		if(subNum() != scomm->neighbT(SComm::fsi,i))
			neighbGlToLocalWImap[i].unpack(pat, scomm->neighbT(SComm::fsi,i), subNum());
}


void
FetiBaseSub::sendNeighbGrbmInfo(FSCommPattern<int> *pat)
{
	// send number of group GRBMs and the group GRBM offset to each potential contact neighbor
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		FSSubRecInfo<int> sInfo = pat->getSendBuffer(subNum(), neighb);
		sInfo.data[0] = numGroupRBM;
		sInfo.data[1] = groupRBMoffset;
	}
}

void
FetiBaseSub::receiveNeighbGrbmInfo(FSCommPattern<int> *pat)
{
	if(neighbNumGroupGrbm) delete [] neighbNumGroupGrbm;
	neighbNumGroupGrbm = new int[scomm->numT(SComm::mpc)];
	if(neighbGroupGrbmOffset) delete [] neighbGroupGrbmOffset;
	neighbGroupGrbmOffset = new int[scomm->numT(SComm::mpc)];
	// get number of group GRBMs and the group GRBM offset for each potential contact neighbor
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		FSSubRecInfo<int> rInfo = pat->recData(neighb, subNum());
		neighbNumGroupGrbm[i] = rInfo.data[0];
		neighbGroupGrbmOffset[i] = rInfo.data[1];
	}
}


void
FetiBaseSub::sendNumNeighbGrbm(FSCommPattern<int> *pat)
{
	// send Number of RBMs for each neighbor, used for augmentation
	for(int i = 0; i < scomm->numT(SComm::std); ++i) {
		FSSubRecInfo<int> sInfo = pat->getSendBuffer(subNum(), scomm->neighbT(SComm::std,i));
		sInfo.data[0] = nGrbm;
	}
}

void
FetiBaseSub::recvNumNeighbGrbm(FSCommPattern<int> *pat)
{
	neighbNumGRBMs = new int[scomm->numT(SComm::std)];
	// get Number of RBMs for each neighbor, used for augmentation
	for(int i = 0; i < scomm->numT(SComm::std); ++i) {
		FSSubRecInfo<int> rInfo = pat->recData(scomm->neighbT(SComm::std,i), subNum());
		neighbNumGRBMs[i] = rInfo.data[0];
	}
}

void FetiBaseSub::makeLocalToGroupMPC(Connectivity *groupToMPC)
{
	// PJSA: new version for multi-body mpc compatability
	int i;
	if(numMPC_primal > 0) {
		localToGroupMPC = new int[numMPC_primal];
		int groupOffset = groupToMPC->offset(group);
		for(i=0; i<groupToMPC->num(group); ++i) {
			int glMpcID = groupToMPC->getTargetValue(groupOffset+i);
			int localMpcID = globalToLocalMPC_primal[glMpcID];
			if(localMpcID > -1) localToGroupMPC[localMpcID] = i;
		}
	}
}

template <typename Scalar>
void FetiSub<Scalar>::makeBs() {
	std::vector<Eigen::Triplet<double>> b;
	std::vector<Eigen::Triplet<double>> bw;
	std::vector<Eigen::Triplet<Scalar>> bm;
	std::vector<Eigen::Triplet<Scalar>> bc;
	std::vector<bool> mpcFlag(numMPC, true);
	for (int iDof = 0; iDof < totalInterfSize; ++iDof) {
		switch (boundDofFlag[iDof]) {
			case 0: // note B is used with a - sign.
				b.emplace_back(allBoundDofs[iDof], iDof, 1.0);
				break;
			case 1:  // wet interface Note Bw is used with a + sign
				bw.emplace_back(-1 - allBoundDofs[iDof], iDof, 1.0);
				break;
			case 2: { // dual mpc or contact. Note Bm is used with a - sign.
				int locMpcNb = -1 - allBoundDofs[iDof];
				if (mpcFlag[locMpcNb]) {
					const auto &m = mpc[locMpcNb];
					for (int k = 0; k < m->nterms; k++) {
						int ccdof = (m->terms)[k].ccdof;
						bm.emplace_back(ccdof, iDof, (m->terms)[k].coef);
					}
					mpcFlag[locMpcNb] = false;
				}
			}
				break;
		}
	}
	B.resize(localLen(), totalInterfSize);
	B.setFromTriplets(b.begin(), b.end());
	Bw.resize(numWIdof, totalInterfSize);
	Bw.setFromTriplets(bw.begin(), bw.end());
	Bm.resize(localLen(), totalInterfSize);
	Bm.setFromTriplets(bm.begin(), bm.end());

	mpcFlag.assign(numMPC, true);
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpcFlag[locMpcNb]) {
			const auto &m = mpc[locMpcNb];
			for (int k = 0; k < m->nterms; k++) {
				int dof = (m->terms)[k].dof;
				if ((dof >= 0) && (cornerMap[dof] >= 0))
					bc.emplace_back(cornerMap[dof], scomm->mapT(SComm::mpc, i), (m->terms)[k].coef);
			}
			mpcFlag[locMpcNb] = false;
		}
	}
	Bc.resize(Src->numCol(), totalInterfSize);
	Bc.setFromTriplets(bc.begin(), bc.end());
}

template<typename Scalar>
double FetiSub<Scalar>::getMpcError() const {
	double ret = 0;
	for(int i = 0; i < numMPC; ++i) {
		if(mpcMaster[i]) {
			if(mpc[i]->type == 0) {
				ret += ScalarTypes::sqNorm(mpc[i]->rhs);
			}
			else if(mpc[i]->type == 1 && ScalarTypes::lessThan(mpc[i]->rhs, 0.)) ret += ScalarTypes::sqNorm(mpc[i]->rhs);
		}
	}
	return ret;
}

template<class Scalar>
void
FetiSub<Scalar>::makeLocalMpcToDof() {
	auto &mpc = this->mpc;
	if (mpcToDof) delete mpcToDof;
	mpcToDof = 0;

	// step 1: make mpcToDof
	int size = numMPC;
	// step 1.1: find size of target: total number of coefficients involving a different dof
	int numtarget = 0;
	int i, j, jj;
	for (i = 0; i < size; i++) {
		for (j = 0; j < mpc[i]->nterms; j++) {
			int dofj = mpc[i]->terms[j].dof;
			for (jj = 0; jj < j; jj++) {
				int dofjj = mpc[i]->terms[jj].dof;
				if (dofj == dofjj) break;
			}
			if ((jj == j) && (dofj >= 0)) numtarget++;
		}
	}
	// step 1.2: fill target with coefficient dofs
	int *pointer = new int[size + 1];
	int *target = new int[numtarget];
	int count = 0;
	for (i = 0; i < size; i++) {
		pointer[i] = count;
		for (j = 0; j < mpc[i]->nterms; j++) {
			int dofj = mpc[i]->terms[j].dof;
			for (jj = 0; jj < j; jj++) {
				int dofjj = mpc[i]->terms[jj].dof;
				if (dofj == dofjj) break;
			}
			if ((jj == j) && (dofj >= 0)) {
				target[count] = dofj;
				count++;
			}
		}
	}
	pointer[i] = numtarget;
	// step 1.3: construct mpcToDof connectivity
	mpcToDof = new Connectivity(size, pointer, target);
}

template<class Scalar>
void
FetiSub<Scalar>::makeLocalMpcToMpc() {
	// step 1: make mpcToDof
	makeLocalMpcToDof();

	// step 2: make localMpcToMpc connectivity
	Connectivity *dofToMpc = mpcToDof->reverse();
	localMpcToMpc = mpcToDof->transcon(dofToMpc);
	delete dofToMpc;
}


template<class Scalar>
void
FetiSub<Scalar>::applyMpcSplitting()
{
	// adjust discrete masses, forces and mpcs using subdomain multiplicity
	// num = number of subdomains touching a dof
	int cdof, num;

	// mpcs (NOTE: optional kscaling is done later, hhs is not split)
	if(getFetiInfo().mpc_scaling == FetiInfo::tscaling) {
		for(int iMPC = 0; iMPC < numMPC; ++iMPC) { // dual mpcs
			if(mpc[iMPC]->type == 2) continue; // bmpc
			for(int i = 0; i < mpc[iMPC]->nterms; ++i) {
				if((cdof = mpc[iMPC]->terms[i].cdof) > -1 && (num = weightPlus[cdof]) > 1)
					mpc[iMPC]->terms[i].coef /= double(num);
			}
		}
	}
	// XXXX kscaling currently not supported for primal mpcs
	for(int iMPC = 0; iMPC < numMPC_primal; ++iMPC) { // primal mpcs
		for(int i = 0; i < mpc_primal[iMPC]->nterms; ++i) {
			if((cdof = mpc_primal[iMPC]->terms[i].cdof) > -1 && (num = weightPlus[cdof]) > 1)
				mpc_primal[iMPC]->terms[i].coef /= double(num);
		}
	}

}


template<class Scalar>
void
FetiSub<Scalar>::subtractMpcRhs(Scalar *interfvec)
{
	for(int i=0; i < scomm->lenT(SComm::mpc); ++i) {
		interfvec[scomm->mapT(SComm::mpc,i)] -= mpc[scomm->mpcNb(i)]->rhs;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::assembleGtGsolver(GenSparseMatrix<Scalar> *GtGsolver)
{
	if(numGroupRBM == 0) return;
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool)*numMPC);
	for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int numGroupRBM2 = neighbNumGroupGrbm[i];
		int groupRBMoffset2 = neighbGroupGrbmOffset[i];
		GenVector<Scalar> d(scomm->lenT(SComm::mpc,i));
		for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
			int locMpcNb = scomm->mpcNb(i,j);
			d[j] = (mpc[locMpcNb]->active) ? 0.0 : 1.0;
		}
		if((numGroupRBM2 > 0) && (subNum() != scomm->neighbT(SComm::mpc,i))) {
			GenFullM<Scalar> tmp2(numGroupRBM, numGroupRBM2);  // coupling term
			G[i]->transposeMultD(*neighbG[i], d, tmp2); // tmp2 = G^T * D * neighbG
			GtGsolver->add(tmp2, groupRBMoffset, groupRBMoffset2);
		}
		for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
			int locMpcNb = scomm->mpcNb(i,j);
			if(!mpcFlag[locMpcNb]) d[j] = 0.0; // prevents duplication for mpc shared between more than 2 subdomains
			else mpcFlag[locMpcNb] = false;
		}
		GenFullM<Scalar> tmp(numGroupRBM, numGroupRBM);
		G[i]->transposeMultD(*G[i], d, tmp); // tmp = G^T * D * G
		GtGsolver->add(tmp, groupRBMoffset, groupRBMoffset);
	}
}


template<class Scalar>
void
FetiSub<Scalar>::getLocalMpcForces(double *mpcLambda, DofSetArray *cornerEqs,
                                        int mpcOffset, GenVector<Scalar> &uc) {
// XXXX needs some work to map both dual and primal into single mpcLambda array
	if (numMPC > 0 && numMPC_primal > 0) std::cerr << "unsupported feature in FetiSub::getLocalMpcForces \n";
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) { // dual mpcs
		int locMpcNb = scomm->mpcNb(i);
		if (localLambda) mpcLambda[locMpcNb] = localLambda[scomm->mapT(SComm::mpc, i)];
		else mpcLambda[locMpcNb] = 0;
	}
	for (int i = 0; i < numMPC_primal; ++i) {
		int glMpcNb = localToGlobalMPC_primal[i];
		int dof = cornerEqs->firstdof(mpcOffset + glMpcNb);
		mpcLambda[i] = -ScalarTypes::Real(uc[dof]);
	}
//	if (salinasFlag)
//		for (int i = 0; i < numMPC + numMPC_primal; ++i)
//			mpcLambda[i] = -mpcLambda[i];  // different sign convention
}


template<class Scalar>
void
FetiSub<Scalar>::useKrrNullspace()
{
	// EXPERMENTAL... use alaebraic null space of Krr with no corners
	int neq = this->Krr->neqs();
	int nzem = this->Krr->numRBM();
	Rstar.setNewSize(neq, nzem);
	if(nzem > 0) {
		std::vector<Scalar> rbmv(neq*nzem);
		this->Krr->getNullSpace(rbmv.data());

		// Copy rigid body modes (rbm) into Rstar
		for(int m=0; m<nzem; ++m) {
			for(int i=0; i<neq; ++i)
				Rstar[i][m] = rbmv[i+m*neq];
		}
	}
}

template<>
void
FetiSub<double>::makeLocalRstar(FullM **Qtranspose)
{
	FullM &R = rigidBodyModesG->R;
	FullM *Rc = rigidBodyModesG->Rc;
	if(numMPC_primal > 0) {
		FullM Qbody(Qtranspose[group]->transpose(), R.numCol(), bodyRBMoffset, Qtranspose[group]->numRow(), 0);
		Rstar = R * Qbody;
	}
	else {
		Rstar = R % *(Qtranspose[group]);
	}
}

template<>
void
FetiSub<DComplex>::makeLocalRstar(FullM **Qtranspose)
{
	std::cerr << "FetiSub<DComplex>::makeLocalRstar(FullM **Qtranspose) is not implemented\n";
}

template<class Scalar>
void
FetiSub<Scalar>::addRalpha(Scalar *u, GenVector<Scalar> &alpha) const
{
	int i, j;
	for(i=0; i<Rstar.numRow(); ++i)
		for(j=0; j<Rstar.numCol(); ++j)
			u[i] += Rstar[i][j] * alpha[groupRBMoffset + j];
}

template<class Scalar>
void
FetiSub<Scalar>::assembleE(GenVector<Scalar> &e, Scalar *f) const
{
	if(numGroupRBM > 0) {
		GenVector<Scalar> local_e(numGroupRBM, 0.0);
		GenVector<Scalar> fvec(f, Rstar.numRow());
		local_e = Rstar ^ fvec; // = Rtranspose * fvec
		e.add(local_e, groupRBMoffset);
	}
}

template<class Scalar>
void
FetiSub<Scalar>::makeG()
{
	// make G for each potential contact/mpc neighbour
	G.resize(scomm->numT(SComm::mpc));
	neighbG.resize(scomm->numT(SComm::mpc));
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {  // loop through all potential contact/mpc neighbours
		neighbG[i] = 0;
		G[i] = std::make_unique<GenFullM<Scalar>>(scomm->lenT(SComm::mpc,i), numGroupRBM);
		if(numGroupRBM > 0) {
			G[i]->zero();
			for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {  // loop through potential contact/mpc nodes
				int locMpcNb = scomm->mpcNb(i,j);
				const auto &m = mpc[locMpcNb];
				for(int k = 0; k < m->nterms; k++) {
					int cDof = (m->terms)[k].cdof;
					if(cDof > -1) {
						for(int iRbm = 0; iRbm < numGroupRBM; ++iRbm)
							(*G[i])[j][iRbm] += Rstar[cDof][iRbm]*(m->terms)[k].coef;
					}
				}
			}
		}
	}
}


template<class Scalar>
void
FetiSub<Scalar>::makeTrbmG(Scalar *rbms, int nrbm, int size)
{
	auto &Src = this->Src;
	auto &mpc = this->mpc;
	// rbms is the null space of the global Kcc^* matrix
	// nrbm is the nullity of the global Kcc^* matrix
	// size is the number of rows and columns of the global Kcc^* matrix
	// TODO what about augmentation
	int numCDofs = (Src) ? Src->numCol() : 0;
	Scalar *localc = (Scalar *) dbg_alloca(sizeof(Scalar)*numCDofs);

	int nrbms_local = 0; int first = 0;
	std::map<int, int> localToGlobalRBM;
	std::map<int, int> globalToLocalRBM;
	for(int iRbm = 0; iRbm < nrbm; ++iRbm) {
		Scalar *rbm = rbms + size*iRbm;
		Scalar dot = 0.0;
		for(int i=0; i<numCDofs; ++i) {
			if(cornerEqNums[i] > -1) { dot += rbm[cornerEqNums[i]]*rbm[cornerEqNums[i]]; }
		}
		if(dot != 0.0) {
			localToGlobalRBM[nrbms_local] = iRbm;
			globalToLocalRBM[iRbm] = nrbms_local;
			if(nrbms_local == 0) first = iRbm;
			nrbms_local++;
		}
	}

	numGroupRBM = nrbms_local; // TODO this isn't general since global trbms may not be grouped like grbms
	groupRBMoffset = first;    // TODO this isn't general

	G.resize(scomm->numT(SComm::mpc));
	neighbG.resize(scomm->numT(SComm::mpc));
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		neighbG[i] = 0;
		G[i] = std::make_unique<GenFullM<Scalar>>(scomm->lenT(SComm::mpc,i), nrbms_local);
		G[i]->zero();
	}

	if(nrbms_local == 0 || numMPC == 0) return;

	Scalar *localr = new Scalar[localLen()];

	for(int iRbm = 0; iRbm < nrbms_local; ++iRbm) {
		int glRbmId = localToGlobalRBM[iRbm];
		Scalar *rbm = rbms + size*glRbmId;
		for(int i=0; i<numCDofs; ++i) {
			if(cornerEqNums[i] > -1) { localc[i] = rbm[cornerEqNums[i]]; }
			else localc[i] = 0.0;
		}

		// G = (-Br^(s)Krr^{-1}Krc + Bc)Lcc Nc
		for(int i=0; i<localLen(); ++i) localr[i] = 0.0;
		if(Src) Src->transposeMultSubtract(localc, localr);
		if(this->Krr) this->Krr->reSolve(localr);
		for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
			for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
				int locMpcNb = scomm->mpcNb(i,j);
				const auto &m = mpc[locMpcNb];
				for(int k = 0; k < m->nterms; k++) {
					int cc_dof = (m->terms)[k].ccdof;
					if(cc_dof >= 0) (*G[i])[j][iRbm] += localr[cc_dof]*(m->terms)[k].coef;
					else {
						int dof = (m->terms)[k].dof;
						if((dof >= 0) && (cornerMap[dof] >= 0))
							(*G[i])[j][iRbm] += localc[cornerMap[dof]] * (m->terms)[k].coef;
					}
				}
			}
		}
	}
	delete [] localr;
}

// ****************************************************************************************************

template<class Scalar>
void
FetiSub<Scalar>::setGCommSize(FSCommStructure *pat) const
{
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int nRow = G[i]->numRow();
		int nCol = numGroupRBM;
		pat->setLen(subNum(), scomm->neighbT(SComm::mpc, i), nRow*nCol);
	}
}

template<class Scalar>
void
FetiSub<Scalar>::sendG(FSCommPattern<Scalar> *rbmPat)
{
	if(numGroupRBM == 0) return;
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i)
		rbmPat->sendData(subNum(), scomm->neighbT(SComm::mpc, i),  G[i]->data());
}

template<class Scalar>
void
FetiSub<Scalar>::receiveG(FSCommPattern<Scalar> *rbmPat)
{
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		FSSubRecInfo<Scalar> rInfo = rbmPat->recData(scomm->neighbT(SComm::mpc, i), subNum());
		int nRow = G[i]->numRow();  // number of potential contact dofs on interface with neighb
		int nCol = neighbNumGroupGrbm[i];  //number of rbms for neighb's group
		neighbG[i] = std::make_unique<GenFullM<Scalar>>(rInfo.data, nRow, nCol);
	}
}

template<class Scalar>
void FetiSub<Scalar>::zeroG()
{
	if(G.size() != 0) {
		for(int i=0; i<scomm->numT(SComm::mpc); ++i)
			if(G[i]) G[i]->zero();
	}
	if(neighbG.size() != 0) {
		for(int i=0; i<scomm->numT(SComm::mpc); ++i)
			if(neighbG[i]) neighbG[i]->zero();
	}
}

template<class Scalar>
void FetiSub<Scalar>::deleteG()
{
	G.clear();
	neighbG.clear();
	if(neighbNumGroupGrbm) { delete [] neighbNumGroupGrbm; neighbNumGroupGrbm = 0; }
	if(neighbGroupGrbmOffset) { delete [] neighbGroupGrbmOffset; neighbGroupGrbmOffset = 0; }
}

template<class Scalar>
void
FetiSub<Scalar>::multG(const GenVector<Scalar> &x, Scalar *y, Scalar alpha) const
{
	// y += alpha * G * x
	Scalar *mpcvec = new Scalar[numMPC];
	for(int i = 0; i < numMPC; ++i) mpcvec[i] = 0.0;
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		for(int j = 0; j < scomm->lenT(SComm::mpc, i); ++j) {
			int locMpcNb = scomm->mpcNb(i,j);
			if(mpcvec[locMpcNb] == 0.0)
				for(int k = 0; k < numGroupRBM; ++k)
					mpcvec[locMpcNb] += (*G[i])[j][k] * x[k + groupRBMoffset];
			if(subNum() != neighb)
				for(int k = 0; k < neighbNumGroupGrbm[i]; ++k)
					mpcvec[locMpcNb] += (*neighbG[i])[j][k] * x[k + neighbGroupGrbmOffset[i]];
		}
	}
	for(int i = 0; i < scomm->lenT(SComm::mpc); ++i)
		y[scomm->mapT(SComm::mpc,i)] += alpha*mpcvec[scomm->boundDofT(SComm::mpc,i)];
	delete [] mpcvec;
}

template<class Scalar>
void
FetiSub<Scalar>::trMultG(const Scalar *x, GenVector<Scalar> &y, Scalar alpha) const
{
	// compute y += alpha * G^t * x
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool)*numMPC);
	for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
			int locMpcNb = scomm->mpcNb(i,j);
			int iDof = scomm->mapT(SComm::mpc,i,j);
			if(mpcFlag[locMpcNb]) {
				for(int k = 0; k < numGroupRBM; ++k)
					y[k + groupRBMoffset] += alpha * (*G[i])[j][k] * x[iDof];
				mpcFlag[locMpcNb] = false;
			}
		}
	}
}


// ************************************************************************************************

template<class Scalar>
void
FetiSub<Scalar>::buildGlobalRBMs(GenFullM<Scalar> &Xmatrix, const Connectivity *cornerToSub)
{
	int i,j,k;
	if(numGroupRBM == 0) numGlobalRBMs = 0;
	else {
		numGlobalRBMs = Xmatrix.numCol();
		GenFullM<Scalar> groupX(Xmatrix, numGroupRBM, groupRBMoffset, numGlobalRBMs, 0);
		Rstar_g = Rstar * groupX;

		if(!cornerToSub) return;
		std::vector<double> sharedUse(Rstar.numRow(), 1.0);
		// if i is a shared dof set sharedUse[i] = 0.0
		// for all but one of the subdomains sharing it in this body

		for(i=0; i<scomm->numT(SComm::std); ++i) {  // check non-corner dofs
			if(subNum() > scomm->neighbT(SComm::std,i))
				for(j=0; j<scomm->lenT(SComm::std,i); ++j)
					sharedUse[ccToC[scomm->boundDofT(SComm::std,i,j)]] = 0.0;
		}
		for(i=0; i<numCRN; ++i) { // check corner dofs
			if(subNum() != (*cornerToSub)[glCornerNodes[i]][0]) {
				int lDof[6];
				get_c_dsa()->number(cornerNodes[i], cornerDofs[i].list(), lDof);
				for(j=0; j<cornerDofs[i].count(); ++j)
					if(lDof[j] >= 0) sharedUse[lDof[j]] = 0.0;
			}
		}

		// if i is a shared dof sharedRstar_g[i] is set to zero 
		// for all but one of the subdomains sharing it in this body
		// (used to prevent duplication in construction of RtR)
		sharedRstar_g = std::make_unique<GenFullM<Scalar>>(Rstar_g.numRow(), Rstar_g.numCol());
		for(i=0; i<Rstar_g.numRow(); ++i)
			for(j=0; j<Rstar_g.numCol(); ++j) (*sharedRstar_g)[i][j] = Rstar_g[i][j] * sharedUse[i];

		// if i is a shared dof tmpRstar_g[i] is set to Rstar_g[i]/n 
		// where n is the number of subdomains (is this body) sharing this dof
		// used to compute a distributed vector, since distvec[i]*n = actual value at dof i  
		tmpRstar_g = std::make_unique<GenFullM<Scalar>>(Rstar_g);
		for(i=0; i<scomm->numT(SComm::all); ++i) {
			for(j=0; j<scomm->lenT(SComm::all,i); ++j) {
				int bdof = scomm->boundDofT(SComm::all,i,j);
				if(bdof >= 0)  // not a contact dof 
					for(k=0; k<numGlobalRBMs; ++k)
						(*tmpRstar_g)[ccToC[bdof]][k] /= getWeights()[bdof];
			}
		}
		for(i=0; i<numCRN; ++i) {
			int lDof[6];
			get_c_dsa()->number(cornerNodes[i], cornerDofs[i].list(), lDof);
			for(j=0; j<cornerDofs[i].count(); ++j)
				if(lDof[j] >= 0)
					for(k=0; k<numGlobalRBMs; ++k)
						(*tmpRstar_g)[lDof[j]][k] /= (double) (cornerToSub->num(glCornerNodes[i]));
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::getGlobalRBM(int iRBM, Scalar *Rvec) const
{
	if(numGlobalRBMs > 0)
		for(int iRow=0; iRow<Rstar_g.numRow(); ++iRow)
			Rvec[iRow] = Rstar_g[iRow][iRBM];
}

template<class Scalar>
void
FetiSub<Scalar>::subtractRstar_g(Scalar *u, GenVector<Scalar> &beta) const
{
	int i;
	if(numGlobalRBMs > 0) {
		// compute u = u - Rstar_g * beta  (second part of displacement projection)
		GenVector<Scalar> tmpu(Rstar_g.numRow());
		tmpu = Rstar_g * beta;
		for(i=0; i<Rstar_g.numRow(); ++i)
			u[i] = u[i] - tmpu[i];
	}
}

template<class Scalar>
void
FetiSub<Scalar>::addRstar_gT(Scalar *u, GenVector<Scalar> &beta) const
{
	if(numGlobalRBMs > 0) {
		// compute beta += Rstar_g^t * u  (first part of displacement projection)
		GenVector<Scalar> uvec(u, Rstar_g.numRow());
		beta += *tmpRstar_g ^ uvec;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::assembleRtR(GenFullM<Scalar> &RtR)
{
	// builds RtR for u projection
	if(numGlobalRBMs > 0) {
		GenFullM<Scalar> tmp(numGlobalRBMs, numGlobalRBMs);
		sharedRstar_g->transposeMult((*sharedRstar_g), tmp);
		RtR.add(tmp, 0, 0);
	}
}


template<class Scalar>
void
FetiSub<Scalar>::multAddCT(const Scalar *interfvec, Scalar *localvec) const {
	auto &mpc = this->mpc;
	// localvec += C^T * interfvec
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
	for (int i = 0; i < numMPC; ++i) mpcFlag[i] = true;

	for (int i = 0; i < scomm->lenT(SComm::mpc); i++) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpcFlag[locMpcNb]) {
			const auto &m = mpc[locMpcNb];
			int iDof = scomm->mapT(SComm::mpc, i);
			for (int k = 0; k < m->nterms; k++) {
				int c_dof = (m->terms)[k].cdof;
				if (c_dof >= 0) localvec[c_dof] += interfvec[iDof] * (m->terms)[k].coef;
			}
			mpcFlag[locMpcNb] = false;
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::multC(const Scalar *localvec, Scalar *interfvec) const {
	auto &mpc = this->mpc;
	// interfvec = C * localvec
	for (int i = 0; i < scomm->lenT(SComm::mpc); i++) {
		int locMpcNb = scomm->mpcNb(i);
		const auto &m = mpc[locMpcNb];
		int iDof = scomm->mapT(SComm::mpc, i);
		interfvec[iDof] = 0;
		for (int k = 0; k < m->nterms; k++) {
			int c_dof = (m->terms)[k].cdof;
			if (c_dof >= 0) interfvec[iDof] += localvec[c_dof] * (m->terms)[k].coef;
		}
	}
}


template<class Scalar>
void
FetiSub<Scalar>::addTrbmRalpha(Scalar *rbms, int nrbms, int glNumCDofs, Scalar *alpha, Scalar *ur) const
{
	auto &Src = this->Src;
	int numCDofs = (Src) ? Src->numCol() : 0;
	Scalar *localc = (Scalar *) dbg_alloca(sizeof(Scalar)*numCDofs);
	Scalar *localr = new Scalar[localLen()];

	for(int iRbm = 0; iRbm < nrbms; ++iRbm) {
		Scalar *rbm = rbms + glNumCDofs*iRbm;
		for(int i=0; i<numCDofs; ++i) {
			if(cornerEqNums[i] > -1) { localc[i] = rbm[cornerEqNums[i]]*alpha[iRbm]; }
			else localc[i] = 0.0;
		}

		for(int i = 0; i < localLen(); ++i) localr[i] = 0.0;
		if(Src) Src->transposeMultAdd(localc, localr);
		if(this->Krr) this->Krr->reSolve(localr);

		for(int i = 0; i < localLen(); ++i) ur[i] -= localr[i];
	}

	delete [] localr;
}

template<class Scalar>
void
FetiSub<Scalar>::assembleTrbmE(Scalar *rbms, int nrbms, int size, Scalar *e, Scalar *fr) const
{
	auto &Src = this->Src;
	int numCDofs = (Src) ? Src->numCol() : 0;
	Scalar *localc = (Scalar *) dbg_alloca(sizeof(Scalar)*numCDofs);
	for(int i=0; i<numCDofs; ++i) localc[i] = 0.0;

	Scalar *localr = new Scalar[localLen()];
	for(int i = 0; i < localLen(); ++i) localr[i] = -fr[i];
	if(this->Krr) this->Krr->reSolve(localr);
	if(Src) Src->multAdd(localr, localc); // localc = - (Krr^-1 Krc)^T fr
	delete [] localr;

	for(int iRbm = 0; iRbm < nrbms; ++iRbm) {
		Scalar *rbm = rbms + size*iRbm;
		for(int i=0; i<numCDofs; ++i) {
			if(cornerEqNums[i] > -1) { e[iRbm] += rbm[cornerEqNums[i]]*localc[i]; } // e += -N^T (Krr^-1 Krc)^T fr
		}
	}
}

// ****************************************************************************************************

template<class Scalar>
Scalar
FetiSub<Scalar>::getMpcRhs(int glMPCnum) const {
	return mpc[globalToLocalMPC[glMPCnum]]->rhs;
}

template<class Scalar>
Scalar
FetiSub<Scalar>::getMpcRhs_primal(int glMPCnum) const {
	return mpc_primal[globalToLocalMPC_primal[glMPCnum]]->rhs;
}

/**
* @param[out] fr Force on remainder DOFs.
* @param[in] uc  Corner displacements.
*/
template<class Scalar>
void
FetiSub<Scalar>::multKrc(Scalar *fr, const Scalar *uc) const {
	const auto &Src = this->Src;
	int numCDofs = Src->numCol();
	Scalar *ucLocal = (Scalar *) dbg_alloca(sizeof(Scalar) * numCDofs);

	for (int i = 0; i < numCDofs; ++i) {
		if (cornerEqNums[i] > -1) ucLocal[i] = uc[cornerEqNums[i]];
		else ucLocal[i] = 0.0;
	}

	// Perform fr = fr - Krc uc
	if (Src) Src->transposeMultSubtract(ucLocal, fr);
}


template<class Scalar>
void
FetiSub<Scalar>::getFr(const Scalar *f, Scalar *fr) const {
	Scalar v[Ave.cols()];
	VectorView<Scalar> t(v, Ave.cols(), 1);
	VectorView<Scalar> l(fr, Ave.rows(), 1);


	for(int dof = 0; dof < ccToC.size(); ++dof)
		fr[dof] = f[ccToC[dof]];

	if (Ave.cols() > 0) { // Averages to zero 072513 JAT
		t = Ave.transpose() * l;
		l.noalias() -= Ave * t;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::getFc(const Scalar *f, Scalar *Fc) const {
	int i, j;
	int dNum[DofSet::max_known_dof];
	int iOff = 0;
	const auto &c_dsa = get_c_dsa();
	for (i = 0; i < numCRN; ++i) {
		int nd = c_dsa->number(cornerNodes[i], cornerDofs[i], dNum);
		for (j = 0; j < nd; ++j) Fc[iOff + j] = f[dNum[j]];
		iOff += nd;
	}

	if (Ave.cols() > 0) { // Average corners 072513 JAT
		int numEquations = this->Krr->neqs();
		Scalar fr[numEquations];
		for(int dof = 0; dof < ccToC.size(); ++dof)
			fr[dof] = f[ccToC[dof]];

		int nAve, nCor, k;
		Scalar s;
		nCor = this->Krc ? this->Krc->numCol() : 0;
		nAve = Src->numCol() - nCor;
		for (int i = 0; i < nAve; ++i) {
			s = 0.0;
			for (int k = 0; k < numEquations; ++k)
				s += Ave[i][k] * fr[k];
			Fc[nCor + i] = s;
		}
	}
}


template<class Scalar>
void
FetiSub<Scalar>::projectActiveIneq(Scalar *v) const {
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpc[locMpcNb]->type == 1 && mpc[locMpcNb]->active)
			v[scomm->mapT(SComm::mpc, i)] = 0.0;
	}
}


// ref: Dostal, Horak and Stefanica (IMACS 2005) "A scalable FETI-DP algorithm for coercive variational inequalities"
// every row of C matrix should have a unit norm to improve the condition number of the Feti operator
template<class Scalar>
void
FetiSub<Scalar>::normalizeCstep1(Scalar *cnorm) {
	auto &mpc = this->mpc;
	for (int i = 0; i < numMPC; ++i)
		for (int j = 0; j < mpc[i]->nterms; ++j)
			cnorm[localToGlobalMPC[i]] += mpc[i]->terms[j].coef * mpc[i]->terms[j].coef;
}

template<class Scalar>
void
FetiSub<Scalar>::normalizeCstep2(Scalar *cnorm) {
	auto &mpc = this->mpc;
	for (int i = 0; i < numMPC; ++i)
		for (int j = 0; j < mpc[i]->nterms; ++j)
			mpc[i]->terms[j].coef /= cnorm[localToGlobalMPC[i]];
}

template<typename Scalar>
void FetiSub<Scalar>::getFw(const Scalar *f, Scalar *fw) const {
	// By default we have no wet interface treatment.
}

template<class Scalar>
void
FetiSub<Scalar>::recvMpcStatus(FSCommPattern<int> *mpcPat, int flag, bool &statusChange) {
	auto &mpc = this->mpc;
	// this function is to make sure that the status of an mpc is the same in all subdomains which share it
	// needed to due to possible roundoff error
	// if flag == 1 then make dual constraint not active in all subdomains if not active in at least one (use in proportioning step)
	// if flag == 0 then make dual constraint active in all subdomains if it is active in at least one (use in expansion step)
	// if flag == -1 use mpc master status in all subdomains
	// note: could use SComm::ieq list
	int i, j;
	bool *tmpStatus = (bool *) alloca(sizeof(bool) * numMPC);
	for (int i = 0; i < numMPC; ++i) tmpStatus[i] = !mpc[i]->active;
	for (i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		if (subNum() != neighb) {
			FSSubRecInfo<int> rInfo = mpcPat->recData(neighb, subNum());
			for (j = 0; j < scomm->lenT(SComm::mpc, i); ++j) {
				int locMpcNb = scomm->mpcNb(i, j);
				if (flag == -1) tmpStatus[locMpcNb] = (rInfo.data[j] > -1) ? bool(rInfo.data[j]) : tmpStatus[locMpcNb];
				else // XXXX
					tmpStatus[locMpcNb] = (flag == 1) ? (tmpStatus[locMpcNb] || bool(rInfo.data[j])) : (
							tmpStatus[locMpcNb] && bool(rInfo.data[j]));
			}
		}
	}

	bool print_debug = false;
	statusChange = false;
	for (i = 0; i < numMPC; ++i) {
		if (getFetiInfo().contactPrintFlag && mpcMaster[i]) {
			if (!mpc[i]->active && !tmpStatus[i]) {
				std::cerr << "-";
				if (print_debug)
					std::cerr << " recvMpcStatus: sub = " << subNum() << ", mpc = " << localToGlobalMPC[i]
					          << std::endl;
			}
			else if (mpc[i]->active && tmpStatus[i]) {
				std::cerr << "+";
				if (print_debug)
					std::cerr << " recvMpcStatus: sub = " << subNum() << ", mpc = " << localToGlobalMPC[i]
					          << std::endl;
			}
		}
		mpc[i]->active = !tmpStatus[i];
		if (mpcStatus2[i] == mpc[i]->active) statusChange = true;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::updateActiveSet(Scalar *v, double tol, int flag, bool &statusChange) {
	// flag = 0 : dual planing
	// flag = 1 : primal planing
	int *chgstatus = (int *) alloca(numMPC * sizeof(int));
	for (int i = 0; i < numMPC; ++i) chgstatus[i] = -1;  // set to 0 to remove, 1 to add

	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpc[locMpcNb]->type == 1) { // inequality constraint requiring planing

			if (flag == 0) { // dual planing
				if (mpcStatus1[locMpcNb] ==
				    1) { // active set expansion only: if constraint was initially active then it will not change status
					// if active and lambda < 0 then remove from active set
					if (mpc[locMpcNb]->active && ScalarTypes::lessThan(v[scomm->mapT(SComm::mpc, i)], tol))
						chgstatus[locMpcNb] = 1;
					// if not active and lambda >= 0 then add to the active set
					if (!mpc[locMpcNb]->active &&
					    ScalarTypes::greaterThanEq(v[scomm->mapT(SComm::mpc, i)], tol))
						chgstatus[locMpcNb] = 0;
				}
			} else { // primal planing
				if (mpcStatus1[locMpcNb] ==
				    0) { // active set contraction only: if constraint was initially inactive then it will not change status
					// if not active and w <= 0 then add to active set
					if (!mpc[locMpcNb]->active && ScalarTypes::lessThanEq(v[scomm->mapT(SComm::mpc, i)], tol))
						chgstatus[locMpcNb] = 0;
					// if active and w > 0 then remove from the active set
					if (mpc[locMpcNb]->active && ScalarTypes::greaterThan(v[scomm->mapT(SComm::mpc, i)], tol))
						chgstatus[locMpcNb] = 1;
				}
			}

		}
	}

	statusChange = false;
	for (int i = 0; i < numMPC; ++i) {
		if (chgstatus[i] > -1) {
			statusChange = true;
			mpc[i]->active = !bool(chgstatus[i]);
			if (getFetiInfo().contactPrintFlag && mpcMaster[i]) {
				if (chgstatus[i] == 0)
					std::cerr << "-";
				else std::cerr << "+";
			}
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::mergeUr(Scalar *ur, Scalar *uc, Scalar *u, Scalar *lambda) {
	int i, iNode;
	int rDofs[DofSet::max_known_dof];
	int oDofs[DofSet::max_known_dof];
	for(int dof = 0; dof < ccToC.size(); ++dof)
		u[ccToC[dof]] = ur[dof];

	int j;
	int iOff = 0;
	for (i = 0; i < numCRN; ++i) {
		int nd = get_c_dsa()->number(cornerNodes[i], cornerDofs[i], oDofs);
		for (j = 0; j < nd; ++j) {
			if (cornerEqNums[iOff + j] > -1)
				u[oDofs[j]] = uc[cornerEqNums[iOff + j]];
			else
				u[oDofs[j]] = 0.0;
		}
		iOff += nd;
	}

	// Primal augmentation 030314 JAT
	if(Ave.cols() > 0) {
		int nCor = this->Krc?this->Krc->numCol() : 0;
		int nAve = Ave.cols();
		for(int i = 0; i < ccToC.size(); ++i)
			for(int j = 0; j < nAve; ++j)
				u[ccToC[i]] = Ave[j][i]*uc[cornerEqNums[nCor+j]];
	}

	// extract uw
	Scalar *uw = (Scalar *) dbg_alloca(numWIdof * sizeof(Scalar));
	for (i = 0; i < scomm->lenT(SComm::wet); ++i)
		uw[scomm->wetDofNb(i)] = lambda[scomm->mapT(SComm::wet, i)];

//	for (i = 0; i < numWInodes; ++i) {
//		DofSet thisDofSet = wetInterfaceDofs[i]; // (*c_dsa)[wetInterfaceNodes[i]];
//		int nd = thisDofSet.count();
//		dsa->number(wetInterfaceNodes[i], thisDofSet, rDofs);
//		c_dsa->number(wetInterfaceNodes[i], thisDofSet, oDofs);
//		for (j = 0; j < nd; ++j) {
//			u[oDofs[j]] = uw[wetInterfaceMap[rDofs[j]]];
//		}
//	}
	if(numWInodes != 0)
		throw "Wet interface is not supported anymore. Work on the above lines if you want it.";

	// keep a local copy of the lagrange multipliers
	setLocalLambda(lambda);
}


template<class Scalar>
void
FetiSub<Scalar>::setLocalLambda(Scalar *_localLambda) {
	if (localLambda) delete[] localLambda;
	localLambda = new double[totalInterfSize];
	for (int i = 0; i < totalInterfSize; ++i) localLambda[i] = ScalarTypes::Real(_localLambda[i]);
}

template<class Scalar>
void
FetiSub<Scalar>::multfc(const VectorView<Scalar> &fr, /*Scalar *fc,*/ const VectorView<Scalar> &lambda) const {
	Scalar v[Ave.cols()];
	VectorView<Scalar> t(v, Ave.cols(), 1);
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> force(localLen());

	force = -fr;
	force -= B*lambda;
	force -= Bm*lambda;


	if (numWIdof) Krw->multAddNew(localw.data(), force.data());  // coupled_dph: force += Krw * uw

	// Primal augmentation 072513 JAT
	if (Ave.cols() > 0) {
		t = this->Eve.transpose() * force;
		force.noalias() -= Ave * t;
	}

	if (this->Krr) this->Krr->solveInPlace(force);

	// Extra orthogonaliztion for stability  072216 JAT
	if (Ave.cols() > 0) {
		t = Ave.transpose() * force;
		force.noalias() -= Ave * t;
	}

	// re-initialization required for mpc/contact
	fcstar.assign(Src->numCol(), 0.0);

	// fcstar = - (Krr^-1 Krc)^T fr
	//        = - Krc^T (Krr^-1 fr)
	//        = Src force
	if (Src) Src->multAdd(force.data(), fcstar.data());

	// for coupled_dph add fcstar -= Kcw Bw uw
	if (numWIdof) {
		if (Kcw) Kcw->mult(localw.data(), fcstar.data(), -1.0, 1.0);
		if (Kcw_mpc) Kcw_mpc->multSubWI(localw.data(), fcstar.data());
	}

	VectorView<Scalar> fcs(fcstar.data(), fcstar.size());
	// add Bc^(s)^T lambda
	fcs += Bc*lambda;
}

template<class Scalar>
void
FetiSub<Scalar>::multAddBrT(const Scalar *interfvec, Scalar *localvec, Scalar *uw) const {
	VectorView<Scalar> locF(localvec, localLen());
	VectorView<const Scalar> lambda(interfvec, totalInterfSize);
	VectorView<Scalar> w(uw, numWIdof);
	locF += B*lambda;
	w -= Bw*lambda;
	locF += Bm*lambda;

	// coupled_dph: localvec -= Krw * uw
	if (Krw) Krw->multAddNew(uw, localvec);
}

template<class Scalar>
void
FetiSub<Scalar>::multBr(const Scalar *localvec, Scalar *interfvec, const Scalar *_uc, const Scalar *uw) const {
	VectorView<const Scalar> loc_u(localvec, localLen());
	VectorView<Scalar> u_interf(interfvec, totalInterfSize);
	VectorView<const Scalar> w(uw, numWIdof);
	VectorView<const Scalar> uc(_uc, Src->numCol());

	u_interf = B.transpose()*loc_u;
	u_interf -= Bw.transpose()*w;
	u_interf += Bc.transpose()*uc;
}

template<class Scalar>
void
FetiSub<Scalar>::fetiBaseOp(Scalar *uc, GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec) const {
	Scalar v[this->Ave.cols()];
	VectorView<Scalar> t(v, this->Ave.cols(), 1);
	VectorView<Scalar> l(localvec, this->Ave.rows(), 1);

	// localvec += Br^T * interfvec
	multAddBrT(interfvec, localvec);

	// Primal augmentation 072513 JAT
	if (this->Ave.cols() > 0) {
		t = this->Eve.transpose() * l;
		l.noalias() -= this->Ave * t;
	}

	// localvec = Krr^-1 * localvec
	if (s) s->reSolve(localvec);

	// Extra orthogonaliztion for stability  072216 JAT
	if (this->Ave.cols() > 0) {
		t = this->Ave.transpose() * l;
		l.noalias() -= this->Ave * t;
	}

	// interfvec = Br * localvec
	multBr(localvec, interfvec, uc);
}

template<class Scalar>
void
FetiSub<Scalar>::fetiBaseOpCoupled1(GenSolver<Scalar> *s, Scalar *localvec, const Scalar *interfvec,
                                         FSCommPattern<Scalar> *wiPat) const {
	auto &localw = this->localw;
	// localvec += Br^T * interfvec
	multAddBrT(interfvec, localvec, localw.data());

	// solve for localvec
	if (s) s->reSolve(localvec);

	if (numWIdof) {
		int i, j;
		// compute Kww uw for this subdomain
		for (i = 0; i < numWIdof; ++i) localw_copy[i] = 0.0;
		this->Kww->mult(localw.data(), localw_copy.data());  // localw_copy = - Kww * uw

		// compute Kww uw to send to neighbors
//		if (getFetiInfo().fsi_corner == 0)
		// TODO Bring this back with wweight
		if(false)
			for (i = 0; i < scomm->numT(SComm::fsi); ++i) {
				if (subNum() != scomm->neighbT(SComm::fsi, i)) {
					FSSubRecInfo<Scalar> sInfo = wiPat->getSendBuffer(subNum(), scomm->neighbT(SComm::fsi, i));
					for (j = 0; j < numNeighbWIdof[i]; ++j) sInfo.data[j] = 0.0;
					neighbKww->multAdd(localw.data(), sInfo.data, glToLocalWImap, neighbGlToLocalWImap[i]);
				} else {
					neighbKww->multAdd(localw.data(), localw_copy.data(), glToLocalWImap);
				}
			}
	}
}


template<class Scalar>
void
FetiSub<Scalar>::fetiBaseOpCoupled2(const Scalar *uc, const Scalar *localvec, Scalar *interfvec,
                                         FSCommPattern<Scalar> *wiPat, const Scalar *fw) const {
	// coupled_dph
	if (numWIdof) {
		int i, j;
		auto &Krw = this->Krw;
		auto &Kcw = this->Kcw;
		auto &Kcw_mpc = this->Kcw_mpc;


//		// TODO Bring this back with wweight
//		if (getFetiInfo().fsi_corner == 0)
//			for (i = 0; i < scomm->numT(SComm::fsi); ++i) {
//				if (subNum() != scomm->neighbT(SComm::fsi, i)) {
//					FSSubRecInfo<Scalar> rInfo = wiPat->recData(scomm->neighbT(SComm::fsi, i), subNum());
//					for (j = 0; j < numWIdof; ++j) localw_copy[j] += rInfo.data[j] / wweight[j];
//				}
//			}

		if (Krw) Krw->transposeMultSubNew(localvec, localw_copy.data()); // localw_copy -= Krw^T * localvec
		int numCDofs = Src->numCol();
		Scalar *ucLocal = (Scalar *) dbg_alloca(sizeof(Scalar) * numCDofs);
		for (i = 0; i < numCDofs; ++i) {
			if (cornerEqNums[i] > -1) ucLocal[i] = uc[cornerEqNums[i]];
			else ucLocal[i] = 0.0;
		}
		if (Kcw) Kcw->trMult(ucLocal, localw_copy.data(), -1.0, 1.0);  // localw_copy -= Kcw^T Bc uc
		if (Kcw_mpc) Kcw_mpc->transposeMultSubtractWI(ucLocal, localw_copy.data());
		if (fw) for (i = 0; i < numWIdof; ++i) localw_copy[i] += fw[i];  // localw_copy += fw
	}

	// interfvec = Br * localvec
	multBr(localvec, interfvec, uc, localw_copy.data());
}

template<class Scalar>
void
FetiSub<Scalar>::multFcB(Scalar *p) {
	int i, k;
	if (Src->numCol() == 0) return;
	if ((totalInterfSize == 0) || (localLen() == 0)) {
		for (i = 0; i < Src->numCol(); ++i) fcstar[i] = 0.0;
		return;
	}

	// fcstar = - (Krr^-1 Krc)^T p
	//        = - Krc^T Krr^-1 p
	//        = - Acr p
	// TODO Change to fully use Eigen.
	GenStackFullM<Scalar> Acr(Src->numCol(), totalInterfSize, BKrrKrc.data());
	fcstar.resize(Src->numCol());
	Acr.mult(p, fcstar.data(), -1.0, 0.0);

	// for coupled_dph add fcstar += Kcw Bw uw
	if (numWIdof) {
		for (i = 0; i < scomm->lenT(SComm::wet); ++i)
			localw[scomm->wetDofNb(i)] = p[scomm->mapT(SComm::wet, i)];
		if (Kcw) Kcw->mult(localw.data(), fcstar.data(), -1.0, 1.0);
		if (Kcw_mpc) Kcw_mpc->multSubWI(localw.data(), fcstar.data());
	}

	// fcstar += Bc^(s)^T p
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
	for (i = 0; i < numMPC; ++i) mpcFlag[i] = true;
	for (i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpcFlag[locMpcNb]) {
			const auto &m = mpc[locMpcNb];
			for (k = 0; k < m->nterms; k++) {
				int dof = (m->terms)[k].dof;
				if ((dof >= 0) && (cornerMap[dof] >= 0))
					fcstar[cornerMap[dof]] += p[scomm->mapT(SComm::mpc, i)] * (m->terms)[k].coef;
			}
			mpcFlag[locMpcNb] = false;
		}
	}
}


template<class Scalar>
void
FetiSub<Scalar>::multKcc() {
	auto &BKrrKrc = this->BKrrKrc;
	auto &Krw = this->Krw;
	auto &localw = this->localw;
	// resize Kcc if necessary to add space for "primal" mpcs and augmentation
	if (this->Src->numCol() != this->Kcc->dim()) {
		std::unique_ptr<GenAssembledFullM<Scalar>> Kcc_copy = std::move(this->Kcc);
		this->Kcc = std::make_unique<GenAssembledFullM<Scalar>>(this->Src->numCol(), cornerMap);
		for (int i = 0; i < numCRNdof; ++i) for (int j = 0; j < numCRNdof; ++j) (*this->Kcc)[i][j] = (*Kcc_copy)[i][j];
	}


	// add in MPC coefficient contributions for Kcc^(s).
	if (numMPC_primal > 0)
		assembleMpcIntoKcc();

	if (this->Krr == 0) return;

	// Kcc* -> Kcc - Krc^T Krr^-1 Krc

	// first, perform Krc I = iDisp
	// which extracts the correct rhs vectors for the forward/backwards

	int nRHS = this->Src->numCol();
	Scalar **iDisp = new Scalar *[nRHS];
	Scalar *firstpointer = new Scalar[nRHS * nRHS];

	int numEquations = this->Krr->neqs();
	BKrrKrc.resize(totalInterfSize, nRHS);
	Scalar *thirdpointer = new Scalar[nRHS * numEquations];
	Scalar **KrrKrc = (Scalar **) dbg_alloca(nRHS * sizeof(Scalar *));
	//if(nRHS*numEquations == 0)
	//  fprintf(stderr, "We have a zero size %d %d %d\n",numEquations,totalInterfSize,nRHS);

	int iRHS, iDof;
	for (iRHS = 0; iRHS < nRHS; ++iRHS) {
		iDisp[iRHS] = firstpointer + iRHS * nRHS;
		KrrKrc[iRHS] = thirdpointer + iRHS * numEquations;
		for (iDof = 0; iDof < numEquations; ++iDof)
			KrrKrc[iRHS][iDof] = 0.0;
	}
	if (this->Src) this->Src->multIdentity(KrrKrc);

	// 070213 JAT
	if (this->Src && (getFetiInfo().augmentimpl == FetiInfo::Primal)) {
		int nAve, nCor;
		nCor = this->Krc ? this->Krc->numCol() : 0;
		nAve = this->Src->numCol() - nCor;
		if (nAve) {
			int i, j, k, nz;
			Scalar *pKve = new Scalar[nAve * numEquations];
			Scalar *pv = new Scalar[nAve];
			GenFullM<Scalar> AKA(nAve);
			Scalar s, *pAKA;
			Scalar **Kve = new Scalar *[nAve];
			this->Ave.resize(numEquations, nAve);
			for (i = 0; i < nAve; ++i) {
				Kve[i] = pKve + i * numEquations;
			}
			for (i = 0; i < nAve; ++i) {
				s = 0.0;
				for (j = 0; j < numEquations; ++j) {
					this->Ave[i][j] = KrrKrc[nCor + i][j];
					s += this->Ave[i][j] * this->Ave[i][j];
				}
				s = 1.0 / sqrt(s);
				for (j = 0; j < numEquations; ++j)
					this->Ave[i][j] *= s;
			}
			for (i = 0; i < nAve; ++i) {
				for (j = 0; j < numEquations; ++j)
					this->Eve[i][j] = this->Ave[i][j];
				this->Krr->reSolve(nAve, this->Eve[i]);
			}
			pAKA = AKA.data();

			Tgemm('T', 'N', nAve, nAve, numEquations, 1.0, this->Eve[0], numEquations,
			      this->Ave[0], numEquations, 0.0, pAKA, nAve);

			AKA.factor();
			for (j = 0; j < numEquations; ++j) {
				for (i = 0; i < nAve; ++i)
					pv[i] = this->Eve[i][j];
				AKA.reSolve(pv);
				for (i = 0; i < nAve; ++i)
					this->Eve[i][j] = pv[i];
			}

			for (i = 0; i < nAve; ++i)
				for (j = 0; j < nCor; ++j) {
					s = 0.0;
					for (k = 0; k < numEquations; ++k)
						s += KrrKrc[j][k] * this->Ave[i][k];
					(*this->Kcc)[nCor + i][j] = s;
					(*this->Kcc)[j][nCor + i] = s;
				}
			for (i = 0; i < nAve; ++i) {
				this->KrrSparse->mult(this->Ave[i], Kve[i]);
				for (j = 0; j < nAve; ++j) {
					s = 0.0;
					for (k = 0; k < numEquations; ++k)
						s += Kve[i][k] * this->Ave[j][k];
					(*this->Kcc)[nCor + i][nCor + j] = s;
				}
			}

			nz = 0;
			for (i = 0; i < nAve; ++i)
				for (k = 0; k < numEquations; ++k)
					if (std::abs(Kve[i][k]) > 0.0) nz++;

			int *KACount = new int[nAve];
			int *KAList = new int[nz];
			Scalar *KACoefs = new Scalar[nz];
			nz = 0;
			for (i = 0; i < nAve; ++i) {
				KACount[i] = 0;
				for (k = 0; k < numEquations; ++k)
					if (std::abs(Kve[i][k]) > 0.0) {
						KACount[i]++;
						KAList[nz] = k;
						KACoefs[nz] = Kve[i][k];
						nz++;
					}
			}

			this->Grc = std::make_unique<GenCuCSparse<Scalar>>(nAve, numEquations, KACount, KAList, KACoefs);

			if (this->Src->num() == 2)
				this->Src->setSparseMatrix(1, this->Grc.get());
			else if (this->Src->num() == 1 && nCor == 0)
				this->Src->setSparseMatrix(0, this->Grc.get());
			else {
				fprintf(stderr, "unsupported number of blocks in Src\n");
				exit(1);
			}

			delete[] KACount;

			for (i = 0; i < nRHS; ++i)
				for (j = 0; j < numEquations; ++j)
					KrrKrc[i][j] = 0.0;

			this->Src->multIdentity(KrrKrc);

			Scalar vt[nAve];
			VectorView<Scalar> v{vt, nAve};
			for (j = 0; j < nRHS; j++) {
				v = this->Eve * VectorView<Scalar>{KrrKrc[j], numEquations};
				VectorView<Scalar>{KrrKrc[j], numEquations} -= this->Ave * v;
			}
		}
	}

	// KrrKrc <- Krr^-1 Krc
	if (this->Krr) this->Krr->reSolve(nRHS, KrrKrc); // this can be expensive when nRHS is large eg for coupled

	// -Krc^T KrrKrc
	for (iRHS = 0; iRHS < nRHS; ++iRHS)
		for (iDof = 0; iDof < nRHS; ++iDof)
			iDisp[iRHS][iDof] = 0.0;

	// Multiple RHS version of multSub: iDisp <- -Krc^T KrrKrc
	if (this->Src) this->Src->multSub(nRHS, KrrKrc, iDisp);

	if (this->Kcc) this->Kcc->add(iDisp);

	delete[] iDisp;
	delete[] firstpointer;
	auto &mpc = this->mpc;

	int k;
	for (iRHS = 0; iRHS < nRHS; ++iRHS) {
		bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
		for (int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
		bool *wiFlag = (bool *) dbg_alloca(sizeof(bool) * numWIdof);
		for (int i = 0; i < numWIdof; ++i) wiFlag[i] = true;

		if (Krw) Krw->transposeMultNew(KrrKrc[iRHS], localw.data());

		for (iDof = 0; iDof < totalInterfSize; iDof++) {
			switch (boundDofFlag[iDof]) {
				case 0:
					BKrrKrc(iDof, iRHS) = KrrKrc[iRHS][allBoundDofs[iDof]];
					break;
				case 1: { // wet interface
					int windex = -1 - allBoundDofs[iDof];
					if (wiFlag[windex]) {
						BKrrKrc(iDof, iRHS) = -localw[-1 - allBoundDofs[iDof]];
						wiFlag[windex] = false;
					} else BKrrKrc(iDof, iRHS) = 0.0;
				}
					break;
				case 2: { // dual mpc
					int locMpcNb = -1 - allBoundDofs[iDof];
					const auto &m = mpc[locMpcNb];
					BKrrKrc(iDof, iRHS) = 0.0;
					if (mpcFlag[locMpcNb]) {
						for (k = 0; k < m->nterms; k++) {
							int cc_dof = (m->terms)[k].ccdof;
							if (cc_dof >= 0) BKrrKrc(iDof, iRHS) += KrrKrc[iRHS][cc_dof] * (m->terms)[k].coef;
						}
						mpcFlag[locMpcNb] = false;
					}
				}
					break;
			}
		}
	}
	delete[] thirdpointer;
}


template<class Scalar>
void
FetiSub<Scalar>::assembleMpcIntoKcc() {
	// compute mpcOffset for my Subdomain
	int mpcOffset = numCRNdof;

	int i, iMPC;
	for (iMPC = 0; iMPC < numMPC_primal; ++iMPC) {
		for (i = 0; i < mpc_primal[iMPC]->nterms; ++i) {
			int d = mpc_primal[iMPC]->terms[i].dof;
			int dof = mpc_primal[iMPC]->terms[i].ccdof;
			if ((dof < 0) && (d >= 0) && !isWetInterfaceDof(d)) {
				int column = cornerMap[d];
				int row = mpcOffset + iMPC;
				if (row > this->Kcc->dim())
					std::cout << " *** ERROR: Dimension Error Row = " << row << " > " << this->Kcc->dim() << std::endl;
				if (column > this->Kcc->dim())
					std::cout << " *** ERROR: Dimension Error Col = " << column << " > " << this->Kcc->dim()
					          << std::endl;
				// i.e. an mpc touches a corner node that also has DBCs
				if (column >= 0) {
					(*this->Kcc)[row][column] += mpc_primal[iMPC]->terms[i].coef;
					(*this->Kcc)[column][row] += mpc_primal[iMPC]->terms[i].coef;
				}
			}
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::constructKcc() {
	Kcc = std::make_unique<GenAssembledFullM<Scalar>>(numCRNdof, cornerMap);
//	memK += numCRNdof * numCRNdof; // TODO Move/duplicate memory use variable ???
}

template<class Scalar>
void
FetiSub<Scalar>::initMpcStatus() {
	for (int i = 0; i < numMPC; ++i) {
		mpc[i]->active = false;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::saveMpcStatus() {
	// this saves the status before first update iteration so it can be reset if nonmonotic
	if (!mpcStatus) mpcStatus = new int[numMPC];
	for (int i = 0; i < numMPC; ++i) {
		mpcStatus[i] = int(!mpc[i]->active);
	}
}

template<class Scalar>
void
FetiSub<Scalar>::restoreMpcStatus() {
	for (int i = 0; i < numMPC; ++i) {
		if (getFetiInfo().contactPrintFlag && mpcMaster[i]) {
			if (!mpc[i]->active && !mpcStatus[i]) std::cerr << "-";
			else if (mpc[i]->active && mpcStatus[i]) std::cerr << "+";
		}
		mpc[i]->active = bool(!mpcStatus[i]);
	}
}

template<class Scalar>
void
FetiSub<Scalar>::saveMpcStatus1() const {
	auto &mpc = this->mpc;
	if (!mpcStatus1) mpcStatus1 = new bool[numMPC];
	for (int i = 0; i < numMPC; ++i) mpcStatus1[i] = !mpc[i]->active;
}

template<class Scalar>
void
FetiSub<Scalar>::saveMpcStatus2() {
	auto &mpc = this->mpc;
	if (!mpcStatus2) mpcStatus2 = new bool[numMPC];
	for (int i = 0; i < numMPC; ++i) mpcStatus2[i] = !mpc[i]->active;
}

template<class Scalar>
void
FetiSub<Scalar>::cleanMpcData() {
	if (mpcStatus) {
		delete[] mpcStatus;
		mpcStatus = 0;
	}
	if (mpcStatus1) {
		delete[] mpcStatus1;
		mpcStatus1 = 0;
	}
	if (mpcStatus2) {
		delete[] mpcStatus2;
		mpcStatus2 = 0;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::makeKccDofsExp2(int nsub, FetiBaseSub **sd,
                                 int augOffset, Connectivity *subToEdge) {
	int numC = numCoarseDofs();
	cornerEqNums.resize(numC);

	// numbers the corner equations
	int offset = 0;
	for (int i = 0; i < numCRN; ++i) {
		int offset2 = 0;
		for (int j = 0; j < nsub; ++j) {
			auto &nodeMap = sd[j]->getGlobalToLocalNode();
			ConstrainedDSA *cornerEqs = sd[j]->get_c_dsa();
			if (nodeMap[glCornerNodes[i]] > -1) {
				int count = cornerEqs->number(nodeMap[glCornerNodes[i]], cornerDofs[i].list(),
				                              cornerEqNums.data() + offset);
				for (int k = 0; k < count; ++k) cornerEqNums[offset + k] += offset2;
				offset += count;
				break;
			}
			offset2 += cornerEqs->size();
		}
	}

	if (getFetiInfo().augmentimpl == FetiInfo::Primal) {
		int iEdgeN = 0;
		for (int iNeighb = 0; iNeighb < scomm->numNeighb; ++iNeighb) {
			if (scomm->isEdgeNeighb[iNeighb]) {
				int k = augOffset + (*subToEdge)[subNum()][iEdgeN];
				int offset2 = 0;
				for (int iSub = 0; iSub < nsub; ++iSub) {
					GlobalToLocalMap &nodeMap = sd[iSub]->getGlobalToLocalNode();
					ConstrainedDSA *cornerEqs = sd[iSub]->get_c_dsa();
					if (nodeMap[k] > -1) {
						int fDof = cornerEqs->firstdof(nodeMap[k]);
						int count = edgeDofSize[iNeighb];
						for (int k = 0; k < count; ++k)
							cornerEqNums[offset + k] = fDof + k + offset2;
						offset += count;
						break;
					}
					offset2 += cornerEqs->size();
				}
				iEdgeN++;
			}
		}
	}
}

template class FetiSub<double>;
template class FetiSub<std::complex<double>>;