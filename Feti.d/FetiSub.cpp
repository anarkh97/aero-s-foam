//
// Created by Michel Lesoinne on 12/6/17.
//
#include <complex>
#include <Utils.d/dofset.h>
#include "FetiSub.h"
#include <Driver.d/Mpc.h>
#include <Utils.d/SolverInfo.h>

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

template class FetiSub<double>;
template class FetiSub<std::complex<double>>;