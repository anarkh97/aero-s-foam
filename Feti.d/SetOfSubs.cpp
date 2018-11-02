//
// Created by Michel Lesoinne on 2018-11-01.
//

#include <complex>
#include "Threads.d/PHelper.h"
#include "SetOfSubs.h"

template<typename Scalar>
void SetOfSubs<Scalar>::getSharedDOFs()
{
	auto myCPU = communicator->cpuNum();
	auto numSub = communicator->size();
	FSCommPattern<int> nodeIntPat(communicator, cpuToSub.get(), myCPU, FSCommPattern<int>::CopyOnSend);
	for (auto &sub: subDomain)
		sub->setNodeCommSize(&nodeIntPat);
	nodeIntPat.finalize();

	paralApplyToAll(numSub, subDomain, &FetiBaseSub::sendDOFList, &nodeIntPat);
	nodeIntPat.exchange();
	paralApply(subDomain, &FetiSub<Scalar>::gatherDOFList, &nodeIntPat);
	paralApply(subDomain, &FetiBaseSub::gatherDOFListPlus, &nodeIntPat);
}

template
class SetOfSubs<double>;