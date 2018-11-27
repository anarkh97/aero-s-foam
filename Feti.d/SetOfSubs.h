//
// Created by Michel Lesoinne on 2018-11-01.
//

#ifndef FEM_SETOFSUBS_H
#define FEM_SETOFSUBS_H

#include <vector>
#include "FetiSub.h"

template <typename Scalar>
class SetOfSubs {
public:
private:
public:
	SetOfSubs(FSCommunicator *communicator, std::vector<std::unique_ptr<FetiSub<Scalar> >> subDomain,
	          const std::shared_ptr<Connectivity> &cpuToSub);

private:
	void getSharedDOFs();

	FSCommunicator *communicator;
	std::vector<std::unique_ptr<FetiSub<Scalar> >>  subDomain;
	std::shared_ptr<Connectivity> cpuToSub;
};




#endif //FEM_SETOFSUBS_H
