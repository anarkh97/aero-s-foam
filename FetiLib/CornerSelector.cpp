//
// Created by Michel Lesoinne on 4/3/18.
//

#include "CornerSelector.h"
namespace FetiLib {

CornerSelector::CornerSelector(int nGlobSub, std::vector<FetiLib::SubCornerHelper> scHelpers,
                               FSCommPattern<int> *cpat, FSCommunicator *fsCommunicator) :
	cornerHelpers(std::move(scHelpers)),
	glNumSub(nGlobSub),
	cpat(cpat),
	communicator(fsCommunicator)
{

}

}
