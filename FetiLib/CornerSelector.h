//
// Created by Michel Lesoinne on 4/3/18.
//

#ifndef FEM_CORNERSELECTOR_H
#define FEM_CORNERSELECTOR_H

#include <Utils.d/Connectivity.h>
#include <Driver.d/Communicator.h>

namespace FetiLib {



class SubCornerHelper {

};

class CornerSelector {
public:
	CornerSelector(int nGlobSub, std::vector<SubCornerHelper> scHelpers,
	               FSCommPattern<int> *cpat, FSCommunicator *communicator);
	~CornerSelector() = default;
	
	int makeCorners();
	Connectivity *getGrToSub() { return grToSub; }
private:
	int glNumSub;
	Connectivity *grToSub;
	std::vector<int> glSubGroup;
	std::vector<SubCornerHelper> cornerHelpers;
	FSCommPattern<int> *cpat;
	FSCommunicator *communicator;
	int dim;
	int dims[4];

	void chooseCorners(char *glCornerList, double (*xyz)[3],
	                   Connectivity &cNConnect, Connectivity &subToRotCrn,
	                   int *glCrnGroup);
};
}


#endif //FEM_CORNERSELECTOR_H
