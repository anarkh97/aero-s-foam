#ifndef _CORNER_MAKER_H_
#define _CORNER_MAKER_H_

#include <vector>
#include <Feti.d/CornerSelector.h>

class Connectivity;
class FSCommunicator;
class Elemset;
class CoordSet;
class DofSetArray;
class ConstrainedDSA;
class BaseSub;
template <class Type> class FSCommPattern;

class CornerMaker
{
	int glNumSub;
	int nSub;
	Connectivity *grToSub;
	int *glSubGroup;
	FetiSubCornerHandler **cornerHandler;
	FSCommPattern<int> *cpat;
	FSCommunicator *communicator;
	int dim;
	int dims[4];

	void chooseCorners(char *glCornerList, double (*xyz)[3],
	                   Connectivity &cNConnect, Connectivity &subToRotCrn,
	                   int *glCrnGroup);
public:
	CornerMaker(int nGlobSub, int nLocSub, FetiSubCornerHandler **,
	            FSCommPattern<int> *, FSCommunicator *);
	~CornerMaker();
	int makeCorners();
	Connectivity *getGrToSub() { return grToSub; }
};

#endif
