//
// Created by Michel Lesoinne on 4/12/18.
//

#ifndef FEM_CORNERSELECTOR_H
#define FEM_CORNERSELECTOR_H


#include <vector>
#include "FetiSub.h"

class FetiSubCornerHandler;
class Connectivity;
class FSCommunicator;
class CoordSet;
class DofSetArray;
class ConstrainedDSA;
class BaseSub;
template <class Type> class FSCommPattern;

class CornerSelector
{
	int glNumSub;
	int nSub;
	Connectivity *grToSub;
	int *glSubGroup;
	std::vector<FetiSubCornerHandler*> cornerHandler;
	FSCommPattern<int> *cpat;
	FSCommunicator *communicator;
	int dim;
	int dims[4];

	void chooseCorners(char *glCornerList, double (*xyz)[3],
	                   Connectivity &cNConnect, Connectivity &subToRotCrn,
	                   int *glCrnGroup);
public:
	CornerSelector(int nGlobSub, int nLocSub, std::vector<FetiSubCornerHandler *> handlers,
	               FSCommPattern<int> *commPattern, FSCommunicator *communicator);
	~CornerSelector();
	int makeCorners();
	Connectivity *getGrToSub() { return grToSub; }
	auto &handlers() { return cornerHandler; }
};

// This class is the algorithmic class to select corner nodes
class FetiSubCornerHandler
{

public:
	FetiSubCornerHandler(int sub, int nn, CoordSet &n, Connectivity &nTn, DofSetArray &d, Connectivity &sh, int *nsb,
		                     ConstrainedDSA *c_dsa, FetiBaseSub *_subPre);
	~FetiSubCornerHandler();
	void markMultiDegNodes();
	void dispatchSafeNodes(FSCommPattern<int> *);
	void markSafeNodes(FSCommPattern<int> *);
	void dispatchRotCorners(FSCommPattern<int> *);
	void markRotCorners(FSCommPattern<int> *);
	void pickAnyCorners();
	void countAndMarkCornerCand(int *mync, int *totnc);
	void getCornerXYZ(int *, double (*)[3], char *essential, int *cTsP, int *cTsT);
	void dispatchNumbering(FSCommPattern<int> *pat, char *crnMrk,
	                       int *allOrigFC, int *allNewFC, int, int *cntOff);
	void dispatchInitialNumbering(FSCommPattern<int> *pat, int *firstC);
	void recNumbering(FSCommPattern<int> *, int *fM);
	void recInitialNumbering(FSCommPattern<int> *pat, int *numRotCrn);
	void listRotCorners(int *fN, int *crnNum);
	void countContact(int *, char *crnMrk);
	void markDims(int *_dims);

	void resendNumbers(FSCommPattern<int> *pat);
	void checkNumbers(FSCommPattern<int> *pat);

	int *getCorners() { return crnList; }
	int getNumCorners() { return totNC; }

protected:
	int glSubNum;
	int nnodes;
	int *deg;
	int *crnList;
	int totNC;
	int *weight;
	std::vector<bool> isCorner;
	std::vector<bool> isSafe;
	std::vector<bool> glSafe;
	int nNeighb;
	Connectivity &sharedNodes;
	Connectivity &nToN;
	int *neighbSubs;
	CoordSet &nodes;
	DofSetArray &dsa;
	std::vector<bool> isRotMidSideNode;
	int dim;
	int dims[4];
	bool allSafe;
	int nTC; // Total number of corner candidates for this sub
	bool checkForColinearCrossPoints(int numCornerPoints,
	                                 int *localCornerPoints);
	using bit_iterator = std::vector<bool>::const_iterator;
	bool addPotCornerPoints(int numShared, int *allNodes,
	                        bit_iterator isSafe);
	bool mixed; // true if subdomain has active fluid and structure dofs
	FetiBaseSub *subPre;
};

#endif //FEM_CORNERSELECTOR_H
