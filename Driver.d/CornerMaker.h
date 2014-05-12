#ifndef _CORNER_MAKER_H_
#define _CORNER_MAKER_H_

class SubCornerHandler;
class Connectivity;
class FSCommunicator;
class ElemSet;
class CoordSet;
class DofSetArray;
class BaseSub;
template <class Type> class FSCommPattern;

class CornerMaker 
{
  int glNumSub;
  int nSub;
  Connectivity *grToSub;
  int *glSubGroup;
  SubCornerHandler **cornerHandler;
  FSCommPattern<int> *cpat;
  FSCommunicator *communicator;
  int dim; 
  int dims[4]; 

  void chooseCorners(char *glCornerList, double (*xyz)[3],
                     Connectivity &cNConnect, Connectivity &subToRotCrn,
                     int *glCrnGroup);
 public:
  CornerMaker(int nGlobSub, int nLocSub, SubCornerHandler **,
              FSCommPattern<int> *, FSCommunicator *);
  ~CornerMaker();
  int makeCorners();
  Connectivity *getGrToSub() { return grToSub; }
};

// This class is the algorithmic class to select corner nodes
class SubCornerHandler 
{
  int glSubNum;
  int nnodes, numEle;
  int *deg;
  int *crnList;
  int totNC;
  int *weight;
  bool *isCorner;
  bool *isSafe;
  bool *glSafe;
  int nNeighb;
  Connectivity &sharedNodes;
  Connectivity &nToN;
  int *neighbSubs;
  CoordSet &nodes;
  Elemset &eles;
  DofSetArray &dsa;
  bool *isRotMidSideNode; 
  int dim;
  int dims[4]; 
  bool allSafe;
  int nTC; // Total number of corner candidates for this sub
  bool checkForColinearCrossPoints(int numCornerPoints,
                                   int *localCornerPoints);
  bool addPotCornerPoints(int numShared, int *allNodes, bool *isSafe);
  bool mixed; // true if subdomain has active fluid and structure dofs
  BaseSub *subPre;
 public:
  SubCornerHandler(int sub, int nn, CoordSet &n, int nele, Elemset &ele,
                   Connectivity &nTn, DofSetArray &d,
		   Connectivity &sh, int *nsb, ConstrainedDSA *c_dsa, BaseSub *_subPre);
  ~SubCornerHandler();
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
};

#endif
