#ifndef _SCOMM_H_
#define _SCOMM_H_

#include <Utils.d/MyComplex.h>
#include <Utils.d/Connectivity.h>

class Connectivity;
class BaseSub;
template<class Scalar> class GenSubDomain;

// Subdomain Communication class
class SComm
{
 public:
  int locSubNum;     // when running in distributed mode
  int *glSubToLocal; // mapping for distributed mode

  BaseSub **neighb;  // pointers to neighbors
  int *remoteId;     // id of this subdomain in the corresponding neighbor
  int numNeighb;     // Number of neighbors with shared nodes
  int *subNums;      // identification numbers of neighbors with shared nodes
  Connectivity *sharedNodes; // nodes shared with the neighbors (wet and dry but no virtual)
  void **exchangeData;

  int numEdgeNeighb; // number of neighbors that have at least one non-virtual node on the interface
  bool *isEdgeNeighb;

  Connectivity *sharedDOFs;  // DOFs # shared (all of types 0, 1 and 2)
  Connectivity *sharedDOFsPlus; // also includes corner dofs

  template <class Scalar>
    SComm(int numNeighb, int *subNums, GenSubDomain<Scalar> **neighb = 0, int *remoteId = 0, Connectivity *sharedNodes = 0);
  ~SComm();

  void setExchangeData(int iSub, void *data);
  void *getExchangeData(int iSub);
  void setEdgeNeighb(int _numEdgeNeighb, bool *_isEdgeNeighb);

  template <class Scalar>
    GenSubDomain<Scalar> *getNeighb(int i) {
      return static_cast<GenSubDomain<Scalar> *>(neighb[i]);
    }

  bool isLocal(int iSub);  // PJSA 6-3-04

  // new lists for different types of interface dofs
  // type 0: regular lagrange multipliers
  //      1: wet interface (primal interface unknowns)
  //      2: dual mpcs
  //      3: merged list of types 0, 1 and 2
  //      4: inequality constraints, subset of 2
  //      5: free constraints, subset of 2
  //      6: fluid-structure interactions
  //      7: inter-group multipliers of types 0, 1 and 2
  enum DofType { std, wet, mpc, all, ieq, free, fsi, interg };

 private:
  int numDofType;  // current default is 6, set in SComm constructor
  int *NumNeighb;
  int **SubNums;
  int **TypeMap;  // maps from type-specific interface to general interface
  Connectivity **SharedDOFs;

 public:
  void setNumDofType(int _numDofType) { numDofType = _numDofType; }
  // function to make type-specific lists from *sharedDOFs combined list using boundDofFlag
  //void makeTypeSpecificLists(int *boundDofFlag);
  void print(DofType t) { 
    cerr << "NumNeighb = " << NumNeighb[t] << endl;
    cerr << "SubNums = "; for(int i=0; i<NumNeighb[t]; ++i) cerr << SubNums[t][i] << " "; cerr << endl;
    cerr << "SharedDOFs = \n"; SharedDOFs[t]->print(); 
  }
  Connectivity *getTypeSpecificList(DofType type) { return SharedDOFs[type]; }

  // function to set any one individual type-specific list
  // if this function is used to set the individual types then mergeTypeSpecificLists() must be called
  void setTypeSpecificList(DofType type, int *_subNums, Connectivity *_sharedDOFs);
  int* mergeTypeSpecificLists();
  void setTypeMap(DofType t, int *map);

  // functions to access any of the type-specific lists
  int numT(DofType type) { return SharedDOFs[type]->csize(); }
  int neighbT(DofType type, int iNeighb) { return SubNums[type][iNeighb]; }
  int mapT(DofType type, int iDof) { return TypeMap[type][iDof]; }
  int mapT(DofType type, int iNeighb, int jDof) { return TypeMap[type][SharedDOFs[type]->offset(iNeighb)+jDof]; }
  int lenT(DofType type) { return SharedDOFs[type]->numConnect(); }
  int lenT(DofType type, int iNeighb) { return SharedDOFs[type]->num(iNeighb); }
  int boundDofT(DofType type, int iDof) { return (*SharedDOFs[type])[0][iDof]; }
  int boundDofT(DofType type, int iNeighb, int jDof) { return (*SharedDOFs[type])[iNeighb][jDof]; }
  int offsetT(DofType type, int iNeighb) { return SharedDOFs[type]->offset(iNeighb); }
  int offsetT(DofType type, int iNeighb, int jDof) { return SharedDOFs[type]->offset(iNeighb)+jDof; }
  int* boundDofsT(DofType type) { return (*SharedDOFs[type])[0]; }
  int* boundDofsT(DofType type, int iNeighb) { return (*SharedDOFs[type])[iNeighb]; }
  int* neighbsT(DofType type) { return SubNums[type]; }

  // standard (boundary) dofs helper functions
  int stdDofNb(int i) { return boundDofT(std,i); }
  int stdDofNb(int i, int j) { return boundDofT(std,i,j); }

  // mpc helper functions
  int mpcNb(int i) { return boundDofT(mpc,i); }
  int mpcNb(int i, int j) { return boundDofT(mpc,i,j); }

  // wet interface helper functions
  int wetDofNb(int i) { return boundDofT(wet,i); }
  int wetDofNb(int i, int j) { return boundDofT(wet,i,j); }

  // all helper function
  int *allBoundDofs() {  return (SharedDOFs[all]->numConnect() > 0) ? (*SharedDOFs[all])[0] : 0; }
  int *allBoundDofs(int i) { return (*SharedDOFs[all])[i]; }
  int totalInterfSize() { return SharedDOFs[all]->numConnect(); }
};

template<class Scalar>
SComm::SComm(int nN, int *subIds, GenSubDomain<Scalar> **subs, int *ids, Connectivity *con)
{
  glSubToLocal = 0;
  sharedDOFs = 0;
  numNeighb = nN;
  subNums = subIds;
  neighb = new BaseSub*[numNeighb];
  for(int i =0; i<numNeighb; ++i)
    neighb[i] = dynamic_cast<BaseSub *>(subs[i]);
  remoteId = ids;
  sharedNodes = con;
  exchangeData = new void*[numNeighb];
  for(int i=0; i<numNeighb; ++i) exchangeData[i] = 0;  // PJSA
  numEdgeNeighb = 0;  
  isEdgeNeighb = 0; 
  sharedDOFsPlus = 0;

  // type specific lists
  numDofType = 8;
  NumNeighb = 0;
  SubNums = 0;
  SharedDOFs = 0;
  TypeMap = 0;
}
#endif
