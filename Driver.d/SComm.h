#ifndef _SCOMM_H_
#define _SCOMM_H_

#include <Utils.d/Connectivity.h>

class Connectivity;

/// \brief Subdomain Communication class.
class SComm
{
public:
	int locSubNum = -1;     //!< when running in distributed mode
	int *glSubToLocal = nullptr; //!< mapping for distributed mode

	int *remoteId = nullptr;     //!< id of this subdomain in the corresponding neighbor
	int numNeighb = 0;     //!< Number of neighbors with shared nodes
	int *subNums = nullptr;      //!< identification numbers of neighbors with shared nodes
	Connectivity *sharedNodes = nullptr; //!< nodes shared with the neighbors (wet and dry but no virtual)
	std::vector<void *> exchangeData;

	int numEdgeNeighb = 0; //!< number of neighbors that have at least one non-virtual node on the interface
	bool *isEdgeNeighb = nullptr;

	Connectivity *sharedDOFs = nullptr;  //<! DOFs # shared (all of types 0, 1 and 2)
	Connectivity *sharedDOFsPlus = nullptr; //<! also includes corner dofs

	SComm(int numNeighb, int *subNums, int *remoteId, Connectivity *);
	~SComm();

	void setExchangeData(int iSub, void *data);
	void *getExchangeData(int iSub);
	void setEdgeNeighb(int _numEdgeNeighb, bool *_isEdgeNeighb);

	int neighborIndex(int i) const { return subNums[i]; }

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
	int *NumNeighb = nullptr;
	int **SubNums = nullptr;
	int **TypeMap = nullptr;  // maps from type-specific interface to general interface
	Connectivity **SharedDOFs = nullptr;

public:
	void setNumDofType(int _numDofType) { numDofType = _numDofType; }
	// function to make type-specific lists from *sharedDOFs combined list using boundDofFlag
	//void makeTypeSpecificLists(int *boundDofFlag);
	void print(DofType t);
	Connectivity *getTypeSpecificList(DofType type) { return SharedDOFs[type]; }

	// function to set any one individual type-specific list
	// if this function is used to set the individual types then mergeTypeSpecificLists() must be called
	void setTypeSpecificList(DofType type, int *_subNums, Connectivity *_sharedDOFs);
	void deleteTypeSpecificList(DofType type);
	int* mergeTypeSpecificLists();
	void setTypeMap(DofType t, int *map);

	// functions to access any of the type-specific lists
	int numT(DofType type) const { return SharedDOFs[type]->csize(); }
	int neighbT(DofType type, int iNeighb) const { return SubNums[type][iNeighb]; }
	int mapT(DofType type, int iDof) const { return TypeMap[type][iDof]; }
	int mapT(DofType type, int iNeighb, int jDof) const { return TypeMap[type][SharedDOFs[type]->offset(iNeighb)+jDof]; }
	int lenT(DofType type) const { return SharedDOFs[type]->numConnect(); }
	int lenT(DofType type, int iNeighb) const { return SharedDOFs[type]->num(iNeighb); }
	int boundDofT(DofType type, int iDof) const { return (*SharedDOFs[type])[0][iDof]; }
	int boundDofT(DofType type, int iNeighb, int jDof) const { return (*SharedDOFs[type])[iNeighb][jDof]; }
	int offsetT(DofType type, int iNeighb) const { return SharedDOFs[type]->offset(iNeighb); }
	int offsetT(DofType type, int iNeighb, int jDof) const { return SharedDOFs[type]->offset(iNeighb)+jDof; }
	const int* boundDofsT(DofType type) const { return (*SharedDOFs[type])[0]; }
	const int* boundDofsT(DofType type, int iNeighb) const { return (*SharedDOFs[type])[iNeighb]; }
	const int* neighbsT(DofType type) const { return SubNums[type]; }

	// standard (boundary) dofs helper functions
	int stdDofNb(int i) const { return boundDofT(std,i); }
	int stdDofNb(int i, int j) const { return boundDofT(std,i,j); }

	// mpc helper functions
	int mpcNb(int i) const { return boundDofT(mpc,i); }
	int mpcNb(int i, int j) const { return boundDofT(mpc,i,j); }

	// wet interface helper functions
	int wetDofNb(int i) const { return boundDofT(wet,i); }
	int wetDofNb(int i, int j) const { return boundDofT(wet,i,j); }

	// all helper function
	const int *allBoundDofs() const {  return (SharedDOFs[all]->numConnect() > 0) ? (*SharedDOFs[all])[0] : 0; }
	const int *allBoundDofs(int i) const { return (*SharedDOFs[all])[i]; }
	int totalInterfSize() const { return SharedDOFs[all]->numConnect(); }
};

#endif
