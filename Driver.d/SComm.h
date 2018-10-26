#ifndef _SCOMM_H_
#define _SCOMM_H_

#include <memory>
#include <Types.h>
#include <Utils.d/Connectivity.h>

class Connectivity;

/// \brief Subdomain Communication class.
class SComm
{
public:
	/** \brief Constructor.
	 *
	 * @param numNeighb Number of neighboring subdomains.
	 * @param subNums Indices of the neighboring subdomains
	 * @param remoteId ID in each neighbor for this subdomain.
	 * @param sharedNodes Nodes shared with each neighbor.
	 */
	SComm(int numNeighb, std::vector<gl_sub_idx> subNums, std::vector<lc_sub_idx> remoteId,
	      std::unique_ptr<Connectivity> sharedNodes);
	~SComm();

	int locSubNum = -1;     //!< when running in distributed mode

	std::vector<lc_sub_idx> remoteId;     //!< id of this subdomain in the corresponding neighbor
	int numNeighb = 0;     //!< Number of neighbors with shared nodes
	std::vector<gl_sub_idx > subNums;      //!< identification numbers of neighbors with shared nodes
	std::unique_ptr<Connectivity> sharedNodes; //!< nodes shared with the neighbors (wet and dry but no virtual)
	std::vector<void *> exchangeData;

	int numEdgeNeighb = 0; //!< number of neighbors that have at least one non-virtual node on the interface
	bool *isEdgeNeighb = nullptr;

	Connectivity *sharedDOFs = nullptr;  //<! DOFs # shared (all of types 0, 1 and 2) (non owning pointer)
	Connectivity *sharedDOFsPlus = nullptr; //<! also includes corner dofs


	void setExchangeData(int iSub, void *data);
	void *getExchangeData(int iSub);
	void setEdgeNeighb(int _numEdgeNeighb, bool *_isEdgeNeighb);

	int neighborIndex(int i) const { return subNums[i]; }

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
	std::vector<int> NumNeighb;
	std::vector<std::vector<gl_sub_idx>> SubNums;
	std::vector<std::vector<int>> TypeMap;  // maps from type-specific interface to general interface
	/// \brief Shared degrees of freedom per type.
	std::vector<Connectivity> sharedDOFsPerType;

public:
	void setNumDofType(int _numDofType) { numDofType = _numDofType; }
	void print(DofType t);
	const Connectivity *getTypeSpecificList(DofType type) const { return &sharedDOFsPerType[type]; }

	/** \brief Set any one individual type-specific list.
	 * \details If this function is used to set the individual types then mergeTypeSpecificLists() must be called */
	void setTypeSpecificList(DofType type, std::vector<gl_sub_idx> _subNums, std::unique_ptr<Connectivity> _sharedDOFs);
	void deleteTypeSpecificList(DofType type);
	/** build combined list of all types 0, 1 and 2 shared dofs:
	 * also make **TypeMap and *boundDofFlag (returned)
	 * update **neighb and *remoteId
	 * resize **exchangeData but i don't think it is necessary to add "virtual nodes"
	 * to sharedNodes list. however, check in code where sharedNodes is used and convert to "std"
	 * NumNeighb and SubNums etc. */
	std::vector<int> mergeTypeSpecificLists();

	// functions to access any of the type-specific lists
	int numT(DofType type) const { return sharedDOFsPerType[type].csize(); }
	int neighbT(DofType type, int iNeighb) const { return SubNums[type][iNeighb]; }
	int mapT(DofType type, int iDof) const { return TypeMap[type][iDof]; }
	int mapT(DofType type, int iNeighb, int jDof) const { return TypeMap[type][sharedDOFsPerType[type].offset(iNeighb)+jDof]; }
	int lenT(DofType type) const { return sharedDOFsPerType[type].numConnect(); }
	int lenT(DofType type, int iNeighb) const { return sharedDOFsPerType[type].num(iNeighb); }
	int boundDofT(DofType type, int iDof) const { return sharedDOFsPerType[type].allTargets()[iDof]; }
	int boundDofT(DofType type, int iNeighb, int jDof) const { return sharedDOFsPerType[type][iNeighb][jDof]; }
	int offsetT(DofType type, int iNeighb) const { return sharedDOFsPerType[type].offset(iNeighb); }
	int offsetT(DofType type, int iNeighb, int jDof) const { return sharedDOFsPerType[type].offset(iNeighb)+jDof; }
	const int* boundDofsT(DofType type) const { return sharedDOFsPerType[type][0].data(); }
	const int* boundDofsT(DofType type, int iNeighb) const { return sharedDOFsPerType[type][iNeighb].data(); }
	const auto &neighbsT(DofType type) const { return SubNums[type]; }

	// standard (boundary) dofs helper functions
	int stdDofNb(int i) const { return boundDofT(std,i); }
	int stdDofNb(int i, int j) const { return boundDofT(std,i,j); }

	// mpc helper functions
	int mpcNb(int i) const { return boundDofT(mpc,i); }
	int mpcNb(int i, int j) const { return boundDofT(mpc,i,j); }

	// wet interface helper functions
	int wetDofNb(int i) const { return boundDofT(wet,i); }
	int wetDofNb(int i, int j) const { return boundDofT(wet,i,j); }

	// TODO make these return spans.
	const int *allBoundDofs() const {  return (sharedDOFsPerType[all].numConnect() > 0) ? sharedDOFsPerType[all][0].data() : nullptr; }
	const int *allBoundDofs(int i) const { return sharedDOFsPerType[all][i].data(); }
	int totalInterfSize() const { return sharedDOFsPerType[all].numConnect(); }
};

#endif
