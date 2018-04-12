//
// Created by Michel Lesoinne on 3/28/18.
//

#ifndef FEM_SHAREDNODES_H
#define FEM_SHAREDNODES_H

#include <FetiLib/VectorReference.h>
#include <FetiLib/Types.h>

namespace FetiLib {
/// Information about the nodes that are shared with other subdomains.
class SharedNodes {
public:
/// Create the shared node information.
	SharedNodes(int numNeighbors,
	            VectorReference<const global_subdomain_index> neighbors, VectorReference<const size_t> neighborStart,
	            VectorReference<const local_node_index> sharedNodeIndices);

};

}// namespace FetiLib

#endif //FEM_SHAREDNODES_H
