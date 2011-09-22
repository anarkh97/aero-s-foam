#ifndef ROM_DISTRDOMAINUTILS_H
#define ROM_DISTRDOMAINUTILS_H

#include <Driver.d/SubDomain.h>

#include <vector>
#include <algorithm>

namespace Rom {

// Master nodes of the subdomain (in local indexing)
template <typename BoolOutIt>
BoolOutIt
master_node_flags(const SubDomain &subDom, BoolOutIt result) {
  SubDomain &sd = const_cast<SubDomain &>(subDom); // Fix constness
  const int mySubId = sd.subNum();

  const int nodeCount = sd.numNode();
  std::vector<bool> work(nodeCount, true); // vector<bool> specialization acceptable here

  Connectivity *sharedNodes = sd.getSComm()->sharedNodes;
  
  const int neighborCount = sd.getSComm()->numNeighb;
  for (int iNeighbor = 0; iNeighbor < neighborCount; ++iNeighbor) {
    const int neighborId = sd.getSComm()->subNums[iNeighbor];
    if (mySubId > neighborId) {
      const int sharedNodeCount = sharedNodes->num(iNeighbor);
      for (int iNode = 0; iNode < sharedNodeCount; ++iNode) {
        const int nodeId = (*sharedNodes)[iNeighbor][iNode];
        work[nodeId] = false;
      }
    }
  }

  return std::copy(work.begin(), work.end(), result);
}

} /* end namespace Rom */

#endif /* ROM_DISTRDOMAINUTILS_H */
