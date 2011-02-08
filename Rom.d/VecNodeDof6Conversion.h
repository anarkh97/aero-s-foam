#ifndef ROM_VECNODEDOF6CONVERSION_H
#define ROM_VECNODEDOF6CONVERSION_H

class DofSetArray;

#include <vector>
#include <cassert>

struct NodeDof {
  typedef int DofType;

  int nodeRank;
  DofType dofId;

  NodeDof(int n, DofType d) :
    nodeRank(n),
    dofId(d)
  {}

  NodeDof() {}
};

class VecNodeDof6Conversion {
public:
  int nodeCount() const { return nodeCount_; }
  int vectorSize() const { return vectorSize_; }

  template <typename NodeDofs6Type, typename VecType>
  const NodeDofs6Type &nodeDof6(const VecType &origin, NodeDofs6Type &target) const;

  template <typename NodeDofs6Type, typename VecType>
  const VecType &vector(const NodeDofs6Type &origin, VecType &target) const;

  NodeDof nodeDof(int vecLoc) const;

  template <typename IndexOut>
  IndexOut locations(int nodeRank, IndexOut result) const;

  explicit VecNodeDof6Conversion(const DofSetArray &);

  ~VecNodeDof6Conversion();

private:
  int nodeCount_;
  int vectorSize_;

  std::vector<NodeDof> locationId_;
  
  typedef int (*DofLocation)[6];
  DofLocation dofLocation_;

  // Disallow copy and assignment
  VecNodeDof6Conversion(const VecNodeDof6Conversion &);
  VecNodeDof6Conversion &operator=(const VecNodeDof6Conversion &);
};

template <typename NodeDofs6Type, typename VecType>
const NodeDofs6Type &
VecNodeDof6Conversion::nodeDof6(const VecType &origin, NodeDofs6Type &target) const {
  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    for (int iDof = 0; iDof < 6; ++iDof) {
      const int loc = dofLocation_[iNode][iDof];
      target[iNode][iDof] = (loc >= 0) ? origin[loc] : 0.0;
    }
  }

  return target;
}

template <typename NodeDofs6Type, typename VecType>
const VecType &
VecNodeDof6Conversion::vector(const NodeDofs6Type &origin, VecType &target) const {
  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    for (int iDof = 0; iDof < 6; ++iDof) {
      const int loc = dofLocation_[iNode][iDof];
      if (loc >= 0) {
        target[loc] = origin[iNode][iDof];
      }
    }
  }

  return target;
}

template <typename IndexOut>
IndexOut
VecNodeDof6Conversion::locations(int nodeRank, IndexOut result) const {
  assert(nodeRank >= 0 && nodeRank < nodeCount());
  const int *nodeLoc = dofLocation_[nodeRank];

  for (int iDof = 0; iDof < 6; ++iDof) {
    const int loc = nodeLoc[iDof];
    if (loc >= 0) {
      *result++ = loc;
    }
  }

  return result;
}

#endif /* ROM_VECNODEDOF6CONVERSION_H */
