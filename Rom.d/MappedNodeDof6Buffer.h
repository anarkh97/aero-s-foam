#ifndef ROM_DISTRNODEDOF6BUFFER_H
#define ROM_DISTRNODEDOF6BUFFER_H

#include "NodeDof6Buffer.h"

#include <vector>
#include <map>

namespace Rom {

class MappedNodeDof6Buffer {
public:
  int nodeCount() const { return nodeIndices_.size(); }

  typedef std::vector<int>::const_iterator NodeItConst;
  NodeItConst nodeIndexBegin() const { return nodeIndices_.begin(); }
  NodeItConst nodeIndexEnd()   const { return nodeIndices_.end(); }

  const double *operator[](int iNode) const;
  double *operator[](int iNode) {
    const MappedNodeDof6Buffer &self = *this;
    return const_cast<double *>(self[iNode]);
  }

  double *array() { return buffer_.array(); }
  const double *array() const { return buffer_.array(); }

  const NodeDof6Buffer &underlyingBuffer() const { return buffer_; }
  NodeDof6Buffer &underlyingBuffer() { return buffer_; }

  // Range [first, last) should not have duplicated elements
  template <typename IdxInIt>
  MappedNodeDof6Buffer(IdxInIt first, IdxInIt last);

private:
  void initialize();

  std::map<int, int> underlyingNodeIndices_;
  std::vector<int> nodeIndices_;
  NodeDof6Buffer buffer_;
};

template <typename IdxInIt>
MappedNodeDof6Buffer::MappedNodeDof6Buffer(IdxInIt first, IdxInIt last) :
  underlyingNodeIndices_(),
  nodeIndices_(first, last),
  buffer_()
{
  initialize();
}

} /* end namespace Rom */

#endif /* ROM_DISTRNODEDOF6BUFFER_H */
