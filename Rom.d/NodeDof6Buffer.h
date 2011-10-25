#ifndef ROM_NODEDOF6BUFFER_H
#define ROM_NODEDOF6BUFFER_H

#include "SimpleBuffer.h"

#include <cstddef>

namespace Rom {

class NodeDof6Buffer {
public:
  explicit NodeDof6Buffer(size_t nodeCount = 0) :
    nodeCount_(nodeCount),
    buffer_(6 * nodeCount)
  {}

  size_t size() const { return nodeCount_; }
  void sizeIs(size_t nodeCount) { buffer_.sizeIs(6 * nodeCount); nodeCount_ = nodeCount; }

  const double *operator[](size_t iNode) const { return buffer_.array() + (6 * iNode); }
  double *operator[](size_t iNode) { return buffer_.array() + (6 * iNode); }

  double *array() { return buffer_.array(); }
  const double *array() const { return buffer_.array(); }

private:
  size_t nodeCount_;
  SimpleBuffer<double> buffer_;

  // Disallow copy and assignment
  NodeDof6Buffer(const NodeDof6Buffer &);
  NodeDof6Buffer &operator=(const NodeDof6Buffer &);
};

} /* end namespace Rom */

#endif /* ROM_NODEDOF6BUFFER_H */
