#ifndef ROM_DISTRNODEDOF6BUFFER_H
#define ROM_DISTRNODEDOF6BUFFER_H

#include "NodeDof6Buffer.h"

#include <vector>
#include <map>

namespace Rom {

template<int DOFS_PER_NODE>
class DistrNodeDofBuffer {
public:
  int localNodeCount() const { return globalNodeIndices_.size(); }

  typedef std::vector<int>::const_iterator NodeItConst;
  NodeItConst globalNodeIndexBegin() const { return globalNodeIndices_.begin(); }
  NodeItConst globalNodeIndexEnd()   const { return globalNodeIndices_.end(); }

  const double *operator[](int globalNodeIdx) const;
  double *operator[](int globalNodeIdx) {
    const DistrNodeDofBuffer &self = *this;
    return const_cast<double *>(self[globalNodeIdx]);
  }

  template <typename IdxInIt>
  DistrNodeDofBuffer(IdxInIt first, IdxInIt last);
  DistrNodeDofBuffer(int size);

private:
  void initialize();

  std::map<int, int> localNodeIndices_;
  std::vector<int> globalNodeIndices_;
  NodeDofBuffer<DOFS_PER_NODE> buffer_;
};

typedef DistrNodeDofBuffer<6> DistrNodeDof6Buffer;
typedef DistrNodeDofBuffer<1> DistrNodeDof1Buffer;

template <int DOFS_PER_NODE>
template <typename IdxInIt>
DistrNodeDofBuffer<DOFS_PER_NODE>::DistrNodeDofBuffer(IdxInIt first, IdxInIt last) :
  localNodeIndices_(),
  globalNodeIndices_(first, last),
  buffer_()
{
  initialize();
}

template <int DOFS_PER_NODE>
DistrNodeDofBuffer<DOFS_PER_NODE>::DistrNodeDofBuffer(int size) :
  localNodeIndices_(),
  globalNodeIndices_(),
  buffer_()
{
  for(int i=0; i<size; ++i) globalNodeIndices_.push_back(i);
  initialize();
}

} /* end namespace Rom */

#endif /* ROM_DISTRNODEDOF6BUFFER_H */
