#ifndef ROM_DISTRNODEDOF6BUFFER_H
#define ROM_DISTRNODEDOF6BUFFER_H

#include "NodeDof6Buffer.h"

#include <vector>
#include <map>

namespace Rom {

class DistrNodeDof6Buffer {
public:
  int localNodeCount() const { return globalNodeIndices_.size(); }

  typedef std::vector<int>::const_iterator NodeItConst;
  NodeItConst globalNodeIndexBegin() const { return globalNodeIndices_.begin(); }
  NodeItConst globalNodeIndexEnd()   const { return globalNodeIndices_.end(); }

  const double *operator[](int globalNodeIdx) const;
  double *operator[](int globalNodeIdx) {
    const DistrNodeDof6Buffer &self = *this;
    return const_cast<double *>(self[globalNodeIdx]);
  }

  template <typename IdxInIt>
  DistrNodeDof6Buffer(IdxInIt first, IdxInIt last);

private:
  void initialize();

  std::map<int, int> localNodeIndices_;
  std::vector<int> globalNodeIndices_;
  NodeDof6Buffer buffer_;
};

template <typename IdxInIt>
DistrNodeDof6Buffer::DistrNodeDof6Buffer(IdxInIt first, IdxInIt last) :
  localNodeIndices_(),
  globalNodeIndices_(first, last),
  buffer_()
{
  initialize();
}

} /* end namespace Rom */

#endif /* ROM_DISTRNODEDOF6BUFFER_H */
