#ifndef ROM_NODALRESTRICTIONMAPPING_H
#define ROM_NODALRESTRICTIONMAPPING_H

class DofSetArray;

#include <vector>

#include <cassert>

class NodalRestrictionMapping {
public:
  typedef int InfoType;

  const InfoType &originInfo()     const { return originInfo_;     }
  const InfoType &restrictedInfo() const { return restrictedInfo_; }

  template <typename VecType>
  const VecType &restriction(const VecType &origin, VecType &target) const;

  template <typename VecType>
  typename VecType::DataType dotProduct(const VecType &originVec,
                                        const VecType &restrictedVec) const;

  template <typename InputIterator>
  NodalRestrictionMapping(const DofSetArray &,
                          InputIterator sampleNodesBegin,
                          InputIterator sampleNodesEnd);

private:
  static InfoType extractOriginalInfo(const DofSetArray &);

  void addSampleNode(int, const DofSetArray &);

  InfoType originInfo_;
  InfoType restrictedInfo_;

  typedef int IndexType;
  std::vector<IndexType> originIndex_;

  // Disallow copy and assignment
  NodalRestrictionMapping(const NodalRestrictionMapping &);
  NodalRestrictionMapping &operator=(const NodalRestrictionMapping &);
};

template <typename InputIterator>
NodalRestrictionMapping::NodalRestrictionMapping(const DofSetArray &dsa,
                                                 InputIterator snBegin,
                                                 InputIterator snEnd) :
  originInfo_(extractOriginalInfo(dsa))
{
  for (InputIterator itNode = snBegin; itNode != snEnd; ++itNode) {
    addSampleNode(*itNode, dsa);
  }
  restrictedInfo_ = originIndex_.size();
}

template <typename VecType>
const VecType &
NodalRestrictionMapping::restriction(const VecType &origin, VecType &target) const {
  assert(origin.info() == originInfo());  
  assert(target.info() == restrictedInfo());  

  typedef std::vector<IndexType>::const_iterator Iterator;
  int targetIdx = 0;
  for (Iterator it = originIndex_.begin(); it != originIndex_.end(); ++it) {
    target[targetIdx++] = origin[*it];
  }

  return target;
}

template <typename VecType>
typename VecType::DataType
NodalRestrictionMapping::dotProduct(const VecType &originVec, const VecType &restrictedVec) const {
  assert(originVec.info() == originInfo());
  assert(restrictedVec.info() == restrictedInfo());

  typename VecType::DataType result = typename VecType::DataType();

  typedef std::vector<IndexType>::const_iterator Iterator;
  int restrictedIdx = 0;
  for (Iterator it = originIndex_.begin(); it != originIndex_.end(); ++it) {
    result += restrictedVec[restrictedIdx++] * originVec[*it];
  }

  return result;
}

#endif /* ROM_NODALRESTRICTIONMAPPING_H */
