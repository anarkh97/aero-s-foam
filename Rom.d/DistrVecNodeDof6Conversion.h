#ifndef ROM_DISTRVECNODEDOF6CONVERSION_H
#define ROM_DISTRVECNODEDOF6CONVERSION_H

#include "RestrictedVecNodeDof6Conversion.h"
#include "VecNodeDof6Conversion.h"
#include "DistrDomainUtils.h"

#include <Driver.d/SubDomain.h> 

#include <vector>
#include <cassert>

namespace Rom {

class DistrVecNodeDof6Conversion {
public:
  int subDomainCount() const { return subDomains_.size(); }

  template <typename NodeDof6Type, typename VecType>
  const NodeDof6Type &paddedNodeDof6(const VecType &origin, NodeDof6Type &target) const;

  template <typename NodeDof6Type, typename VecType>
  const VecType &vector(const NodeDof6Type &origin, VecType &target) const;
  
  template <typename NodeDof6Type, typename VecType>
  const VecType &paddedMasterVector(const NodeDof6Type &origin, VecType &target) const;

  template <typename SubDomPtrFwdIt>
  DistrVecNodeDof6Conversion(SubDomPtrFwdIt first, SubDomPtrFwdIt last);

  ~DistrVecNodeDof6Conversion();

private:
  typedef std::vector<const SubDomain *> SubDomContainer;
  typedef std::vector<const VecNodeDof6Conversion *> ConversionContainer;
  typedef std::vector<const RestrictedVecNodeDof6Conversion *> RestrictedConversionContainer;
  SubDomContainer subDomains_;
  ConversionContainer subConversions_;
  RestrictedConversionContainer subRestrictedConversions_;

  // Disallow copy and assignment
  DistrVecNodeDof6Conversion(const DistrVecNodeDof6Conversion &);
  DistrVecNodeDof6Conversion &operator=(const DistrVecNodeDof6Conversion &);
};

template <typename SubDomPtrFwdIt>
DistrVecNodeDof6Conversion::DistrVecNodeDof6Conversion(SubDomPtrFwdIt first, SubDomPtrFwdIt last) :
  subDomains_(first, last)
{
  for (SubDomPtrFwdIt it = first; it != last; ++it) {
    SubDomain * const s = *it;

    subConversions_.push_back(new VecNodeDof6Conversion(*s->getCDSA()));

    std::vector<bool> masterFlags;
    master_node_flags(*s, std::back_inserter(masterFlags));
    subRestrictedConversions_.push_back(new RestrictedVecNodeDof6Conversion(*s->getCDSA(),
                                                                            masterFlags.begin(),
                                                                            masterFlags.end()));
  }
}

template <typename NodeDof6Type>
class SubNodeDof6Adapter {
public:
  typedef double Scalar; 

  explicit SubNodeDof6Adapter(const NodeDof6Type &nodeDof6, const SubDomain &subDom) :
    nodeDof6_(nodeDof6), subDomain_(subDom)
  {}

  const Scalar *operator[](int locIdx) const {
    const int globIdx = const_cast<SubDomain &>(subDomain_).localToGlobal(locIdx);
    const Scalar *result = nodeDof6_[globIdx];
    assert(result);
    return result;
  }

  Scalar *operator[](int locIdx) {
    const SubNodeDof6Adapter &self = *this;
    return const_cast<Scalar *>(self[locIdx]);
  }

private:
  const NodeDof6Type &nodeDof6_;
  const SubDomain &subDomain_;
};

template <typename NodeDof6Type, typename VecType>
const NodeDof6Type &
DistrVecNodeDof6Conversion::paddedNodeDof6(const VecType &origin, NodeDof6Type &target) const {
  for (int iSub = 0; iSub < subDomainCount(); ++iSub) {
    SubNodeDof6Adapter<NodeDof6Type> subNodeDof(target, *subDomains_[iSub]);
    const GenStackVector<double> subVector(const_cast<VecType &>(origin).subData(iSub),
                                           const_cast<VecType &>(origin).subLen(iSub));
    subRestrictedConversions_[iSub]->paddedNodeDof6(subVector, subNodeDof);
  }

  return target;
}

template <typename NodeDof6Type, typename VecType>
const VecType &
DistrVecNodeDof6Conversion::vector(const NodeDof6Type &origin, VecType &target) const {
  for (int iSub = 0; iSub < subDomainCount(); ++iSub) {
    const SubNodeDof6Adapter<NodeDof6Type> subNodeDof(origin, *subDomains_[iSub]);
    GenStackVector<double> subVector(target.subData(iSub), target.subLen(iSub));
    subConversions_[iSub]->vector(subNodeDof, subVector);
  }

  return target;
};

template <typename NodeDof6Type, typename VecType>
const VecType &
DistrVecNodeDof6Conversion::paddedMasterVector(const NodeDof6Type &origin, VecType &target) const {
  for (int iSub = 0; iSub < subDomainCount(); ++iSub) {
    const SubNodeDof6Adapter<NodeDof6Type> subNodeDof(origin, *subDomains_[iSub]);
    GenStackVector<double> subVector(target.subData(iSub), target.subLen(iSub));
    subRestrictedConversions_[iSub]->paddedVector(subNodeDof, subVector);
  }

  return target;
};

} // end namespace Rom

#endif /* ROM_DISTRVECNODEDOF6CONVERSION_H */
