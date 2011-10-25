#ifndef ROM_RENUMBERINGUTILS_H
#define ROM_RENUMBERINGUTILS_H

#include <utility>
#include <map>
#include <vector>

class Connectivity;

namespace Rom {

class SampledMeshRenumbering {
public:
  typedef std::vector<int> Restriction;
  typedef std::map<int, int> Extension;

  const Restriction &sampleNodeIds() const { return sampleNodeIds_; }
  const Restriction &reducedNodeIds() const { return reducedNodeIds_; }
  const Restriction &reducedSampleNodeIds() const { return reducedSampleNodeIds_; }
  const Extension &nodeRenumbering() const { return nodeRenumbering_; }

  const Restriction &reducedElemIds() const { return reducedElemIds_; }
  const Extension &elemRenumbering() const { return elemRenumbering_; }

  template <typename IdxInIt>
  SampledMeshRenumbering(IdxInIt firstSample, IdxInIt lastSample,
                         const Connectivity &nodeToNode,
                         const Connectivity &nodeToElem);

private:
  Restriction sampleNodeIds_;
  Restriction reducedNodeIds_;
  Restriction reducedSampleNodeIds_;
  Extension nodeRenumbering_;
  
  Restriction reducedElemIds_;
  Extension elemRenumbering_;

  void init(const Connectivity &, const Connectivity &);
};

template <typename IdxInIt>
SampledMeshRenumbering::SampledMeshRenumbering(IdxInIt firstSample,
                                               IdxInIt lastSample,
                                               const Connectivity &nodeToNode,
                                               const Connectivity &nodeToElem) :
  sampleNodeIds_(firstSample, lastSample)
{
  init(nodeToNode, nodeToElem);
}

// Precondition: The range [first, last) must not contain duplicate elements
template <typename InputIterator, typename OutputIterator>
OutputIterator
inverse_numbering(InputIterator first, InputIterator last,
                  OutputIterator result) {
  int index = 0;

  for (InputIterator it = first; it != last; ++it) {
    *result++ = std::pair<const int, int>(*it, index++);
  }

  return result;
}

template <typename InputIterator, typename OutputIterator>
OutputIterator
renumber(const std::map<int, int> &renum, InputIterator first, InputIterator last, OutputIterator result) {
  for (InputIterator it = first; it != last; ++it) {
    std::map<int, int>::const_iterator r = renum.find(*it);
    if (r != renum.end()) {
      *result++ = r->second;
    }
  }

  return result;
}

} /* end namespace Rom */

#endif /* ROM_RENUMBERINGUTILS_H */
