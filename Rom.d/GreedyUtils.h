#ifndef ROM_GREEDYUTILS_H
#define ROM_GREEDYUTILS_H

#include "SimpleBuffer.h"

#include <algorithm>
#include <iterator>
#include <functional>
#include <cassert>

// Helper class: Allows the ordering of candidates
template <typename V, typename W>
struct SampleCandidate {
  typedef V ValueType;
  typedef W WeightType;

  V value;
  W weight;

  SampleCandidate() {}

  explicit SampleCandidate(V v, W w = W()) :
    value(v), weight(w)
  {}

  bool operator<(const SampleCandidate &other) const {
    return std::less<W>()(weight, other.weight);
  }
  
  bool operator>(const SampleCandidate &other) const {
    return std::greater<W>()(weight, other.weight);
  }
};

// Helper class: Generates a enumerated sequence of SampleCandidates
template <typename V, typename W>
class SampleCandidateGenerator {
public:
  typedef V ValueType;
  typedef W WeightType;
  typedef SampleCandidate<V, W> GeneratedType;

  GeneratedType operator()() { return GeneratedType(value_++); }

  const V &nextValue() const { return value_; }

  explicit SampleCandidateGenerator(V start = V()) :
    value_(start)
  {}

private:
  V value_;
};

// Determine the indices with the highest sum of squares across all vectors
template <typename VecFwdIt, typename IdxOutIt>
IdxOutIt
max_magnitude_indices(VecFwdIt first, VecFwdIt last, IdxOutIt result, int maxResultCount = 1) {
  // Exit early when trivial case is detected
  if (first == last || maxResultCount <= 0) {
    return result;
  }
  
  typedef typename std::iterator_traits<VecFwdIt>::value_type VecType;
  typedef typename VecType::DataType Scalar;
 
  // Determine the number of candidates and initialize
  const int indexCount = first->size();
  maxResultCount = std::min(indexCount, maxResultCount);

  typedef SampleCandidate<int, Scalar> Candidate;
  SimpleBuffer<Candidate> candidates(indexCount);
  std::generate_n(candidates.array(), indexCount, SampleCandidateGenerator<int, Scalar>());

  // Compute the magnitude for each candidate
  for (VecFwdIt it = first; it != last; ++it) {
    const VecType &vec = *it;
    assert(vec.size() == indexCount);
    for (int i = 0; i != indexCount; ++i) {
      const Scalar entry = vec[i];
      candidates[i].weight += entry * entry;
    }
  }

  // Find the candidates with the highest magnitude
  Candidate *lastSelected = candidates.array() + maxResultCount;
  std::partial_sort(candidates.array(),
                    lastSelected,
                    candidates.array() + indexCount,
                    std::greater<Candidate>());

  // Output the selected indices (Ordered by decreasing magnitude)
  for (const Candidate *c = candidates.array(); c != lastSelected; ++c) {
    *result++ = c->value;
  }

  return result;
}

#endif /* ROM_GREEDYUTILS_H */
