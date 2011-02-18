#ifndef ROM_GREEDYUTILS_H
#define ROM_GREEDYUTILS_H

#include "VecNodeDof6Conversion.h"
#include "SamplingErrorEvaluation.h"

#include "SimpleBuffer.h"

#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <functional>
#include <cmath>
#include <cassert>

// debug
#include <iostream>

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


struct VectorCountLess {
  template <typename BasisType>
  bool operator()(const BasisType &a,
                  const BasisType &b) const {
    return a.vectorCount() < b.vectorCount();
  }
};

template <typename BasisFwdIt>
int vector_count_min(BasisFwdIt first, BasisFwdIt last) {
  const BasisFwdIt it = std::min_element(first, last, VectorCountLess());
  return it != last ? it->vectorCount() : 0;
}

template <typename BasisFwdIt>
int vector_count_max(BasisFwdIt first, BasisFwdIt last) {
  const BasisFwdIt it = std::max_element(first, last, VectorCountLess());
  return it != last ? it->vectorCount() : 0;
}

template <typename BasisRanIt, typename IndexRanIt, typename IndexOutIt>
IndexOutIt greedy_sampling(BasisRanIt podFirst, BasisRanIt podLast,
                           IndexRanIt sampleFirst, IndexRanIt sampleLast,
                           IndexOutIt result,
                           const VecNodeDof6Conversion &conversion,
                           double aspectRatio) {
  if (podFirst == podLast) {
    throw std::domain_error("Must have at least one pod basis");
  }

  // Greedy iterations parameters
  const int podVectorCountMax = vector_count_max(podFirst, podLast);
  const int targetSampleSize = std::max(static_cast<int>(aspectRatio * podVectorCountMax), podVectorCountMax);
  const int dofsPerNodeEstimate = 6;

  std::vector<int> sampleLocations;

  if (sampleFirst != sampleLast) {
    // Add all the dofs attached to the user-specified nodes
    std::cerr << "User-specified nodes: ";
    for (std::vector<int>::const_iterator it = sampleFirst; it != sampleLast; ++it) {
      std::cerr << *it + 1 << " ";
      *result++ = *it;
      conversion.locations(*it, std::back_inserter(sampleLocations));
    }
    std::cerr << std::endl;
  }
 
  typedef typename std::iterator_traits<BasisRanIt>::value_type BasisType;
  typedef typename BasisType::VecType VecType; 
  std::vector<VecType> reconstructionError(std::distance(podFirst, podLast));
  SamplingErrorEvaluation errorEval;

  // Start greedy iterations
  int handledVectorCount = 0;
  int greedyIter = 0;
  const int podVectorCountMin = vector_count_min(podFirst, podLast);
  while (sampleLocations.size() < targetSampleSize && handledVectorCount < podVectorCountMin) {
    std::cerr << sampleLocations.size() << " sample locations at greedy iteration: " << greedyIter++ << std::endl;

    // Update the reconstruction error
    handledVectorCount++;
    typename std::vector<VecType>::iterator errorIt = reconstructionError.begin();
    for (typename std::vector<BasisType>::const_iterator podIt = podFirst; podIt != podLast; ++podIt) {
      errorEval(podIt->begin(), podIt->begin() + handledVectorCount,
                sampleLocations.begin(), sampleLocations.end(),
                *errorIt++);
    }

    // Zero out the entries corresponding to the already selected indices to avoid selecting them again
    for (typename std::vector<VecType>::iterator vecIt = reconstructionError.begin(); vecIt != reconstructionError.end(); ++vecIt) {
      VecType &vec = *vecIt;
      for (std::vector<int>::const_iterator it = sampleLocations.begin(); it != sampleLocations.end(); ++it) {
        vec[*it] = 0.0;
      }
    }

    // Greedily determine the best entries and find the corresponding nodes
    const int remainingLocations = targetSampleSize - sampleLocations.size();
    const int remainingPodVectors = (podVectorCountMin - handledVectorCount) + 1;
    const double ratio = static_cast<double>(remainingLocations) / (remainingPodVectors * dofsPerNodeEstimate);
    const int locationsToSelect = static_cast<int>(std::ceil(ratio));

    std::vector<int> selectedLocations;
    max_magnitude_indices(reconstructionError.begin(), reconstructionError.end(),
                          std::back_inserter(selectedLocations), locationsToSelect);

    std::set<int> selectedNodes; // To ensure unicity, since several locations can correspond to one node
    for (std::vector<int>::const_iterator it = selectedLocations.begin(); it != selectedLocations.end(); ++it) {
      selectedNodes.insert(conversion.nodeDof(*it).nodeRank);
    }

    // Add all the dofs attached to the newly selected nodes
    std::cerr << "Selected nodes: ";
    for (std::set<int>::const_iterator it = selectedNodes.begin(); it != selectedNodes.end(); ++it) {
      std::cerr << *it + 1 << " ";
      *result++ = *it; // Newly selected nodes are not already selected
      conversion.locations(*it, std::back_inserter(sampleLocations));
    }
    std::cerr << std::endl;
  }

  return result;
}

#endif /* ROM_GREEDYUTILS_H */
