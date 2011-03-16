#ifndef ROM_SAMPLINGERROREVALUATION_H
#define ROM_SAMPLINGERROREVALUATION_H

#include "LeastSquaresSolver.h"

#include <algorithm>
#include <iterator>
#include <stdexcept>

namespace Rom {

class SamplingErrorEvaluation {
public:
  SamplingErrorEvaluation() {}

  template <typename VecRandomIt, typename IndexRandomIt>
  const typename std::iterator_traits<VecRandomIt>::value_type &
  operator()(VecRandomIt firstVec, VecRandomIt lastVec,
             IndexRandomIt firstSampleIdx, IndexRandomIt lastSampleIdx, 
             typename std::iterator_traits<VecRandomIt>::value_type &result);

private:
  GenLeastSquaresSolver<double> solver_;
};

template <typename VecRandomIt, typename IndexRandomIt>
const typename std::iterator_traits<VecRandomIt>::value_type &
SamplingErrorEvaluation::operator()(VecRandomIt firstVec, VecRandomIt lastVec,
                                    IndexRandomIt firstSampleIdx, IndexRandomIt lastSampleIdx,  
                                    typename std::iterator_traits<VecRandomIt>::value_type &result) {
  if (firstVec == lastVec) {
    throw std::domain_error("Must have at least one vector");
  }

  // Initialize with the value of the last vector
  result = lastVec[-1];
  
  typedef typename std::iterator_traits<VecRandomIt>::difference_type vector_count_type;
  vector_count_type vectorCount = std::distance(firstVec, lastVec);
 
  // No minimization required when only one vector 
  if (vectorCount > 1) {
    typedef typename std::iterator_traits<IndexRandomIt>::difference_type sample_count_type;
    sample_count_type sampleCount = std::distance(firstSampleIdx, lastSampleIdx);
   
    // Determine a linear combination of all vectors except the last
    vector_count_type coefCount = vectorCount - 1; 

    // Must solve an overdetermined minimization problem:
    // One unknown for each sought-after coefficient
    // One equation for each sample index
    if (coefCount > sampleCount) {
      throw std::domain_error("Must have #vectors <= #samples + 1");
    }

    solver_.problemSizeIs(sampleCount, coefCount);

    // Assemble the sampled matrix (restriction of all vectors except the last)
    for (int col = 0; col != solver_.unknownCount(); ++col) {
      typedef typename std::iterator_traits<VecRandomIt>::value_type VecType;
      const VecType &vec = firstVec[col];
      for (int row = 0; row != solver_.equationCount(); ++row) {
        solver_.matrixEntry(row, col) = vec[firstSampleIdx[row]];
      }
    }

    // Assemble the sampled rhs (restriction of the last vector)
    for (int row = 0; row != solver_.equationCount(); ++row) {
      solver_.rhsEntry(row) = result[firstSampleIdx[row]];
    }

    // Solve the minimization problem to determine the coefficients
    solver_.statusIs(LeastSquares::SOLVED);

    // Substract the linear combination of all vectors except the last
    for (int iVec = 0; iVec != solver_.unknownCount(); ++iVec) {
      result.linAdd(-solver_.rhsEntry(iVec), firstVec[iVec]);
    }
    
    solver_.statusIs(LeastSquares::READY);
  }

  return result;
}

} /* end namespace Rom */

#endif /* ROM_SAMPLINGERROREVALUATION_H */
