#ifndef ROM_BASISORTHOGONALIZATION_H
#define ROM_BASISORTHOGONALIZATION_H

#include "SvdOrthogonalization.h"

#include <algorithm>
#include <utility>

class BasisOrthogonalization {
public:
  template <typename InputRange, typename OutputIterator>
  OutputIterator basisNew(InputRange &input, OutputIterator output);

  BasisOrthogonalization() {}

private:
  SvdOrthogonalization solver_;

  // Disallow copy and assignment
  BasisOrthogonalization(const BasisOrthogonalization &);
  BasisOrthogonalization &operator=(const BasisOrthogonalization &);
};

template <typename InputRange, typename OutputIterator>
OutputIterator
BasisOrthogonalization::basisNew(InputRange &input, OutputIterator output) {
  solver_.matrixSizeIs(input.vectorSize(), input.size());

  typedef typename InputRange::const_iterator InputIterator;
  const InputIterator inputEnd = input.end();
  
  int iCol = 0;
  for (InputIterator it = input.begin(); it != inputEnd; ++it) {
    (*it)(solver_.matrixCol(iCol));
    ++iCol;
  }

  solver_.solve();

  OutputIterator it = output;
  for (int iVec = 0; iVec < solver_.singularValueCount(); ++iVec) {
    *it = std::make_pair(solver_.singularValue(iVec), solver_.matrixCol(iVec));
    ++it;
  }

  return it;
}

#endif /* ROM_BASISORTHOGONALIZATION_H */
