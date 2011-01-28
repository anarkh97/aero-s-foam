#ifndef ROM_BASISORTHOGONALIZATION_H
#define ROM_BASISORTHOGONALIZATION_H

#include "SvdOrthogonalization.h"

#include <algorithm>
#include <utility>

class BasisOrthogonalization {
public:
  template <typename InputRange, typename OutputRange>
  OutputRange &basisNew(InputRange &input, OutputRange &output);

  BasisOrthogonalization() {}

private:
  SvdOrthogonalization solver_;

  // Disallow copy and assignment
  BasisOrthogonalization(const BasisOrthogonalization &);
  BasisOrthogonalization &operator=(const BasisOrthogonalization &);
};

template <typename InputRange, typename OutputRange>
OutputRange &
BasisOrthogonalization::basisNew(InputRange &input, OutputRange &output) {
  solver_.matrixSizeIs(input.vectorSize(), input.size());

  int iCol = 0;
  while (input) {
    input >> solver_.matrixCol(iCol++);
  }

  solver_.solve();

  for (int iVec = 0; iVec < solver_.singularValueCount(); ++iVec) {
    output << std::make_pair(solver_.singularValue(iVec), solver_.matrixCol(iVec));
  }

  return output;
}

#endif /* ROM_BASISORTHOGONALIZATION_H */
