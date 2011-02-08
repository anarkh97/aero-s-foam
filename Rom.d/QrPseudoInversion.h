#ifndef ROM_QRPSEUDOINVERSION_H
#define ROM_QRPSEUDOINVERSION_H

#include "VecBasis.h"

#include "SimpleBuffer.h"

class QrPseudoInversion {
public:
  // Overwrites the provided full-column-rank basis
  // with the transpose of its pseudo-inverse.
  const VecBasis &operator()(VecBasis &basis) const;

  QrPseudoInversion() {}

private:
  typedef SimpleBuffer<double> ScalarBuffer;

  // Disallow copy and assignment
  QrPseudoInversion(const QrPseudoInversion &);
  QrPseudoInversion &operator=(const QrPseudoInversion &);
};

#endif /* ROM_QRPSEUDOINVERSION_H */
