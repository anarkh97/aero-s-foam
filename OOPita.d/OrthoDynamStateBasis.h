#ifndef PITA_ORTHODYNAMSTATEBASIS_H
#define PITA_ORTHODYNAMSTATEBASIS_H

#include "Fwk.h"

#include "DynamState.h"
#include "DynamStateOps.h"

#include "DynamStatePlainBasis.h"

namespace Pita {

class OrthoDynamStateBasis {
public:
  size_t vectorSize() const { return vectorSize_; }

  const DynamOps * metric() const { return metric_.ptr(); }

  const DynamStateBasis * orthoBasis() const { return orthoBasis_.ptr(); }
  const DynamStateBasis * dualBasis() const { return dualBasis_.ptr(); }

  double tolerance() const { return tolerance_; }

  // Direct access
  //DynamStatePlainBasis * orthoBasis() { return orthoBasis_.ptr(); }
  //DynamStatePlainBasis * dualBasis() { return dualBasis_.ptr(); }

  // Pass-by-value
  void originalBasisIs(const DynamStateBasis & ob);
  void originalBasisInc(const DynamStateBasis & ob);

  explicit OrthoDynamStateBasis(size_t vectorSize, const DynamOps * metric, double tolerance);

protected:
  void checkSize(const DynamStateBasis & basis);
  void addBasis(const DynamStateBasis & originalBasis);

private:
  size_t vectorSize_;
  DynamOps::PtrConst metric_;
  DynamStatePlainBasis::Ptr orthoBasis_;
  DynamStatePlainBasis::Ptr dualBasis_;
  double tolerance_;
};

} /* end namespace Pita */

#endif /* PITA_ORTHODYNAMSTATEBASIS_H */
