#ifndef ROM_GAPPYNONLINDYNAMIC_H
#define ROM_GAPPYNONLINDYNAMIC_H

#include "GappyProjectionSolver.h"
#include "VecBasis.h"
#include "NodalRestrictionMapping.h"
#include "VecNodeDof6Conversion.h"

#include <Problems.d/NonLinDynam.h>
#include <Driver.d/Domain.h>

#include <string>
#include <memory>

namespace Rom {

class GappyNonLinDynamic : public NonLinDynamic {
public:
  explicit GappyNonLinDynamic(Domain *);

  // Required additional pre- and post-processing
  virtual void preProcess();

  GappyProjectionSolver *getSolver(); // Hiding NonLinDynamic::getSolver
  
  // Direct hook in NLDynamSolver (relies on function hiding)
  double getResidualNorm(const Vector &);

private:
  virtual bool factorWhenBuilding() const; // Overriden

  std::auto_ptr<VecNodeDof6Conversion> vecNodeDof6Conversion_;

  std::auto_ptr<NodalRestrictionMapping> restrictionMapping_;  
  VecBasis reducedBasis_;
  VecBasis jacobianProjection_;
  VecBasis residualProjection_;
  
  void fillBasisFromInput(const std::string &, VecBasis &) const;
  void fillRestrictedBasisFromInput(const std::string &, VecBasis &) const;

  // Disallow copy and assignment
  GappyNonLinDynamic(const GappyNonLinDynamic &);
  GappyNonLinDynamic &operator=(const GappyNonLinDynamic &);
};

} /* end namespace Rom */

#endif /* ROM_GAPPYNONLINDYNAMIC_H */
