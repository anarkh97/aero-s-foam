#ifndef PITA_DYNAMOPS_H
#define PITA_DYNAMOPS_H

#include "Fwk.h"

template <typename> class GenSparseMatrix;

namespace Pita {

class DynamOps : public Fwk::PtrInterface<DynamOps> {
public:
  typedef Fwk::Ptr<DynamOps> Ptr;
  typedef Fwk::Ptr<const DynamOps> PtrConst;

  virtual const GenSparseMatrix<double> * massMatrix() const = 0;
  virtual const GenSparseMatrix<double> * stiffnessMatrix() const = 0;
};
  
} // end namespace Pita

#endif /* PITA_DYNAMOPS_H */
