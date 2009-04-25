#ifndef PITA_DYNAMSTATEBASIS_H
#define PITA_DYNAMSTATEBASIS_H

#include "Fwk.h"
#include "DynamState.h"

namespace Pita {

class DynamStateBasis : public Fwk::PtrInterface<DynamStateBasis> {
public:
  typedef Fwk::Ptr<DynamStateBasis> Ptr;
  typedef Fwk::Ptr<const DynamStateBasis> PtrConst;

  size_t vectorSize() const { return vectorSize_; }
  virtual size_t stateCount() const = 0;
  virtual DynamState state(size_t index) const = 0; // Unsafe

  //class IteratorConst;
  
protected:
  explicit DynamStateBasis(size_t vectorSize) : vectorSize_(vectorSize) {}
  
  DynamStateBasis(const DynamStateBasis &); // No implementation
  DynamStateBasis & operator=(const DynamStateBasis &); // No implementation
  
private:
  size_t vectorSize_;
};

/*class DynamStateBasis::IteratorConst {
public:
  virtual const DynamState & operator*() const = 0;
  virtual IteratorConst & operator++() = 0;
  virtual IteratorConst operator++(int) = 0;
  virtual operator bool() = 0;
};*/

}

#endif /* PITA_DYNAMSTATEBASIS_H */
