#ifndef PITA_DYNAMSTATEPLAINBASIS_H
#define PITA_DYNAMSTATEPLAINBASIS_H

#include "DynamStateBasis.h"
#include <deque>

namespace Pita {

class DynamStatePlainBasis : public DynamStateBasis {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStatePlainBasis);

  virtual size_t stateCount() const { return state_.size(); }
  using DynamStateBasis::state;
  DynamState state(size_t index) const { return state_[index]; } // Unsafe

  void stateIs(size_t index, const DynamState & newState) { state_[index] = newState; } // Unsafe

  virtual void lastStateIs(const DynamState & ds);
  virtual void lastStateBasisIs(const DynamStateBasis * dsb);
  virtual void lastStateBasisIs(const DynamStatePlainBasis * dsb);

  void stateBasisDel();
  
  static DynamStatePlainBasis::Ptr New(size_t vectorSize) { return new DynamStatePlainBasis(vectorSize); }
  
protected:
  explicit DynamStatePlainBasis(size_t vectorSize) : DynamStateBasis(vectorSize) {}
 
  void addState(const DynamState & ds);
  void addStateBasis(const DynamStatePlainBasis * dsb);
  
private:
  std::deque<DynamState> state_;
};

inline
void
DynamStatePlainBasis::stateBasisDel() {
  this->state_.clear();
}

inline
void
DynamStatePlainBasis::addState(const DynamState & ds) {
  this->state_.push_back(ds);
}

inline 
void
DynamStatePlainBasis::addStateBasis(const DynamStatePlainBasis * dsb) {
  this->state_.insert(this->state_.end(), dsb->state_.begin(), dsb->state_.end());
}

} // end namespace Pita

#endif /* PITA_DYNAMSTATEPLAINBASIS_H */
