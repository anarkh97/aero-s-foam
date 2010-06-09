#ifndef PITA_ZEROSHAREDSTATE_H
#define PITA_ZEROSHAREDSTATE_H

#include "NamedTask.h"
#include "SharedState.h"

namespace Pita {

template <typename S>
class ZeroSharedState : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(ZeroSharedState);

  SharedState<S> * target() const { return target_.ptr(); }
  void targetIs(SharedState<S> * t) { target_ = t; }

  size_t stateSize() const { return stateSize_; }
  void stateSizeIs(size_t ss) { stateSize_ = ss; }

  // overriden
  virtual void iterationIs(IterationRank i) {
    target_->stateIs(S(stateSize(), 0.0));
    target_->statusIs(Seed::INACTIVE);
    target_->iterationIs(i);
  }

  static Ptr New(const String & name) {
    return new ZeroSharedState(name);
  }

private:
  typename SharedState<S>::Ptr target_;
  size_t stateSize_;

  explicit ZeroSharedState(const String & name) :
    NamedTask(name),
    stateSize_(0)
  {}

  DISALLOW_COPY_AND_ASSIGN(ZeroSharedState);
};

} /* end namespace Pita */

#endif /* PITA_ZEROSHAREDSTATE_H */
