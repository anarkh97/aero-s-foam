#ifndef PITA_DYNAMTIMEINTEGRATOR_H
#define PITA_DYNAMTIMEINTEGRATOR_H

#include "Fwk.h"
#include "Types.h"
#include "DynamState.h"
#include "DynamOps.h"

namespace Pita {

class DynamTimeIntegrator : public Fwk::PtrInterface<DynamTimeIntegrator> {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamTimeIntegrator);

  size_t vectorSize() const { return vectorSize_; }
  
  Seconds timeStepSize() const { return timeStepSize_; }
  const DynamState & initialState() const { return initialState_; }
  const DynamState & currentState() const { return currentState_; }
  Seconds initialTime() const { return initialTime_; }
  Seconds currentTime() const { return currentTime_; }
  TimeStepCount timeStepCount() const { return timeStepCount_; }
  
  virtual void timeStepSizeIs(Seconds dt) = 0;
  virtual void initialConditionIs(const DynamState & initialState, Seconds initialTime = Seconds(0.0)) = 0;
  virtual void currentTimeInc(Seconds currentTime) = 0;
  virtual void timeStepCountInc(TimeStepCount steps = TimeStepCount(1)) = 0;
  
  SliceRank currentSlice() const { return currentSlice_; } // TODO HACK
  void currentSliceIs(SliceRank slice) { setCurrentSlice(slice); performNotification(&NotifieeConst::onCurrentSlice); } // TODO HACK

  class NotifieeConst : public Fwk::BaseMultiNotifiee<const DynamTimeIntegrator, NotifieeConst> {
  public:
    EXPORT_PTRINTERFACE_TYPES(NotifieeConst);
    
    virtual void onInitialCondition() {}
    virtual void onCurrentCondition() {}
    virtual void onCurrentSlice()     {} // TODO HACK

  protected:
    explicit NotifieeConst(const DynamTimeIntegrator * notifier) :
      Fwk::BaseMultiNotifiee<const DynamTimeIntegrator, NotifieeConst>(notifier)
    {}
  };

  virtual void lastNotifieeIs(NotifieeConst * ln) const { notifierDelegate().lastNotifieeIs(ln); }
  virtual void notifieeDel(NotifieeConst * n) const { notifierDelegate().notifieeDel(n); }
  
protected:
  explicit DynamTimeIntegrator(size_t vectorSize);

  typedef GenNotifierDelegate<NotifieeConst>::NotificationType NotificationType;
  void performNotification(NotificationType ln) { notifierDelegate().lastNotificationIs(ln); }

  void setTimeStepSize(Seconds timeStepSize) { timeStepSize_ = timeStepSize; /*if (lastNotifiee()) lastNotifiee()->onTimeStepSize();*/ }
  void setInitialState(const DynamState & initialState) { initialState_ = initialState; /*if (lastNotifiee()) lastNotifiee()->onInitialState();*/ } 
  void setCurrentState(const DynamState & currentState) { currentState_ = currentState; /*if (lastNotifiee()) lastNotifiee()->onCurrentState();*/ } 
  void setInitialTime(Seconds initialTime) { initialTime_ = initialTime; /*if (lastNotifiee()) lastNotifiee()->onInitialTime();*/ }
  void setCurrentTime(Seconds currentTime) { currentTime_ = currentTime; /*if (lastNotifiee()) lastNotifiee()->onCurrentTime();*/ }
  void setTimeStepCount(TimeStepCount timeStepCount) { timeStepCount_ = timeStepCount; /*if (lastNotifiee()) lastNotifiee()->onTimeStepCount();*/ }
  void setCurrentSlice(SliceRank slice) { currentSlice_ = slice; /*if (lastNotifiee()) lastNotifiee()->onCurrentSlice();*/ }

private:
  size_t vectorSize_;
  DynamState initialState_;
  DynamState currentState_;
  Seconds initialTime_; 
  Seconds currentTime_; 
  Seconds timeStepSize_;
  TimeStepCount timeStepCount_;
  SliceRank currentSlice_; // TODO HACK

  GenNotifierDelegate<NotifieeConst> notifierDelegate_;
  GenNotifierDelegate<NotifieeConst> & notifierDelegate() const { return const_cast<DynamTimeIntegrator *>(this)->notifierDelegate_; }

  DISALLOW_COPY_AND_ASSIGN(DynamTimeIntegrator);
};

} // end namespace Pita

#endif /* PITA_DYNAMTIMEINTEGRATOR_H */
