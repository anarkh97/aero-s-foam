#ifndef PITA_TIMESLICENETWORK_H
#define PITA_TIMESLICENETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "LocalTimeSlice.h"
#include "RemoteTimeSlice.h"
#include "TimeSliceMapping.h"
#include "CommManager.h"
#include "SeedInitializer.h"

#include "deque"

namespace Pita {

class TimeSliceNetwork : public Fwk::PtrInterface<TimeSliceNetwork> {
public:
  typedef Fwk::Ptr<TimeSliceNetwork> Ptr;
  typedef Fwk::Ptr<const TimeSliceNetwork> PtrConst;

  CpuRank localCpuRank() const { return myCpu_; }
  TimeSliceMapping::Ptr sliceMapping() const { return mapping_; }
  LocalTimeSlice::Manager::Ptr localSliceManager() const { return localManager_; }
  RemoteTimeSlice::Manager::Ptr remoteSliceManager() const { return remoteManager_; }

  IterationRank lastIteration() const { return lastIteration_; }
  void lastIterationIs(IterationRank i) { lastIteration_ = i; }
  
  static TimeSliceNetwork::Ptr New(CpuRank myCpu, TimeSliceMapping * mapping, LocalTimeSlice::Manager * localManager, RemoteTimeSlice::Manager * remoteManager, SeedInitializer * initializer) {
    return new TimeSliceNetwork(myCpu, mapping, localManager, remoteManager, initializer);
  }
  
protected:
  TimeSliceNetwork(CpuRank myCpu, TimeSliceMapping * mapping, LocalTimeSlice::Manager * localManager, RemoteTimeSlice::Manager * remoteManager, SeedInitializer * initializer);

  class ConvergenceReactor;
  friend class ConvergenceReactor;
  
private:
  typedef std::deque<TimeSlice::Ptr> SliceMap;
  SliceMap slice_;

  CpuRank myCpu_;
  TimeSliceMapping::Ptr mapping_;
  TimeSliceMapping::SliceIteratorConst sliceIt_;
  
  LocalTimeSlice::Manager::Ptr localManager_;
  RemoteTimeSlice::Manager::Ptr remoteManager_;
  SeedInitializer::Ptr initializer_;

  IterationRank lastIteration_;

  Fwk::Ptr<ConvergenceReactor> convergenceReactor_;

  void scheduleActivities();
  void forceConvergence();
};

class TimeSliceNetwork::ConvergenceReactor : public Activity::Notifiee {
public:
  typedef Fwk::Ptr<ConvergenceReactor> Ptr;
  typedef Fwk::Ptr<const ConvergenceReactor> PtrConst;

  void onStatus();

  ConvergenceReactor(Activity * notifier, TimeSliceNetwork * parent);

private:
  TimeSliceNetwork * parent_;
};

}

#endif
