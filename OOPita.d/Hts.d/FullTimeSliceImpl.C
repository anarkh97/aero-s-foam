#include "FullTimeSliceImpl.h"

#include <Comm.d/Communicator.h>

#include <memory>

namespace Pita { namespace Hts {

/*class ReductorPropagatorHead : public DynamPropagatorHead {
public:
  EXPORT_PTRINTERFACE_TYPES(ReductorPropagatorHead);

  // Overriden
  virtual void initialStateIs(const DynamState & is) {
    reductor_->initialStateIs(is);
    const Vector & c = reductor_->reducedBasisComponents();
    tailProxy_->reducedBasisComponentsIs(c, Seed::ACTIVE); // TODO flag
  }

  ReductorPropagatorHead(DynamStateReductor * r, FullTimeSliceHeadImpl::TailProxy * tp) :
    DynamPropagatorHead(r ? r->vectorSize() : 0),
    reductor_(r),
    tailProxy_(tp)
  {}

private:
  DynamStateReductor::Ptr reductor_;
  FullTimeSliceHeadImpl::TailProxy::Ptr tailProxy_;
};*/

namespace FullTimeSliceProxy {

class LocalHeadProxy;

class LocalTailProxy : public FullTimeSliceHeadImpl::TailProxy {
public:
  EXPORT_PTRINTERFACE_TYPES(LocalTailProxy);

  explicit LocalTailProxy(LocalHeadProxy * headProxy);
  virtual ~LocalTailProxy();

  virtual void reducedBasisComponentsIs(const Vector & c, Seed::Status cvgFlag);
  const Vector & reducedBasisComponents() const { return *components_; }

private:
  LocalHeadProxy * headProxy_;

  Vector defaultComponents_;
  const Vector * components_;

  friend class LocalHeadProxy;
};

class RemoteTailProxy : public FullTimeSliceHeadImpl::TailProxy {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteTailProxy);

  virtual void reducedBasisComponentsIs(const Vector & c, Seed::Status cvgFlag);

  CpuRank targetCpu() const { return targetCpu_; }

  RemoteTailProxy(Communicator * timeComm, CpuRank targetCpu);

private:
  Communicator * timeComm_;
  CpuRank targetCpu_;
};

class LocalHeadProxy : public FullTimeSliceTailImpl::HeadProxy {
public:
  EXPORT_PTRINTERFACE_TYPES(LocalHeadProxy);

  explicit LocalHeadProxy(LocalTailProxy * tailProxy);
  virtual ~LocalHeadProxy();

  virtual const Vector & reducedBasisComponents(Vector & components) const; // overriden  
  void convergenceFlagIs(Seed::Status cvgFlag) { setConvergenceFlag(cvgFlag); }

private:
  LocalTailProxy * tailProxy_;

  friend class LocalTailProxy;
};

class RemoteHeadProxy : public FullTimeSliceTailImpl::HeadProxy {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteHeadProxy);

  virtual const Vector & reducedBasisComponents(Vector & components) const; // overriden

  CpuRank originCpu() const { return originCpu_; }

  RemoteHeadProxy(Communicator * timeComm, CpuRank originCpu);

private:
  Communicator * timeComm_;
  CpuRank originCpu_;
};

} // end FullTimeSliceProxy
  
// FullTimeSliceHeadImpl implementation

FullTimeSliceHeadImpl::FullTimeSliceHeadImpl(HalfSliceRank headRank, CpuRank tailCpu,
                                             FullTimeSliceHeadImpl::Manager * manager,
                                             DynamStateReductor * reductor,
                                             DynamPropagator * coarsePropagator,
                                             Seed * update,
                                             Seed * correction) :
  FullTimeSliceHead(headRank, tailCpu),
  manager_(manager),
  reductor_(reductor),
  coarsePropagator_(coarsePropagator),
  correction_(correction),
  seedReader_(NULL),
  schedulingReactor_(NULL),
  correctionReactor_(NULL),
  update_(update)
{}

void
FullTimeSliceHeadImpl::updatedSeedIs(const Seed * us) {
  setUpdatedSeed(us);
}

void
FullTimeSliceHeadImpl::rightPropagatedSeedIs(const Seed * rps) {
  setRightPropagatedSeed(rps);
  schedulingReactor_->notifierIs(rps);
}

void
FullTimeSliceHeadImpl::phaseIs(PhaseRank p) {
  schedulingReactor_->activity()->phaseIs(p);
  setPhase(p);
}

void
FullTimeSliceHeadImpl::tailCpuIs(CpuRank tc) {
  if (tc == this->tailCpu())
    return;
  seedReader_->targetCpuIs(tc);
  correctionReactor_->tailProxyIs(manager_->getTailProxy(tailHalfSlice(), tc));
  setTailCpu(tc);
}

const FullTimeSliceHeadImpl::TailProxy *
FullTimeSliceHeadImpl::tailProxy() const {
  return correctionReactor_->tailProxy(); 
}

FullTimeSliceHeadImpl::Manager::Manager(DynamStateReductor::Manager * reductorMgr, DynamPropagator * coarsePropagator, Seed::Manager * seedMgr, Communicator * timeComm, CpuRank localCpu) :
  reductorMgr_(reductorMgr),
  coarsePropagator_(coarsePropagator),
  seedMgr_(seedMgr),
  seedReaderMgr_(RemoteSeedReader::Manager<HalfSliceRank>::New(timeComm)),
  timeComm_(timeComm),
  localCpu_(localCpu),
  tailMgr_(NULL)
{}

FullTimeSliceHeadImpl::Manager::~Manager() {
  this->tailManagerIs(NULL);
}

void
FullTimeSliceHeadImpl::Manager::tailManagerIs(FullTimeSliceTailImpl::Manager * tm) {
  if (tm == tailMgr_)
    return;

  FullTimeSliceTailImpl::Manager * previousTailMgr = tailMgr_;
  tailMgr_ = NULL;
  if (previousTailMgr) previousTailMgr->headManagerIs(NULL);

  tailMgr_ = tm;
  if (tm) tm->headManagerIs(this);
}

FullTimeSliceHeadImpl *
FullTimeSliceHeadImpl::Manager::createNewInstance(const HalfSliceRank & r) {

  // TODO Policy shared / unique Reductor
  Fwk::String reductorName("FTSH");
  DynamStateReductor::Ptr reductor = reductorMgr_->instance(reductorName);
  if (!reductor)
    reductor = reductorMgr_->instanceNew(reductorName);

  Fwk::String updateName = String("U") + toString(r);
  Seed::Ptr update = seedMgr_->instance(updateName);
  if (!update)
    update = seedMgr_->instanceNew(updateName);

  Fwk::String correctionName = String("C") + toString(r + HalfSliceCount(1));
  Seed::Ptr correction = seedMgr_->instance(correctionName);
  if (!correction)
   correction = seedMgr_->instanceNew(correctionName);

  RemoteSeedReader::Ptr remoteSeedReader = seedReaderMgr_->instanceNew(r);
  remoteSeedReader->notifierIs(correction.ptr());

  //DynamPropagatorHead::Ptr propagatorHead = new DynamPropagatorHead(reductor.ptr(), NULL); // TODO

  std::auto_ptr<FullTimeSliceHeadImpl> newSlice(new FullTimeSliceHeadImpl(r, CpuRank(), this, reductor.ptr(), coarsePropagator_.ptr(), update.ptr(), correction.ptr()));
 
  newSlice->seedReader_ = remoteSeedReader;
   
  String activityName = String("Projection_") + toString(newSlice->headHalfSlice());
  Activity::Ptr activity = activityManagerInstance()->activityNew(activityName);

  SchedulingReactor::Ptr schedulingReactor = new SchedulingReactor(newSlice->updatedSeed(), activity.ptr()); 
  CorrectionReactor::Ptr correctionReactor = new CorrectionReactor(schedulingReactor->activity(), newSlice.get());
  correctionReactor->tailProxyIs(getTailProxy(newSlice->tailHalfSlice(), newSlice->tailCpu())); // TODO HACK
  
  correctionReactor->notifier()->phaseIs(newSlice->phase());
 
  newSlice->schedulingReactor_ = schedulingReactor;
  newSlice->correctionReactor_ = correctionReactor;

  return newSlice.release();
}

FullTimeSliceHeadImpl::TailProxy *
FullTimeSliceHeadImpl::Manager::getTailProxy(HalfSliceRank tailRank, CpuRank tailCpu) const {
  if (this->localCpu() == tailCpu) {
    FullTimeSliceProxy::LocalHeadProxy * headProxy = NULL;
    if (tailMgr_) {
      FullTimeSliceTailImpl * tail = tailMgr_->instance(tailRank);
      if (tail) {
        headProxy = dynamic_cast<FullTimeSliceProxy::LocalHeadProxy *>(tail->headProxy());
      }
    }
    return new FullTimeSliceProxy::LocalTailProxy(headProxy);
  } else {
    return new FullTimeSliceProxy::RemoteTailProxy(timeComm_, tailCpu);
  }
}

FullTimeSliceHeadImpl::SchedulingReactor::SchedulingReactor(const Seed * notifier, Activity * activity) :
  Seed::NotifieeConst(notifier),
  activity_(activity)
{}

void
FullTimeSliceHeadImpl::SchedulingReactor::onState() {
  if (notifier()->status() == Seed::ACTIVE || notifier()->status() == Seed::CONVERGED) {
    activity_->iterationIs(activityManagerInstance()->currentIteration().next());
    activity_->statusIs(Activity::scheduled);
  }
}

FullTimeSliceHeadImpl::CorrectionReactor::CorrectionReactor(Activity * notifier, FullTimeSliceHeadImpl * parent) :
  Activity::Notifiee(notifier),
  parent_(parent)
{}

void
FullTimeSliceHeadImpl::CorrectionReactor::onStatus() {
  if (notifier()->status() == Activity::executing && parent_) {
    parent_->update_->statusIs(parent_->updatedSeed()->status());
    parent_->update_->stateIs(parent_->updatedSeed()->state() - parent_->rightPropagatedSeed()->state());
    log() << "Reductor basis size = " << parent_->reductor_->reducedBasisSize() << "\n";
    if (parent_->reductor_->reducedBasisSize() > 0) {
      parent_->reductor_->initialStateIs(parent_->update_->state());
      const Vector & c = parent_->reductor_->reducedBasisComponents();
      tailProxy_->reducedBasisComponentsIs(c, Seed::ACTIVE);
    } else {
      log() << "Correction on coarse grid\n";
      parent_->correction_->statusIs(Seed::ACTIVE);
      parent_->coarsePropagator_->initialStateIs(parent_->update_->state());
      parent_->correction_->stateIs(parent_->coarsePropagator_->finalState());
    }
  }
}


// FullTimeSliceTailImpl implemetation

FullTimeSliceTailImpl::FullTimeSliceTailImpl(HalfSliceRank tailRank, CpuRank headCpu, Manager * manager, DynamStateReconstructor * reconstructor, Seed * correction) :
  FullTimeSliceTail(tailRank, headCpu),
  manager_(manager),
  reconstructor_(reconstructor),
  correction_(correction),
  schedulingReactor_(NULL),
  correctionReactor_(NULL),
  headProxy_(manager->getHeadProxy(headHalfSlice(), headCpu))
{}

void
FullTimeSliceTailImpl::nextUpdatedSeedIs(Seed * seed) {
  setNextUpdatedSeed(seed);
}

void
FullTimeSliceTailImpl::nextLeftPropagatedSeedIs(const Seed * seed) {
  schedulingReactor_->notifierIs(seed);
  setNextLeftPropagatedSeed(seed);
}

void
FullTimeSliceTailImpl::headCpuIs(CpuRank hc) {
  headProxy_ = manager_->getHeadProxy(headHalfSlice(), hc);
  seedWriter_->originCpuIs(hc);
  setHeadCpu(hc);
}

void
FullTimeSliceTailImpl::phaseIs(PhaseRank p) {
  // Reactor/Activity
  schedulingReactor_->activity()->phaseIs(p);
  setPhase(p);
}

FullTimeSliceTailImpl::Manager::Manager(DynamStateReconstructor::Manager * reconstructorMgr, Seed::Manager * seedMgr, Communicator * timeComm, CpuRank localCpu) :
  reconstructorMgr_(reconstructorMgr),
  seedMgr_(seedMgr),
  seedWriterMgr_(RemoteSeedWriter::Manager<HalfSliceRank>::New(timeComm)),
  timeComm_(timeComm),
  localCpu_(localCpu),
  headMgr_(NULL)
{}

FullTimeSliceTailImpl::Manager::~Manager() {
  this->headManagerIs(NULL);
}

void
FullTimeSliceTailImpl::Manager::headManagerIs(FullTimeSliceHeadImpl::Manager * hm) {
  if (hm == headMgr_)
    return;

  FullTimeSliceHeadImpl::Manager * previousHeadMgr = headMgr_;
  headMgr_ = NULL;
  if (previousHeadMgr) previousHeadMgr->tailManagerIs(NULL);

  headMgr_ = hm;
  if (hm) hm->tailManagerIs(this);
}

FullTimeSliceTailImpl *
FullTimeSliceTailImpl::Manager::createNewInstance(const HalfSliceRank & r) {
 
  String reconstructorName("FTST");
  DynamStateReconstructor::Ptr reconstructor = reconstructorMgr_->instance(reconstructorName);
  if (!reconstructor)
    reconstructor = reconstructorMgr_->instanceNew(reconstructorName);

  Fwk::String correctionName = String("C") + toString(r);
  Seed::Ptr correction = seedMgr_->instance(correctionName);
  if (!correction)
   correction = seedMgr_->instanceNew(correctionName);
 
  RemoteSeedWriter::Ptr remoteWriter = seedWriterMgr_->instanceNew(r);
  remoteWriter->targetIs(correction.ptr());

  std::auto_ptr<FullTimeSliceTailImpl> newSlice(new FullTimeSliceTailImpl(r, CpuRank(), this, reconstructor.ptr(), correction.ptr()));

  newSlice->seedWriter_ = remoteWriter;

  String activityName = String("Correction_") + toString(newSlice->tailHalfSlice());
  Activity::Ptr activity = activityManagerInstance()->activityNew(activityName);

  SchedulingReactor::Ptr schedulingReactor = new SchedulingReactor(newSlice->nextLeftPropagatedSeed(), activity.ptr());
  CorrectionReactor::Ptr correctionReactor = new CorrectionReactor(activity.ptr(), newSlice.get());

  correctionReactor->notifier()->phaseIs(newSlice->phase());

  newSlice->schedulingReactor_ = schedulingReactor;
  newSlice->correctionReactor_ = correctionReactor;

  return newSlice.release();
}

FullTimeSliceTailImpl::HeadProxy *
FullTimeSliceTailImpl::Manager::getHeadProxy(HalfSliceRank headRank, CpuRank headCpu) const {
  if (this->localCpu() == headCpu) {
    FullTimeSliceProxy::LocalTailProxy * tailProxy = NULL;
    if (headMgr_) {
      FullTimeSliceHeadImpl * head = headMgr_->instance(headRank);
      if (head) {
        tailProxy = dynamic_cast<FullTimeSliceProxy::LocalTailProxy *>(head->tailProxy());
      }
    }
    return new FullTimeSliceProxy::LocalHeadProxy(tailProxy); 
  } else {
    return new FullTimeSliceProxy::RemoteHeadProxy(timeComm_, headCpu);
  }
}


FullTimeSliceTailImpl::SchedulingReactor::SchedulingReactor(const Seed * notifier, Activity * activity) :
  Seed::NotifieeConst(notifier),
  activity_(activity)
{}

void
FullTimeSliceTailImpl::SchedulingReactor::onState() {
  if (notifier()->status() != Seed::SPECIAL) {
    activity_->iterationIs(activity_->currentIteration().next());
    activity_->statusIs(Activity::scheduled);
  }
}

FullTimeSliceTailImpl::CorrectionReactor::CorrectionReactor(Activity * notifier, FullTimeSliceTailImpl * parent) :
  Activity::Notifiee(notifier),
  parent_(parent)
{}

void
FullTimeSliceTailImpl::CorrectionReactor::onStatus() {
  if (notifier()->status() == Activity::executing) {
    if (parent_->nextLeftPropagatedSeed()->status() == Seed::ACTIVE) {
      log() << "Reconstructor basis size = " << parent_->reconstructor_->reducedBasisSize() << "\n";
      if (parent_->reconstructor_->reducedBasisSize() > 0) {
        Vector components(parent_->reconstructor_->reducedBasisSize());
        parent_->headProxy_->reducedBasisComponents(components);
        parent_->reconstructor_->reducedBasisComponentsIs(components); 
        //Seed::Status cvgFlag = headProxy_->convergenceFlag(); // TODO Use cvgFlag
        parent_->correction_->statusIs(Seed::ACTIVE);
        parent_->correction_->stateIs(parent_->reconstructor_->finalState());
      } else {
        log() << "Correction on coarse grid\n";
        parent_->seedWriter_->vectorSizeIs(parent_->nextLeftPropagatedSeed()->state().vectorSize());
        parent_->seedWriter_->statusIs(RemoteSeedWriter::BUSY);
      }
      DynamState nextUpdatedState = parent_->correction_->state()
                                  + parent_->nextLeftPropagatedSeed()->state();

      parent_->nextUpdatedSeed()->statusIs(parent_->correction_->status());
      parent_->nextUpdatedSeed()->stateIs(nextUpdatedState);
    } else if (parent_->nextLeftPropagatedSeed()->status() == Seed::CONVERGED) {
      parent_->nextUpdatedSeed()->statusIs(Seed::CONVERGED);
      parent_->nextUpdatedSeed()->stateIs(parent_->nextLeftPropagatedSeed()->state());
    }
  }
}

namespace FullTimeSliceProxy {

RemoteTailProxy::RemoteTailProxy(Communicator * timeComm, CpuRank targetCpu) :
  timeComm_(timeComm),
  targetCpu_(targetCpu)
{}

void
RemoteTailProxy::reducedBasisComponentsIs(const Vector & c, Seed::Status cvgFlag) {
  // TODO Handle cvgFlag
  // TODO Asynchronous Send
  int messageTag = timeComm_->myID(); // TODO Id
  if (this->targetCpu() != CpuRank(-1)) {
    timeComm_->sendTo(this->targetCpu().value(), messageTag, c.data(), c.size());
  }
  timeComm_->waitForAllReq();
}

LocalTailProxy::LocalTailProxy(LocalHeadProxy * headProxy) :
  headProxy_(headProxy),
  defaultComponents_(),
  components_(&defaultComponents_)
{}

LocalTailProxy::~LocalTailProxy() {
  if (headProxy_)
    headProxy_->tailProxy_ = NULL;
}

void
LocalTailProxy::reducedBasisComponentsIs(const Vector & c, Seed::Status cvgFlag) {
  this->components_ = &c;
  if (headProxy_) {
    headProxy_->convergenceFlagIs(cvgFlag);
  }
}

LocalHeadProxy::LocalHeadProxy(LocalTailProxy * tailProxy) :
  tailProxy_(tailProxy)
{}

LocalHeadProxy::~LocalHeadProxy() {
  if (tailProxy_)
    tailProxy_->headProxy_ = NULL;
}

const Vector &
LocalHeadProxy::reducedBasisComponents(Vector & components) const {
  if (tailProxy_) {
    components = tailProxy_->reducedBasisComponents();
  }
  return components;
}

RemoteHeadProxy::RemoteHeadProxy(Communicator * timeComm, CpuRank originCpu) :
  timeComm_(timeComm),
  originCpu_(originCpu)
{}

const Vector &
RemoteHeadProxy::reducedBasisComponents(Vector & components) const {
  int messageTag = this->originCpu().value();
  timeComm_->recFrom(messageTag, components.data(), components.size());
  const_cast<RemoteHeadProxy *>(this)->setConvergenceFlag(Seed::INACTIVE); // TODO: Handle CVG flag
  return components;
}
  
} // end namespace FullTimeSliceProxy

} /* end namespace Hts */ } // end namespace Pita
