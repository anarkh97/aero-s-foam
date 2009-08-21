#include "ReducedFullTimeSliceImpl.h"

#include <Comm.d/Communicator.h>

#include <memory>

namespace Pita {

namespace ReducedFullTimeSliceProxy {

class LocalHeadProxy;

class LocalTailProxy : public ReducedFullTimeSliceHeadImpl::TailProxy {
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

class RemoteTailProxy : public ReducedFullTimeSliceHeadImpl::TailProxy {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteTailProxy);

  virtual void reducedBasisComponentsIs(const Vector & c, Seed::Status cvgFlag);

  CpuRank targetCpu() const { return targetCpu_; }

  RemoteTailProxy(Communicator * timeComm, CpuRank targetCpu);

private:
  Communicator * timeComm_;
  CpuRank targetCpu_;
};

class LocalHeadProxy : public ReducedFullTimeSliceTailImpl::HeadProxy {
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

class RemoteHeadProxy : public ReducedFullTimeSliceTailImpl::HeadProxy {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteHeadProxy);

  virtual const Vector & reducedBasisComponents(Vector & components) const; // overriden

  CpuRank originCpu() const { return originCpu_; }

  RemoteHeadProxy(Communicator * timeComm, CpuRank originCpu);

private:
  Communicator * timeComm_;
  CpuRank originCpu_;
};

} // end ReducedFullTimeSliceProxy
  
// ReducedFullTimeSliceHeadImpl implementation

ReducedFullTimeSliceHeadImpl::ReducedFullTimeSliceHeadImpl(HalfSliceRank headRank, CpuRank tailCpu,
                                             ReducedFullTimeSliceHeadImpl::Manager * manager,
                                             UpdateProjector * updateProjector,
                                             DynamPropagator * coarsePropagator,
                                             Seed * update,
                                             Seed * nextCorrection) :
  ReducedFullTimeSliceHead(headRank, tailCpu),
  manager_(manager),
  coarsePropagator_(coarsePropagator),
  nextCorrection_(nextCorrection),
  seedReader_(NULL),
  updateProjector_(updateProjector),
  updatedSeedAssembler_(NULL),
  jumpBuilder_(NULL),
  schedulingReactor_(NULL),
  correctionReactor_(NULL),
  update_(update)
{}

void
ReducedFullTimeSliceHeadImpl::tailCpuIs(CpuRank tc) {
  if (tc == this->tailCpu())
    return;
  seedReader_->targetCpuIs(tc);
  correctionReactor_->tailProxyIs(manager_->getTailProxy(tailHalfSlice(), tc));
  setTailCpu(tc);
}

void
ReducedFullTimeSliceHeadImpl::phaseIs(PhaseRank p) {
  // TODO
  setPhase(p);
}

void
ReducedFullTimeSliceHeadImpl::rightPropagatedSeedIs(const Seed * rps) {
  jumpBuilder_->propagatedSeedIs(JumpBuilder::RIGHT, rps);
  schedulingReactor_->notifierIs(rps);
  setRightPropagatedSeed(rps);
}

void
ReducedFullTimeSliceHeadImpl::leftPropagatedSeedIs(const Seed * lps) {
  jumpBuilder_->propagatedSeedIs(JumpBuilder::LEFT, lps);
  updatedSeedAssembler_->propagatedSeedIs(lps);
  setLeftPropagatedSeed(lps);
}

void
ReducedFullTimeSliceHeadImpl::updatedSeedIs(Seed * us) {
  updatedSeedAssembler_->updatedSeedIs(us);
  setUpdatedSeed(us);
}

void
ReducedFullTimeSliceHeadImpl::correctionIs(const ReducedSeed * c) {
  updateProjector_->reducedCorrectionIs(c);
  updatedSeedAssembler_->correctionComponentsIs(c);
  setCorrection(c);
}

const ReducedFullTimeSliceHeadImpl::TailProxy *
ReducedFullTimeSliceHeadImpl::tailProxy() const {
  return correctionReactor_->tailProxy(); 
}

ReducedFullTimeSliceHeadImpl::Manager::Manager(UpdateProjector::Manager * updateProjectorMgr,
                                               DynamPropagator * coarsePropagator,
                                               Seed::Manager * seedMgr,
                                               ReducedSeed::Manager * reducedSeedMgr,
                                               UpdatedSeedAssembler::Manager * updateAssemblerMgr,
                                               Communicator * timeComm,
                                               CpuRank localCpu) :
  updateProjectorMgr_(updateProjectorMgr),
  coarsePropagator_(coarsePropagator),
  seedMgr_(seedMgr),
  reducedSeedMgr_(reducedSeedMgr),
  updateAssemblerMgr_(updateAssemblerMgr),
  seedReaderMgr_(RemoteSeedReader::Manager<HalfSliceRank>::New(timeComm)),
  timeComm_(timeComm),
  localCpu_(localCpu),
  tailMgr_(NULL)
{}

ReducedFullTimeSliceHeadImpl::Manager::~Manager() {
  this->tailManagerIs(NULL);
}

void
ReducedFullTimeSliceHeadImpl::Manager::tailManagerIs(ReducedFullTimeSliceTailImpl::Manager * tm) {
  if (tm == tailMgr_)
    return;

  ReducedFullTimeSliceTailImpl::Manager * previousTailMgr = tailMgr_;
  tailMgr_ = NULL;
  if (previousTailMgr) previousTailMgr->headManagerIs(NULL);

  tailMgr_ = tm;
  if (tm) tm->headManagerIs(this);
}

ReducedFullTimeSliceHeadImpl *
ReducedFullTimeSliceHeadImpl::Manager::createNewInstance(const HalfSliceRank & r) {

  // Policy shared / unique Reductor
  Fwk::String updateProjectorName(String("UP_") + toString(r));
  UpdateProjector::Ptr updateProjector = updateProjectorMgr_->instance(updateProjectorName);
  if (!updateProjector)
    updateProjector = updateProjectorMgr_->instanceNew(updateProjectorName);

  Fwk::String updateName = String("U") + toString(r);
  Seed::Ptr update = seedMgr_->instance(updateName);
  if (!update)
    update = seedMgr_->instanceNew(updateName);

  Fwk::String correctionName = String("C") + toString(r + HalfSliceCount(1));
  Seed::Ptr nextCorrection = seedMgr_->instance(correctionName);
  if (!nextCorrection)
   nextCorrection = seedMgr_->instanceNew(correctionName);

  RemoteSeedReader::Ptr remoteSeedReader = seedReaderMgr_->instanceNew(r);
  remoteSeedReader->notifierIs(nextCorrection.ptr());

  std::auto_ptr<ReducedFullTimeSliceHeadImpl> newSlice(new ReducedFullTimeSliceHeadImpl(r, CpuRank(), this, updateProjector.ptr(), coarsePropagator_.ptr(), update.ptr(), nextCorrection.ptr()));

  newSlice->seedReader_ = remoteSeedReader;

  /* jumpBuilder_ */
  JumpBuilder::Ptr jumpBuilder = JumpBuilder::New();
  
  Fwk::String jumpSeedName = String("J") + toString(r);
  Seed::Ptr jumpSeed = seedMgr_->instance(updateName);
  if (!jumpSeed)
    jumpSeed = seedMgr_->instanceNew(jumpSeedName);

  jumpBuilder->jumpSeedIs(jumpSeed.ptr());
  newSlice->updateProjector_->jumpIs(jumpSeed.ptr()); 
  newSlice->jumpBuilder_ = jumpBuilder;

  /* updatedSeedAssembler */
  Fwk::String seedAssemblerName = String("H_") + toString(r);
  UpdatedSeedAssembler::Ptr updatedSeedAssembler = updateAssemblerMgr_->instance(seedAssemblerName);
  if (!updatedSeedAssembler)
    updatedSeedAssembler = updateAssemblerMgr_->instanceNew(seedAssemblerName);

  newSlice->updatedSeedAssembler_ = updatedSeedAssembler;

  String activityName = String("Projection_") + toString(newSlice->headHalfSlice());
  Activity::Ptr activity = activityManagerInstance()->activityNew(activityName); 

  SchedulingReactor::Ptr schedulingReactor = new SchedulingReactor(newSlice->rightPropagatedSeed(), activity.ptr()); 
  CorrectionReactor::Ptr correctionReactor = new CorrectionReactor(schedulingReactor->activity(), newSlice.get());
  correctionReactor->tailProxyIs(getTailProxy(newSlice->tailHalfSlice(), newSlice->tailCpu())); // TODO HACK
  
  correctionReactor->notifier()->phaseIs(newSlice->phase());
 
  newSlice->schedulingReactor_ = schedulingReactor;
  newSlice->correctionReactor_ = correctionReactor;

  return newSlice.release();
}

ReducedFullTimeSliceHeadImpl::TailProxy *
ReducedFullTimeSliceHeadImpl::Manager::getTailProxy(HalfSliceRank tailRank, CpuRank tailCpu) const {
  if (this->localCpu() == tailCpu) {
    ReducedFullTimeSliceProxy::LocalHeadProxy * headProxy = NULL;
    if (tailMgr_) {
      ReducedFullTimeSliceTailImpl * tail = tailMgr_->instance(tailRank);
      if (tail) {
        headProxy = dynamic_cast<ReducedFullTimeSliceProxy::LocalHeadProxy *>(tail->headProxy());
      }
    }
    return new ReducedFullTimeSliceProxy::LocalTailProxy(headProxy);
  } else {
    return new ReducedFullTimeSliceProxy::RemoteTailProxy(timeComm_, tailCpu);
  }
}

ReducedFullTimeSliceHeadImpl::SchedulingReactor::SchedulingReactor(const Seed * notifier, Activity * activity) :
  Seed::NotifieeConst(notifier),
  activity_(activity)
{}

void
ReducedFullTimeSliceHeadImpl::SchedulingReactor::onState() {
  if (notifier()->status() == Seed::ACTIVE || notifier()->status() == Seed::CONVERGED) {
    activity_->iterationIs(activityManagerInstance()->currentIteration().next());
    activity_->statusIs(Activity::scheduled);
  }
}

ReducedFullTimeSliceHeadImpl::CorrectionReactor::CorrectionReactor(Activity * notifier, ReducedFullTimeSliceHeadImpl * parent) :
  Activity::Notifiee(notifier),
  parent_(parent)
{}

void
ReducedFullTimeSliceHeadImpl::CorrectionReactor::onStatus() {
  if (notifier()->status() == Activity::executing && parent_) {
    size_t reducedBasisSize = parent_->updateProjector_->reducedBasisSize();
    log() << "Reductor basis size = " << reducedBasisSize << "\n";
    if (reducedBasisSize > 0) {
      log() << "Reduced correction on fine time-grid\n";
      const Vector & c = parent_->updateProjector_->updateComponents();
      tailProxy_->reducedBasisComponentsIs(c, Seed::ACTIVE);
      // TODO reconstruct update on full space and commit
      //parent_->update_->statusIs(parent_->updatedSeed()->status());
      //parent_->update_->stateIs(parent_->updatedSeed()->state() - parent_->rightPropagatedSeed()->state());
    } else {
      log() << "Full correction on coarse time-grid\n";
      log() << "Not implemented yet !\n";
      throw "Not implemented";
      //parent_->nextCorrection_->statusIs(Seed::ACTIVE);
      //parent_->coarsePropagator_->initialStateIs(parent_->update_->state());
      //parent_->nextCorrection_->stateIs(parent_->coarsePropagator_->finalState());
    }
  }
}


// ReducedFullTimeSliceTailImpl implemetation

ReducedFullTimeSliceTailImpl::ReducedFullTimeSliceTailImpl(HalfSliceRank tailRank, CpuRank headCpu, Manager * manager, Seed * correction) :
  ReducedFullTimeSliceTail(tailRank, headCpu),
  manager_(manager),
  updatedSeedAssembler_(NULL),
  correction_(correction),
  schedulingReactor_(NULL),
  correctionReactor_(NULL),
  headProxy_(manager->getHeadProxy(headHalfSlice(), headCpu))
{}

void
ReducedFullTimeSliceTailImpl::nextUpdatedSeedIs(Seed * seed) {
  updatedSeedAssembler_->updatedSeedIs(seed);
  setNextUpdatedSeed(seed);
}

void
ReducedFullTimeSliceTailImpl::nextLeftPropagatedSeedIs(const Seed * seed) {
  schedulingReactor_->notifierIs(seed);
  setNextLeftPropagatedSeed(seed);
}

void
ReducedFullTimeSliceTailImpl::headCpuIs(CpuRank hc) {
  headProxy_ = manager_->getHeadProxy(headHalfSlice(), hc);
  seedWriter_->originCpuIs(hc);
  setHeadCpu(hc);
}

void
ReducedFullTimeSliceTailImpl::phaseIs(PhaseRank p) {
  // Reactor/Activity
  schedulingReactor_->activity()->phaseIs(p);
  setPhase(p);
}

ReducedFullTimeSliceTailImpl::Manager::Manager(UpdatedSeedAssembler::Manager * updateAssemblerMgr, Seed::Manager * seedMgr, Communicator * timeComm, CpuRank localCpu) :
  updateAssemblerMgr_(updateAssemblerMgr),
  seedMgr_(seedMgr),
  seedWriterMgr_(RemoteSeedWriter::Manager<HalfSliceRank>::New(timeComm)),
  timeComm_(timeComm),
  localCpu_(localCpu),
  headMgr_(NULL)
{}

ReducedFullTimeSliceTailImpl::Manager::~Manager() {
  this->headManagerIs(NULL);
}

void
ReducedFullTimeSliceTailImpl::Manager::headManagerIs(ReducedFullTimeSliceHeadImpl::Manager * hm) {
  if (hm == headMgr_)
    return;

  ReducedFullTimeSliceHeadImpl::Manager * previousHeadMgr = headMgr_;
  headMgr_ = NULL;
  if (previousHeadMgr) previousHeadMgr->tailManagerIs(NULL);

  headMgr_ = hm;
  if (hm) hm->tailManagerIs(this);
}

ReducedFullTimeSliceTailImpl *
ReducedFullTimeSliceTailImpl::Manager::createNewInstance(const HalfSliceRank & r) {
 
  //String reconstructorName("FTST");
  ///DynamStateReconstructor::Ptr reconstructor = reconstructorMgr_->instance(reconstructorName);
  //if (!reconstructor)
  //  reconstructor = reconstructorMgr_->instanceNew(reconstructorName);
  
  Fwk::String correctionName = String("C") + toString(r);
  Seed::Ptr correction = seedMgr_->instance(correctionName);
  if (!correction)
   correction = seedMgr_->instanceNew(correctionName);
 
  RemoteSeedWriter::Ptr remoteWriter = seedWriterMgr_->instanceNew(r);
  remoteWriter->targetIs(correction.ptr());

  std::auto_ptr<ReducedFullTimeSliceTailImpl> newSlice(new ReducedFullTimeSliceTailImpl(r, CpuRank(), this, correction.ptr()));
  
  /* updatedSeedAssembler */
  Fwk::String seedAssemblerName = String("T_") + toString(r);
  UpdatedSeedAssembler::Ptr updatedSeedAssembler = updateAssemblerMgr_->instance(seedAssemblerName);
  if (!updatedSeedAssembler)
    updatedSeedAssembler = updateAssemblerMgr_->instanceNew(seedAssemblerName);

  newSlice->updatedSeedAssembler_ = updatedSeedAssembler;

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

ReducedFullTimeSliceTailImpl::HeadProxy *
ReducedFullTimeSliceTailImpl::Manager::getHeadProxy(HalfSliceRank headRank, CpuRank headCpu) const {
  if (this->localCpu() == headCpu) {
    ReducedFullTimeSliceProxy::LocalTailProxy * tailProxy = NULL;
    if (headMgr_) {
      ReducedFullTimeSliceHeadImpl * head = headMgr_->instance(headRank);
      if (head) {
        tailProxy = dynamic_cast<ReducedFullTimeSliceProxy::LocalTailProxy *>(head->tailProxy());
      }
    }
    return new ReducedFullTimeSliceProxy::LocalHeadProxy(tailProxy); 
  } else {
    return new ReducedFullTimeSliceProxy::RemoteHeadProxy(timeComm_, headCpu);
  }
}


ReducedFullTimeSliceTailImpl::SchedulingReactor::SchedulingReactor(const Seed * notifier, Activity * activity) :
  Seed::NotifieeConst(notifier),
  activity_(activity)
{}

void
ReducedFullTimeSliceTailImpl::SchedulingReactor::onState() {
  if (notifier()->status() != Seed::SPECIAL) {
    activity_->iterationIs(activity_->currentIteration().next());
    activity_->statusIs(Activity::scheduled);
  }
}

ReducedFullTimeSliceTailImpl::CorrectionReactor::CorrectionReactor(Activity * notifier, ReducedFullTimeSliceTailImpl * parent) :
  Activity::Notifiee(notifier),
  parent_(parent)
{}

void
ReducedFullTimeSliceTailImpl::CorrectionReactor::onStatus() {
  if (notifier()->status() == Activity::executing) {
    if (parent_->nextLeftPropagatedSeed()->status() == Seed::ACTIVE) {
      //log() << "Reconstructor basis size = " << parent_->reconstructor_->reducedBasisSize() << "\n";
      if (false/*parent_->reconstructor_->reducedBasisSize() > 0*/) {
        //Vector components(parent_->reconstructor_->reducedBasisSize());
        //parent_->headProxy_->reducedBasisComponents(components);
        //parent_->reconstructor_->reducedBasisComponentsIs(components); 
        //Seed::Status cvgFlag = headProxy_->convergenceFlag(); // TODO Use cvgFlag
        parent_->correction_->statusIs(Seed::ACTIVE);
        //parent_->correction_->stateIs(parent_->reconstructor_->finalState());
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

namespace ReducedFullTimeSliceProxy {

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
  
} // end namespace ReducedFullTimeSliceProxy

} // end namespace Pita
