#ifndef PITA_REDUCEDFULLTIMESLICEIMPL_H
#define PITA_REDUCEDFULLTIMESLICEIMPL_H

#include "ReducedFullTimeSlice.h"

#include "Seed.h"
#include "Activity.h"

#include "DynamPropagator.h" 

#include "JumpBuilder.h"
#include "UpdateProjector.h"
#include "UpdatedSeedAssembler.h"

#include "RemoteSeedReader.h"
#include "RemoteSeedWriter.h"

class Communicator;

namespace Pita {

class ReducedFullTimeSliceHeadImpl : public ReducedFullTimeSliceHead {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedFullTimeSliceHeadImpl);
  class Manager;

  // Overriden mutators
  virtual void tailCpuIs(CpuRank tc);
  virtual void phaseIs(PhaseRank p);
  virtual void rightPropagatedSeedIs(const Seed * rps);
  virtual void leftPropagatedSeedIs(const Seed * lps);
  virtual void updatedSeedIs(Seed * us);
  virtual void correctionIs(const ReducedSeed * c);

  class TailProxy : public Fwk::PtrInterface<TailProxy> {
  public:
    EXPORT_PTRINTERFACE_TYPES(TailProxy);

    virtual void reducedBasisComponentsIs(const Vector & c, Seed::Status cvgFlag) = 0;

  protected:
    TailProxy() {}

  private:
    DISALLOW_COPY_AND_ASSIGN(TailProxy);
  };

  const TailProxy * tailProxy() const;
  TailProxy * tailProxy() { return const_cast<TailProxy *>(const_cast<const ReducedFullTimeSliceHeadImpl *>(this)->tailProxy()); }

protected:
  ReducedFullTimeSliceHeadImpl(HalfSliceRank headRank, CpuRank tailCpu, Manager * manager,
      UpdateProjector * updateProjector, DynamPropagator * coarsePropagator,
      Seed * update, Seed * nextCorrection);
  
  class SchedulingReactor;
  class CorrectionReactor;

  friend class Manager;
  friend class CorrectionReactor;

private:
  Manager * manager_;

  Seed::Ptr update_;
  Seed::Ptr nextCorrection_;
  RemoteSeedReader::Ptr seedReader_;

  DynamPropagator::Ptr coarsePropagator_;

  Fwk::Ptr<JumpBuilder> jumpBuilder_;
  Fwk::Ptr<UpdateProjector> updateProjector_;
  Fwk::Ptr<UpdatedSeedAssembler> updatedSeedAssembler_;

  Fwk::Ptr<SchedulingReactor> schedulingReactor_;
  Fwk::Ptr<CorrectionReactor> correctionReactor_;
};


class ReducedFullTimeSliceTailImpl : public ReducedFullTimeSliceTail {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedFullTimeSliceTailImpl);
  class Manager;

  // Overriden mutators
  virtual void nextUpdatedSeedIs(Seed * seed);
  virtual void nextLeftPropagatedSeedIs(const Seed * seed);
  virtual void headCpuIs(CpuRank hc);
  virtual void phaseIs(PhaseRank p);

  //const DynamStateReconstructor * reconstructor() const { return reconstructor_.ptr(); }

  class HeadProxy : public Fwk::PtrInterface<HeadProxy> {
  public:
    EXPORT_PTRINTERFACE_TYPES(HeadProxy);

    HeadProxy() : cvgFlag_(Seed::INACTIVE) {}

    virtual const Vector & reducedBasisComponents(Vector & components) const = 0; // components is modified
    Seed::Status convergenceFlag() const { return cvgFlag_; }

  protected:
    void setConvergenceFlag(Seed::Status cf) { cvgFlag_ = cf; }

  private:
    Seed::Status cvgFlag_;

    DISALLOW_COPY_AND_ASSIGN(HeadProxy);
  };

  const HeadProxy * headProxy() const { return headProxy_.ptr(); } 
  HeadProxy * headProxy() { return headProxy_.ptr(); } 

protected:
  class SchedulingReactor;
  friend class SchedulingReactor;
  class CorrectionReactor;
  friend class CorrectionReactor;

  ReducedFullTimeSliceTailImpl(HalfSliceRank tailRank, CpuRank headCpu, Manager * manager, Seed * correction);
  
private:
  Manager * manager_;
  
  Fwk::Ptr<UpdatedSeedAssembler> updatedSeedAssembler_;
  Seed::Ptr correction_;
  RemoteSeedWriter::Ptr seedWriter_;

  Fwk::Ptr<SchedulingReactor> schedulingReactor_; 
  Fwk::Ptr<CorrectionReactor> correctionReactor_;

  Fwk::Ptr<HeadProxy> headProxy_;
};

class ReducedFullTimeSliceHeadImpl::Manager : public ReducedFullTimeSliceHead::Manager, private Fwk::GenManagerImpl<ReducedFullTimeSliceHeadImpl, HalfSliceRank> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  ReducedFullTimeSliceHeadImpl * instance(const HalfSliceRank & r) const { return Fwk::GenManagerImpl<ReducedFullTimeSliceHeadImpl, HalfSliceRank>::instance(r); }
  InstanceCount instanceCount() const { return Fwk::GenManagerImpl<ReducedFullTimeSliceHeadImpl, HalfSliceRank>::instanceCount(); }

  ReducedFullTimeSliceHeadImpl * instanceNew(const HalfSliceRank & r) { return Fwk::GenManagerImpl<ReducedFullTimeSliceHeadImpl, HalfSliceRank>::instanceNew(r); } 
  void instanceDel(const HalfSliceRank & r) { Fwk::GenManagerImpl<ReducedFullTimeSliceHeadImpl, HalfSliceRank>::instanceDel(r); }

  CpuRank localCpu() const { return localCpu_; }

  ReducedFullTimeSliceTailImpl::Manager * tailManager() const { return tailMgr_; }
  void tailManagerIs(ReducedFullTimeSliceTailImpl::Manager * tm);

  static Ptr New(UpdateProjector::Manager * updateProjectorMgr, DynamPropagator * coarsePropagator,
                 Seed::Manager * seedMgr, ReducedSeed::Manager * reducedSeedMgr,
                 UpdatedSeedAssembler::Manager * updateAssemblerMgr, 
                 Communicator * timeComm, CpuRank localCpu) {
    return new Manager(updateProjectorMgr, coarsePropagator, seedMgr, reducedSeedMgr, updateAssemblerMgr, timeComm, localCpu);
  }
  
protected:
  Manager(UpdateProjector::Manager * upm, DynamPropagator * cp, Seed::Manager * sm, ReducedSeed::Manager * rsm, UpdatedSeedAssembler::Manager * uam, Communicator * tc, CpuRank lc);
  virtual ~Manager();

  virtual ReducedFullTimeSliceHeadImpl * createNewInstance(const HalfSliceRank & r);
  
  TailProxy * getTailProxy(HalfSliceRank tailRank, CpuRank tailCpu) const;
  
  friend class ReducedFullTimeSliceHeadImpl;

private:
  UpdateProjector::Manager::Ptr updateProjectorMgr_;
 
  DynamPropagator::Ptr coarsePropagator_;

  Seed::Manager::Ptr seedMgr_;
  ReducedSeed::Manager::Ptr reducedSeedMgr_;
  UpdatedSeedAssembler::Manager::Ptr updateAssemblerMgr_;
  RemoteSeedReader::Manager<HalfSliceRank>::Ptr seedReaderMgr_;

  Communicator * timeComm_; 
  CpuRank localCpu_;
  ReducedFullTimeSliceTailImpl::Manager * tailMgr_;
};

class ReducedFullTimeSliceTailImpl::Manager : public ReducedFullTimeSliceTail::Manager, private Fwk::GenManagerImpl<ReducedFullTimeSliceTailImpl, HalfSliceRank> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  ReducedFullTimeSliceTailImpl * instance(const HalfSliceRank & r) const { return Fwk::GenManagerImpl<ReducedFullTimeSliceTailImpl, HalfSliceRank>::instance(r); }
  InstanceCount instanceCount() const { return Fwk::GenManagerImpl<ReducedFullTimeSliceTailImpl, HalfSliceRank>::instanceCount(); }

  ReducedFullTimeSliceTailImpl * instanceNew(const HalfSliceRank & r) { return Fwk::GenManagerImpl<ReducedFullTimeSliceTailImpl, HalfSliceRank>::instanceNew(r); } 
  void instanceDel(const HalfSliceRank & r) { Fwk::GenManagerImpl<ReducedFullTimeSliceTailImpl, HalfSliceRank>::instanceDel(r); }

  CpuRank localCpu() const { return localCpu_; }

  const ReducedFullTimeSliceHeadImpl::Manager * headManager() const { return headMgr_; }
  void headManagerIs(ReducedFullTimeSliceHeadImpl::Manager * hm);

  static Ptr New(UpdatedSeedAssembler::Manager * updateAssemblerMgr,
                 Seed::Manager * seedMgr,
                 Communicator * timeComm,
                 CpuRank localCpu) {
    return new Manager(updateAssemblerMgr, seedMgr, timeComm, localCpu);
  }

protected:
  Manager(UpdatedSeedAssembler::Manager * updateAssemblerMgr, Seed::Manager * seedMgr, Communicator * timeComm, CpuRank localCpu);
  virtual ~Manager();

  friend class ReducedFullTimeSliceTailImpl;

  virtual ReducedFullTimeSliceTailImpl * createNewInstance(const HalfSliceRank & r);
  HeadProxy * getHeadProxy(HalfSliceRank headRank, CpuRank headCpu) const;

private:
  UpdatedSeedAssembler::Manager::Ptr updateAssemblerMgr_;
  Seed::Manager::Ptr seedMgr_;
  RemoteSeedWriter::Manager<HalfSliceRank>::Ptr seedWriterMgr_;

  Communicator * timeComm_;
  CpuRank localCpu_;
  ReducedFullTimeSliceHeadImpl::Manager * headMgr_;
};


class ReducedFullTimeSliceHeadImpl::SchedulingReactor : public Seed::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(SchedulingReactor);

  virtual void onState(); // overriden

  Activity * activity() const { return activity_; }

  SchedulingReactor(const Seed * notifier, Activity * activity);

private:
  Activity * activity_;
};

class ReducedFullTimeSliceHeadImpl::CorrectionReactor : public Activity::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(CorrectionReactor);

  virtual void onStatus(); // overriden

  ReducedFullTimeSliceHead * parent() const { return parent_; }

  const TailProxy * tailProxy() const { return tailProxy_.ptr(); }
  TailProxy * tailProxy() { return tailProxy_.ptr(); }

  void tailProxyIs(TailProxy * tp) { tailProxy_ = tp; }

  CorrectionReactor(Activity * notifier, ReducedFullTimeSliceHeadImpl * parent);

private:
  ReducedFullTimeSliceHeadImpl * parent_;
  Fwk::Ptr<TailProxy> tailProxy_;
};



class ReducedFullTimeSliceTailImpl::SchedulingReactor : public Seed::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(SchedulingReactor);

  virtual void onState(); // overriden

  Activity * activity() const { return activity_; }

  SchedulingReactor(const Seed * notifier, Activity * activity);

private:
  Activity * activity_;
};

class ReducedFullTimeSliceTailImpl::CorrectionReactor : public Activity::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(CorrectionReactor);

  virtual void onStatus();

  CorrectionReactor(Activity * notifier, ReducedFullTimeSliceTailImpl * parent);

private:
  ReducedFullTimeSliceTailImpl * parent_;
};



} // end namespace Pita

#endif /* PITA_REDUCEDFULLTIMESLICEIMPL_H */
