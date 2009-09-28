#ifndef PITA_HTS_FULLTIMESLICEIMPL_H
#define PITA_HTS_FULLTIMESLICEIMPL_H

#include "FullTimeSliceHeadTail.h"

#include "../Seed.h"
#include "../Activity.h"

#include "../DynamPropagator.h" 

#include "../DynamStateReductor.h"
#include "../DynamStateReconstructor.h"

#include "../RemoteSeedReader.h"
#include "../RemoteSeedWriter.h"

class Communicator;

namespace Pita { namespace Hts {

class FullTimeSliceHeadImpl : public FullTimeSliceHead {
public:
  EXPORT_PTRINTERFACE_TYPES(FullTimeSliceHeadImpl);
  class Manager;

  // Overriden mutators
  virtual void phaseIs(PhaseRank p);
  virtual void tailCpuIs(CpuRank tc);
  virtual void updatedSeedIs(const Seed * us);
  virtual void rightPropagatedSeedIs(const Seed * rps);

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
  TailProxy * tailProxy() { return const_cast<TailProxy *>(const_cast<const FullTimeSliceHeadImpl *>(this)->tailProxy()); }

protected:
  FullTimeSliceHeadImpl(HalfSliceRank headRank, CpuRank tailCpu, Manager * manager,
      DynamStateReductor * reductor, DynamPropagator * coarsePropagator,
      Seed * update, Seed * correction);
  
  class SchedulingReactor;
  class CorrectionReactor;

  friend class Manager;
  friend class CorrectionReactor;

private:
  Manager * manager_;

  DynamStateReductor::Ptr reductor_;
  Seed::Ptr update_;
  Seed::Ptr correction_;
  RemoteSeedReader::Ptr seedReader_;

  DynamPropagator::Ptr coarsePropagator_;

  Fwk::Ptr<SchedulingReactor> schedulingReactor_;
  Fwk::Ptr<CorrectionReactor> correctionReactor_;
};


class FullTimeSliceTailImpl : public FullTimeSliceTail {
public:
  EXPORT_PTRINTERFACE_TYPES(FullTimeSliceTailImpl);
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

  FullTimeSliceTailImpl(HalfSliceRank tailRank, CpuRank headCpu, Manager * manager, DynamStateReconstructor * reconstructor, Seed * correction);
  
private:
  Manager * manager_;
  
  DynamStateReconstructor::Ptr reconstructor_;
  Seed::Ptr correction_;
  RemoteSeedWriter::Ptr seedWriter_;

  Fwk::Ptr<SchedulingReactor> schedulingReactor_; 
  Fwk::Ptr<CorrectionReactor> correctionReactor_;

  Fwk::Ptr<HeadProxy> headProxy_;
};

class FullTimeSliceHeadImpl::Manager : public FullTimeSliceHead::Manager, private Fwk::GenManagerImpl<FullTimeSliceHeadImpl, HalfSliceRank> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  FullTimeSliceHeadImpl * instance(const HalfSliceRank & r) const { return Fwk::GenManagerImpl<FullTimeSliceHeadImpl, HalfSliceRank>::instance(r); }
  InstanceCount instanceCount() const { return Fwk::GenManagerImpl<FullTimeSliceHeadImpl, HalfSliceRank>::instanceCount(); }

  FullTimeSliceHeadImpl * instanceNew(const HalfSliceRank & r) { return Fwk::GenManagerImpl<FullTimeSliceHeadImpl, HalfSliceRank>::instanceNew(r); } 
  void instanceDel(const HalfSliceRank & r) { Fwk::GenManagerImpl<FullTimeSliceHeadImpl, HalfSliceRank>::instanceDel(r); }

  CpuRank localCpu() const { return localCpu_; }

  FullTimeSliceTailImpl::Manager * tailManager() const { return tailMgr_; }
  void tailManagerIs(FullTimeSliceTailImpl::Manager * tm);

  static Ptr New(DynamStateReductor::Manager * reductorMgr, DynamPropagator * coarsePropagator, Seed::Manager * seedMgr, Communicator * timeComm, CpuRank localCpu) {
    return new Manager(reductorMgr, coarsePropagator, seedMgr, timeComm, localCpu);
  }
  
protected:
  Manager(DynamStateReductor::Manager * reductorMgr, DynamPropagator * coarsePropagator, Seed::Manager * seedMgr, Communicator * timeComm, CpuRank localCpu);
  virtual ~Manager();

  virtual FullTimeSliceHeadImpl * createNewInstance(const HalfSliceRank & r);
  
  TailProxy * getTailProxy(HalfSliceRank tailRank, CpuRank tailCpu) const;
  
  friend class FullTimeSliceHeadImpl;

private:
  DynamStateReductor::Manager::Ptr reductorMgr_;
 
  DynamPropagator::Ptr coarsePropagator_;

  Seed::Manager::Ptr seedMgr_;
  RemoteSeedReader::Manager<HalfSliceRank>::Ptr seedReaderMgr_;

  Communicator * timeComm_; 
  CpuRank localCpu_;
  FullTimeSliceTailImpl::Manager * tailMgr_;
};

class FullTimeSliceTailImpl::Manager : public FullTimeSliceTail::Manager, private Fwk::GenManagerImpl<FullTimeSliceTailImpl, HalfSliceRank> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  FullTimeSliceTailImpl * instance(const HalfSliceRank & r) const { return Fwk::GenManagerImpl<FullTimeSliceTailImpl, HalfSliceRank>::instance(r); }
  InstanceCount instanceCount() const { return Fwk::GenManagerImpl<FullTimeSliceTailImpl, HalfSliceRank>::instanceCount(); }

  FullTimeSliceTailImpl * instanceNew(const HalfSliceRank & r) { return Fwk::GenManagerImpl<FullTimeSliceTailImpl, HalfSliceRank>::instanceNew(r); } 
  void instanceDel(const HalfSliceRank & r) { Fwk::GenManagerImpl<FullTimeSliceTailImpl, HalfSliceRank>::instanceDel(r); }

  CpuRank localCpu() const { return localCpu_; }

  const FullTimeSliceHeadImpl::Manager * headManager() const { return headMgr_; }
  void headManagerIs(FullTimeSliceHeadImpl::Manager * hm);

  static Ptr New(DynamStateReconstructor::Manager * reconstructorMgr,
      Seed::Manager * seedMgr,
      Communicator * timeComm,
      CpuRank localCpu) {
    return new Manager(reconstructorMgr, seedMgr, timeComm, localCpu);
  }

protected:
  Manager(DynamStateReconstructor::Manager * reconstructorMgr, Seed::Manager * seedMgr, Communicator * timeComm, CpuRank localCpu);
  virtual ~Manager();

  friend class FullTimeSliceTailImpl;

  virtual FullTimeSliceTailImpl * createNewInstance(const HalfSliceRank & r);
  HeadProxy * getHeadProxy(HalfSliceRank headRank, CpuRank headCpu) const;

private:
  DynamStateReconstructor::Manager::Ptr reconstructorMgr_;
  Seed::Manager::Ptr seedMgr_;
  RemoteSeedWriter::Manager<HalfSliceRank>::Ptr seedWriterMgr_;

  Communicator * timeComm_;
  CpuRank localCpu_;
  FullTimeSliceHeadImpl::Manager * headMgr_;
};

class FullTimeSliceHeadImpl::SchedulingReactor : public Seed::NotifieeConst {
public:
  typedef Fwk::Ptr<SchedulingReactor> Ptr;
  typedef Fwk::Ptr<const SchedulingReactor> PtrConst;

  virtual void onState(); // overriden

  Activity * activity() const { return activity_; }

  SchedulingReactor(const Seed * notifier, Activity * activity);

private:
  Activity * activity_;
};

class FullTimeSliceHeadImpl::CorrectionReactor : public Activity::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(CorrectionReactor);

  virtual void onStatus(); // overriden

  FullTimeSliceHead * parent() const { return parent_; }

  const TailProxy * tailProxy() const { return tailProxy_.ptr(); }
  TailProxy * tailProxy() { return tailProxy_.ptr(); }

  void tailProxyIs(TailProxy * tp) { tailProxy_ = tp; }

  CorrectionReactor(Activity * notifier, FullTimeSliceHeadImpl * parent);

private:
  FullTimeSliceHeadImpl * parent_;
  Fwk::Ptr<TailProxy> tailProxy_;
};



class FullTimeSliceTailImpl::SchedulingReactor : public Seed::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(SchedulingReactor);

  virtual void onState(); // overriden

  Activity * activity() const { return activity_; }

  SchedulingReactor(const Seed * notifier, Activity * activity);

private:
  Activity * activity_;
};

class FullTimeSliceTailImpl::CorrectionReactor : public Activity::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(CorrectionReactor);

  virtual void onStatus();

  CorrectionReactor(Activity * notifier, FullTimeSliceTailImpl * parent);

private:
  FullTimeSliceTailImpl * parent_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_FULLTIMESLICEIMPL_H */
