#ifndef ROM_DISTRSNAPSHOTNONLINDYNAMIC_H
#define ROM_DISTRSNAPSHOTNONLINDYNAMIC_H

#include <Paral.d/MDNLDynam.h>

#include <Corotational.d/DistrGeomState.h>
#include <Feti.d/DistrVector.h>
#include <Driver.d/StateUpdater.h>

#include <memory>

class Domain;

namespace Rom {

// Specialization of the non-linear dynamics problem enabling the collection the snapshots
class DistrSnapshotNonLinDynamic : public MDNLDynamic {
public:
  explicit DistrSnapshotNonLinDynamic(Domain *);

  // Required additional pre- and post-processing
  virtual void preProcess();
  void postProcess() { impl_->postProcess(); }
  
  // Helper class to be used as template parameter in NLDynamSolver 
  class Updater;

protected:
  // Interface to implementation
  class Impl {
  public:
    virtual void lastMidTimeIs(double) = 0;
    virtual void lastDeltaIs(double) = 0;
    virtual void stateSnapshotAdd(const DistrGeomState &) = 0;
    virtual void postProcess() = 0;
    
    virtual ~Impl() {}

  protected:
    Impl() {}

  private:
    // Disallow copy and assigment
    Impl(const Impl &);
    Impl &operator=(const Impl &);
  };

private:
  // Snapshot collection 
  void saveMidTime(double t) { impl_->lastMidTimeIs(t); }
  void saveDelta(double dt) { impl_->lastDeltaIs(dt); }
  void saveStateSnapshot(const DistrGeomState &state) { impl_->stateSnapshotAdd(state); }
 
  std::auto_ptr<Impl> impl_; 

  friend class Updater;
};

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
class DistrSnapshotNonLinDynamic::Updater : public IncrUpdater<DistrSnapshotNonLinDynamic, GenDistrVector<double>, DistrGeomState> {
public:
  static double integrate(DistrSnapshotNonLinDynamic *pbd, DistrGeomState *refState, DistrGeomState *geomState,
                          GenDistrVector<double> *du, GenDistrVector<double> &residual,
                          GenDistrVector<double> &elementInternalForce, GenDistrVector<double> &gRes, GenDistrVector<double> &vel_n,
                          GenDistrVector<double> &accel, double midTime) {
    pbd->saveMidTime(midTime);

    return IncrUpdater<DistrSnapshotNonLinDynamic, GenDistrVector<double>, DistrGeomState>::integrate(
        pbd, refState, geomState, du, residual, elementInternalForce, gRes, vel_n, accel, midTime);
  }

  static void midpointIntegrate(DistrSnapshotNonLinDynamic *pbd, GenDistrVector<double> &velN,
                                double delta, DistrGeomState *refState, 
                                DistrGeomState *geomState, GenDistrVector<double> *dummy1,
                                GenDistrVector<double> &dummy2, GenDistrVector<double> &dummy3,
                                GenDistrVector<double> &dummy4, GenDistrVector<double> &acceleration, bool zeroRot) {
    pbd->saveDelta(delta);

    IncrUpdater<DistrSnapshotNonLinDynamic, GenDistrVector<double>, DistrGeomState>::midpointIntegrate(
        pbd, velN, delta, refState, geomState,
        dummy1, dummy2, dummy3, dummy4, acceleration, zeroRot);
    
    pbd->saveStateSnapshot(*geomState);
  }
};

} /* end namespace Rom */

#endif /* ROM_DISTRSNAPSHOTNONLINDYNAMIC_H */
