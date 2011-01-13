#ifndef ROM_SNAPSHOTNONLINDYNAMIC_H
#define ROM_SNAPSHOTNONLINDYNAMIC_H

#include <Problems.d/NonLinDynam.h>

class Domain;
#include <Corotational.d/GeomState.h>
#include <Math.d/Vector.h>

#include <memory>

// Specialization of the non-linear dynamics problem enabling the collection the snapshots
class SnapshotNonLinDynamic : public NonLinDynamic {
public:
  enum BasisType { RAW, ORTHOGONAL };

  explicit SnapshotNonLinDynamic(Domain * d);

  // Problem configuration 
  BasisType outputBasisType() const { return outputBasisType_; }

  // Required additional pre- and post-processing
  virtual void preProcess();
  void postProcess() { impl_->postProcess(); }
  
  // Helper class to be used as template parameter in NLDynamSolver 
  class Updater;

  // Interface to implementation
  class Impl {
  public:
    virtual void stateSnapshotAdd(const GeomState &) = 0;
    virtual void residualSnapshotAdd(const GenVector<double> &) = 0;
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
  void saveStateSnapshot(const GeomState & state) { impl_->stateSnapshotAdd(state); }
  void saveResidualSnapshot(const GenVector<double> & residual) { impl_->residualSnapshotAdd(residual); }
 
  BasisType outputBasisType_;
  std::auto_ptr<Impl> impl_; 

  friend class Updater;
};

#include <Driver.d/StateUpdater.h>

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
class SnapshotNonLinDynamic::Updater : public IncrUpdater<SnapshotNonLinDynamic, GenVector<double>, GeomState> {
public:
  static void midpointIntegrate(SnapshotNonLinDynamic *pbd, GenVector<double> &velN,
                                double delta, GeomState *refState, 
                                GeomState *geomState, GenVector<double> *dummy1,
                                GenVector<double> &dummy2, GenVector<double> &dummy3,
                                GenVector<double> &dummy4, GenVector<double> &acceleration) {
    pbd->saveStateSnapshot(*geomState);

    IncrUpdater<SnapshotNonLinDynamic, GenVector<double>, GeomState>::midpointIntegrate(
        pbd, velN, delta, refState, geomState,
        dummy1, dummy2, dummy3, dummy4, acceleration);
  }
  
  static double formRHScorrector(SnapshotNonLinDynamic *pbd, GenVector<double> &inc_displac,
                                 GenVector<double> &vel_n, GenVector<double> &accel,
                                 GenVector<double> &residual, GenVector<double> &rhs,
                                 GeomState *geomState) {
    const double result = IncrUpdater<SnapshotNonLinDynamic, GenVector<double>, GeomState>::formRHScorrector(
        pbd, inc_displac, vel_n, accel, residual, rhs, geomState);

    pbd->saveResidualSnapshot(rhs);

    return result;
  }
};

#endif /* ROM_SNAPSHOTNONLINDYNAMIC_H */
