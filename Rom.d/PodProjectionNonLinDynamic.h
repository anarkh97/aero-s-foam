#ifndef ROM_PODPROJECTIONNONLINDYNAMIC_H
#define ROM_PODPROJECTIONNONLINDYNAMIC_H

#include <Problems.d/NonLinDynam.h>
#include <Driver.d/StateUpdater.h>

#include <memory>

namespace Rom {

template <typename Scalar> class GenPodProjectionSolver;

class PodProjectionNonLinDynamic : public NonLinDynamic {
public:
  explicit PodProjectionNonLinDynamic(Domain *);
  virtual ~PodProjectionNonLinDynamic();

  // Required additional pre-processing
  virtual void preProcess();

  // Hiding NonLinDynamic::getSolve
  GenPodProjectionSolver<double> *getSolver();
  const GenPodProjectionSolver<double> *getSolver() const;

  // Helper class to be used as template parameter in NLDynamSolver 
  class Updater;

  // Hooks in NLDynamSolver
  virtual double getResidualNorm(const Vector &, GeomState &, double);
  int checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time); // relies on function hiding

protected:
  class Impl;

private:
  virtual bool factorWhenBuilding() const; // Overriden

  void saveMidTime(double); 
  void saveDelta(double);
  void saveStateSnapshot(const GeomState &);
  void handleResidualSnapshot(const Vector &);

  std::auto_ptr<Impl> impl_;
  std::auto_ptr<Impl> sttImpl_;
  std::auto_ptr<Impl> resImpl_;
  std::auto_ptr<Impl> jacImpl_;
  
  friend class Updater;
  friend class Impl;

  // Disallow copy and assignment
  PodProjectionNonLinDynamic(const PodProjectionNonLinDynamic &);
  PodProjectionNonLinDynamic &operator=(const PodProjectionNonLinDynamic &);
};

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
class PodProjectionNonLinDynamic::Updater : public IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, GeomState> {
public:

  static double integrate(PodProjectionNonLinDynamic *pbd, GeomState *refState, GeomState *geomState,
                          GenVector<double> *du, GenVector<double> &residual,
                          GenVector<double> &elementInternalForce, GenVector<double> &gRes, GenVector<double> &vel_n,
                          GenVector<double> &accel, double midTime) {
    pbd->saveMidTime(midTime);

    return IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, GeomState>::integrate(
        pbd, refState, geomState, du, residual, elementInternalForce, gRes, vel_n, accel, midTime);
  }
  
  static void midpointIntegrate(PodProjectionNonLinDynamic *pbd, GenVector<double> &velN,
                                double delta, GeomState *refState, 
                                GeomState *geomState, GenVector<double> *dummy1,
                                GenVector<double> &dummy2, GenVector<double> &dummy3,
                                GenVector<double> &dummy4, GenVector<double> &acceleration, bool zeroRot) {
    pbd->saveDelta(delta);

    IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, GeomState>::midpointIntegrate(
        pbd, velN, delta, refState, geomState,
        dummy1, dummy2, dummy3, dummy4, acceleration, zeroRot);

    pbd->saveStateSnapshot(*geomState);
  } 

 static double formRHScorrector(PodProjectionNonLinDynamic *pbd, GenVector<double> &inc_displac,
                                 GenVector<double> &vel_n, GenVector<double> &accel,
                                 GenVector<double> &residual, GenVector<double> &rhs,
                                 GeomState *geomState, double delta) {
    const double result = IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, GeomState>::formRHScorrector(
        pbd, inc_displac, vel_n, accel, residual, rhs, geomState, delta);

    pbd->handleResidualSnapshot(rhs);

    return result;
  }
};

} /* end namespace Rom */

#endif /* ROM_PODPROJECTIONNONLINDYNAMIC_H */
