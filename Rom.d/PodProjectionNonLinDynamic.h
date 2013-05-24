#ifndef ROM_PODPROJECTIONNONLINDYNAMIC_H
#define ROM_PODPROJECTIONNONLINDYNAMIC_H

#include <Problems.d/NonLinDynam.h>
#include <Driver.d/StateUpdater.h>
#include <Rom.d/ModalGeomState.h>
#include <Corotational.d/GeomState.h>

#include <memory>

namespace Rom {

template <typename Scalar> class GenPodProjectionSolver;

class PodProjectionNonLinDynamic : public NonLinDynamic {
public:
  explicit PodProjectionNonLinDynamic(Domain *);
  virtual ~PodProjectionNonLinDynamic();

  virtual int solVecInfo() const;
  virtual void preProcess();

  void readRestartFile(Vector &, Vector &, Vector &, Vector &, ModalGeomState &);
  int getInitState(Vector& d_n, Vector& v_n, Vector &a_n, Vector &v_p);
  void updatePrescribedDisplacement(ModalGeomState *);
  void getConstForce(Vector &);
  void getExternalForce(Vector &, Vector &, int, double, ModalGeomState *, Vector &, Vector &, double);
  void getIncDisplacement(ModalGeomState *, Vector &, ModalGeomState *, bool);
  double formRHScorrector(Vector &, Vector &, Vector &, Vector &, Vector &, ModalGeomState *, double);
  void formRHSpredictor(Vector &, Vector &, Vector &, Vector &, ModalGeomState &, double, double);
  void formRHSinitializer(Vector &, Vector &, Vector &, ModalGeomState &, Vector &, ModalGeomState * = NULL);
  ModalGeomState* createGeomState();
  ModalGeomState* copyGeomState(ModalGeomState *);
  virtual void updateStates(ModalGeomState *, ModalGeomState &);
  double getStiffAndForce(ModalGeomState &, Vector &, Vector &, double = -1, ModalGeomState * = NULL);

  void reBuild(ModalGeomState &, int, double, double);
  void dynamCommToFluid(ModalGeomState *, ModalGeomState *, Vector &, Vector &, Vector &, Vector &, int, int, int);
  void dynamOutput(ModalGeomState *, Vector &, Vector &, double, int, Vector &, Vector &, Vector &, ModalGeomState *) const;
  void getConstraintMultipliers(ModalGeomState &);
  void initializeParameters(ModalGeomState *);
  void updateParameters(ModalGeomState *);

  // Hiding NonLinDynamic::getSolve
  GenPodProjectionSolver<double> *getSolver();
  const GenPodProjectionSolver<double> *getSolver() const;

  // Helper class to be used as template parameter in NLDynamSolver 
  class Updater;

  // Hooks in NLDynamSolver
  virtual double getResidualNorm(const Vector &, ModalGeomState &, double);
  int checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time); // relies on function hiding

protected:
  class Impl;
  GeomState *geomState_Big, *refState_Big;

private:
  virtual bool factorWhenBuilding() const; // Overriden

  void saveMidTime(double); 
  void saveDelta(double);
  void saveStateSnapshot(const GeomState &);
  void saveVelocitySnapshot(const Vector &);
  void saveAccelerationSnapshot(const Vector &);
  void handleResidualSnapshot(const Vector &);

  std::auto_ptr<Impl> impl_;
  std::auto_ptr<Impl> sttImpl_;
  std::auto_ptr<Impl> velImpl_;
  std::auto_ptr<Impl> accImpl_;
  std::auto_ptr<Impl> resImpl_;
  std::auto_ptr<Impl> jacImpl_;
  
  friend class Updater;
  friend class Impl;

  // Disallow copy and assignment
  PodProjectionNonLinDynamic(const PodProjectionNonLinDynamic &);
  PodProjectionNonLinDynamic &operator=(const PodProjectionNonLinDynamic &);
};

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
class PodProjectionNonLinDynamic::Updater : public IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, ModalGeomState> {
public:

  static double integrate(PodProjectionNonLinDynamic *pbd, ModalGeomState *refState, ModalGeomState *geomState,
                          GenVector<double> *du, GenVector<double> &residual,
                          GenVector<double> &elementInternalForce, GenVector<double> &gRes, GenVector<double> &vel_n,
                          GenVector<double> &accel, double midTime) {
    pbd->saveMidTime(midTime);

    return IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, ModalGeomState>::integrate(
        pbd, refState, geomState, du, residual, elementInternalForce, gRes, vel_n, accel, midTime);

/*
    geomState->update(*du, 2);
    return pbd->getStiffAndForce(*geomState, residual, elementInternalForce, midTime, refState);
*/
  }
  
  static void midpointIntegrate(PodProjectionNonLinDynamic *pbd, GenVector<double> &velN,
                                double delta, ModalGeomState *refState, 
                                ModalGeomState *geomState, GenVector<double> *dummy1,
                                GenVector<double> &dummy2, GenVector<double> &dummy3,
                                GenVector<double> &dummy4, GenVector<double> &acceleration, bool zeroRot) {
    pbd->saveDelta(delta);

    IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, ModalGeomState>::midpointIntegrate(
        pbd, velN, delta, refState, geomState,
        dummy1, dummy2, dummy3, dummy4, acceleration, zeroRot);

/* FIXME: needs to be full geomState
    pbd->saveStateSnapshot(*geomState);
*/
    pbd->saveVelocitySnapshot(velN);
    pbd->saveAccelerationSnapshot(acceleration);
  } 

 static double formRHScorrector(PodProjectionNonLinDynamic *pbd, GenVector<double> &inc_displac,
                                GenVector<double> &vel_n, GenVector<double> &accel,
                                GenVector<double> &residual, GenVector<double> &rhs,
                                ModalGeomState *geomState, double delta) {
    const double result = IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, ModalGeomState>::formRHScorrector(
        pbd, inc_displac, vel_n, accel, residual, rhs, geomState, delta);

    pbd->handleResidualSnapshot(rhs);

    return result;
  }
};

} /* end namespace Rom */

#endif /* ROM_PODPROJECTIONNONLINDYNAMIC_H */
