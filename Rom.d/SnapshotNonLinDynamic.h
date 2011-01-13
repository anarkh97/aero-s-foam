#ifndef ROM_SNAPSHOTNONLINDYNAMIC_H
#define ROM_SNAPSHOTNONLINDYNAMIC_H

#include <Problems.d/NonLinDynam.h>

class Domain;
class GeomState;
template <typename Scalar> class GenVector;

#include "BasisOutputFile.h"

#include <memory>

// Specialization of the non-linear dynamics problem enabling the collection the snapshots
class SnapshotNonLinDynamic : public NonLinDynamic {
public:
  explicit SnapshotNonLinDynamic(Domain *);
  virtual ~SnapshotNonLinDynamic();

  virtual void preProcess();

  // Additional operations  
  void saveStateSnapshot(const GeomState &);
  void saveResidualSnapshot(const GenVector<double> &);

private:
  enum SnapshotType { STATE_SNAP, RESIDUAL_SNAP };

  std::auto_ptr<BasisOutputFile> snapFile_[2];
  double (*snapBuffer_)[6];
  int (*dofLocation_)[6];
};

#include <Driver.d/StateUpdater.h>

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
template <typename ProbDescr, typename VecType, typename GeomType>
class SnapshotIncrUpdater : public IncrUpdater<ProbDescr, VecType, GeomType> {
public:
  static void midpointIntegrate(ProbDescr *pbd, VecType &velN,
                                double delta, GeomType *refState, 
                                GeomType *geomState, VecType *dummy1,
                                VecType &dummy2, VecType &dummy3,
                                VecType &dummy4, VecType &acceleration) {
    pbd->saveStateSnapshot(*geomState);

    IncrUpdater<ProbDescr, VecType, GeomType>::midpointIntegrate(
        pbd, velN, delta, refState, geomState,
        dummy1, dummy2, dummy3, dummy4, acceleration);
  }
  
  static double formRHScorrector(ProbDescr *pbd, VecType &inc_displac,
                                 VecType &vel_n, VecType &accel,
                                 VecType &residual, VecType &rhs,
                                 GeomType *geomState) {
    const double result = IncrUpdater<ProbDescr, VecType, GeomType>::formRHScorrector(
        pbd, inc_displac, vel_n, accel, residual, rhs, geomState);

    pbd->saveResidualSnapshot(rhs);

    return result;
  }
};

#endif /* ROM_SNAPSHOTNONLINDYNAMIC_H */
