#ifndef ROM_FOM_H
#define ROM_FOM_H

#include <Driver.d/StateUpdater.h>

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

#endif /* ROM_FOM_H */
