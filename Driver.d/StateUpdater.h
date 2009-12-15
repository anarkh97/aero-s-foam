#ifndef _STATE_UPDATER_H_
#define _STATE_UPDATER_H_

template< class ProbDescr, class VecType, class GeomType >
class IncrUpdater  {

public:

  typedef VecType StateIncr;
  typedef int RefState;

  static void update(RefState *, GeomType *unp, StateIncr *, VecType &du) { unp->update(du); }
  static RefState *initRef(GeomType *)  { return 0; }
  static StateIncr *initInc(GeomType *, VecType *v)  { return new VecType(*v); }
  static void copyState(GeomType *, RefState *)  {};
  static void zeroInc(StateIncr *du)  { du->zero(); }
  static double integrate(ProbDescr *pbd, RefState *, GeomType *geomState,
		  StateIncr *du, VecType &residual, 
		  VecType &elementInternalForce, VecType &gRes) {
    geomState->update(*du);
    return pbd->getStiffAndForce(*geomState, residual, elementInternalForce, gRes);
  }
  static double integrate(ProbDescr *pbd, RefState *, GeomType *geomState,
		  StateIncr *du, VecType &residual, 
		  VecType &elementInternalForce, VecType &gRes, VecType& vel_n,
                  VecType &accel, double midTime) {
    geomState->update(*du);
    return pbd->getStiffAndForce(*geomState, residual, elementInternalForce, midTime);
  }
  static void midpointIntegrate(ProbDescr *pbd, 
                  VecType &velN, double delta, GeomType *refState, 
                  GeomType *geomState,
		  StateIncr *, VecType &, 
		  VecType &, VecType &, VecType &acceleration) {

//   This method should return the proper reactions.
//   the computation requires the accelerations.
//
//
//


     geomState->midpoint_step_update(velN, acceleration, delta, *refState);
  }
  static void updateIncr(StateIncr *du, VecType &ddu) { *du = ddu; }

  static void formRHSpredictor(ProbDescr *pbd, VecType &vel_n, VecType &accel,
    VecType &residual, VecType &rhs, GeomType *geomState, double midtime){
    pbd->formRHSpredictor(vel_n, residual, rhs, *geomState, midtime);
  }
  static double formRHScorrector(ProbDescr *pbd, VecType &inc_displac, VecType &vel_n,
    VecType &accel, VecType &residual, VecType &rhs, GeomType *geomState){
    return pbd->formRHScorrector(inc_displac, vel_n, accel, residual, rhs); 
  }

  static void copyTo(RefState *refState, GeomType *geomState, GeomType *stepState, StateIncr *stateIncr, VecType &v, VecType &a, VecType &vp, VecType &force,
                     RefState *refState_bk, GeomType *geomState_bk, GeomType *stepState_bk, StateIncr *stateIncr_bk, VecType &v_bk, VecType &a_bk, VecType &vp_bk, VecType &force_bk){
    if(refState)
      *refState_bk  = *refState; 
    if(geomState)
      *geomState_bk = *geomState;
    if(stepState)
      *stepState_bk = *stepState;
    if(stateIncr)
      *stateIncr_bk = *stateIncr;
    vp_bk = vp;
    a_bk = a;
    v_bk = v;
    force_bk = force;
  }
    

};

//-------------------------------------------------------------------------------------------------

template< class ProbDescr, class VecType, class GeomType >
class NLModalUpdater : public IncrUpdater< ProbDescr, VecType, GeomType> {

public:

  typedef VecType StateIncr;
  typedef int RefState;

  static double integrate(ProbDescr *pbd, RefState *rs, GeomType *geomState,
    StateIncr *du, VecType &residual, VecType &elementInternalForce,
    VecType &gRes, VecType& vel_n, VecType &accel, double midTime){

//************************************************
//    pbd->derivTest(*geomState, vel_n, accel, residual);
//    exit(0);
//************************************************
  
    geomState->update(*du, pbd->getDelta());
    return pbd->getStiffAndForce(*geomState, residual);
  }

  static void formRHSpredictor(ProbDescr *pbd, VecType &vel_n, VecType &accel,
    VecType &residual, VecType &rhs, GeomType *geomState, double midtime){
    pbd->formRHSpredictor(residual, rhs, *geomState);
  }

  static double formRHScorrector(ProbDescr *pbd, VecType &inc_displac, VecType &vel_n,
    VecType &accel, VecType &residual, VecType &rhs, GeomType *geomState){
    return pbd->formRHScorrector(residual, rhs, *geomState); 
  }

  static void midpointIntegrate(ProbDescr *pbd, VecType &vel_n, double delta,
    GeomType *stepState, GeomType *geomState, StateIncr *, VecType &,
    VecType &, VecType &, VecType &accel) {
    geomState->midpoint_step_update(delta, *stepState);
  }

  static void copyTo(RefState *refState, GeomType *geomState, GeomType *stepState, StateIncr *stateIncr, VecType &v, VecType &a, VecType &vp, VecType &force, 
                     RefState *refState_bk, GeomType *geomState_bk, GeomType *stepState_bk, StateIncr *stateIncr_bk, VecType &v_bk, VecType &a_bk, VecType &vp_bk, VecType &force_bk){
    if(refState)
      *refState_bk  = *refState; 
    if(geomState)
      *geomState_bk = *geomState;
    if(stepState)
      *stepState_bk = *stepState;
    if(stateIncr)
      *stateIncr_bk = *stateIncr;
    a_bk = a;
    v_bk = v;
    vp_bk = vp;
    force_bk = force;
  }

};

//-------------------------------------------------------------------------------------------------

template< class ProbDescr, class VecType, class GeomType >
class TotalUpdater  {

public:
  typedef VecType StateIncr;
  typedef GeomType RefState;

  static void update(RefState *rs, GeomType *unp, StateIncr *incr, VecType &du)  {
    *incr += du; unp->update(*rs, *incr);  }

  static VecType *initInc(GeomType *, VecType *v)  { return new VecType(*v);  }

  static RefState *initRef(GeomType *u)  { return new RefState(*u); }
  static void copyState(GeomType *gn, RefState *gp)  { *gp = *gn; }
  static void zeroInc(StateIncr *incr)  { incr->zero(); }
  static double integrate(ProbDescr *pbd, GeomType *sn, GeomType *snp,
		  StateIncr *du, VecType &residual, 
		  VecType &elementInternalForce, VecType &totalRes, VecType=0, VecType=0,
		  double=0) {
    return pbd->integrate(*sn, *snp, *du, residual, 
                          elementInternalForce, totalRes);
  }
  static void updateIncr(StateIncr *du, VecType &ddu) { *du += ddu; }
  static void midpointIntegrate(ProbDescr *pbd, 
                  VecType &velN, double delta, RefState *sn, 
                  GeomType *snp,
		  StateIncr *du, VecType &residual, 
		  VecType &elementInternalForce, VecType &totalRes, VecType &acceleration) {
     snp->midpoint_step_update(velN, delta, *sn, acceleration);
     
     pbd->integrate(*sn, *snp, *du, residual, 
                          elementInternalForce, totalRes);
     *sn = *snp;
  }

  static void formRHSpredictor(ProbDescr *pbd, VecType &vel_n, VecType &accel,
    VecType &residual, VecType &rhs, GeomType *geomState, double midtime){
    pbd->formRHSpredictor(vel_n, residual, rhs, *geomState, midtime);
  }
  static double formRHScorrector(ProbDescr *pbd, VecType &inc_displac, VecType &vel_n,
    VecType &accel, VecType &residual, VecType &rhs, GeomType *geomState){
    return pbd->formRHScorrector(inc_displac, vel_n, residual, rhs); 
  }

  static void copyTo(RefState *refState, GeomType *geomState, GeomType *stepState, StateIncr *stateIncr, VecType &v, VecType &a, VecType &vp, VecType &force,
                     RefState *refState_bk, GeomType *geomState_bk, GeomType *stepState_bk, StateIncr *stateIncr_bk, VecType &v_bk, VecType &a_bk, VecType &vp_bk, VecType &force_bk){
    if(refState)
      *refState_bk  = *refState; 
    if(geomState)
      *geomState_bk = *geomState;
    if(stepState)
      *stepState_bk = *stepState;
    if(stateIncr)
      *stateIncr_bk = *stateIncr;
    a_bk = a;
    v_bk = v;
    vp_bk = vp;
    force_bk = force;
  }

};

//-------------------------------------------------------------------------------------------------


#endif
