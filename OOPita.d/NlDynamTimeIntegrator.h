#ifndef PITA_NLDYNAMTIMEINTEGRATOR_H
#define PITA_NLDYNAMTIMEINTEGRATOR_H

#include "Fwk.h"
#include "Types.h"
#include "DynamState.h"
#include "DynamTimeIntegrator.h"
#include "NlDynamOps.h"

class GeomState;

namespace Pita {

class PitaNonLinDynamic;

class NlDynamTimeIntegrator : public DynamTimeIntegrator {
public:
  typedef Fwk::Ptr<NlDynamTimeIntegrator> Ptr;
  typedef Fwk::Ptr<const NlDynamTimeIntegrator> PtrConst;

  virtual void timeStepSizeIs(Seconds dt);
  virtual void initialConditionIs(const DynamState & initialState, Seconds initialTime = Seconds(0.0));
  virtual void currentTimeInc(Seconds increment);
  virtual void timeStepCountInc(TimeStepCount steps = TimeStepCount(1));

  static NlDynamTimeIntegrator::Ptr New(PitaNonLinDynamic * pbDesc) {
    return new NlDynamTimeIntegrator(pbDesc); 
  }

  NlDynamOps::Ptr nlDynamOpsNew() const { return NlDynamOps::New(probDesc_); }
  
  const PitaNonLinDynamic * probDesc() const { return probDesc_; }
  
  // Accessors for output
  GeomState * geomState() const { return geomState_; }
  const GenVector<double> & externalForce() const { return externalForce_; }
  
protected:
  explicit NlDynamTimeIntegrator(PitaNonLinDynamic * pbDesc);
  ~NlDynamTimeIntegrator();

private:
  void integrate(unsigned int stepCount);
  void updateDelta(double dt);
  void updateForce();
  
  PitaNonLinDynamic * probDesc_;
  
  GeomState * geomState_;
  GeomState * stepState_;
  GeomState * refState_; 
  GenVector<double> displacement_;
  GenVector<double> velocity_;
  GenVector<double> incDisplac_;
  GenVector<double> gravityForce_;
  
  GenVector<double> elementInternalForce_;
  GenVector<double> residual_;
  GenVector<double> rhs_;
  GenVector<double> aeroForce_;
  GenVector<double> prevIntForce_;
  GenVector<double> externalForce_;
  GenVector<double> stateIncr_;
  GenVector<double> dummyVp_; 
  
  double currTime_;
  double midTime_;
  
  int currStep_;
  
  double localDt_;
  double localDelta_;
  
  int numStages_;
  int maxNumIter_;
  double dlambda_;
};
  
} // end namespace Pita

#endif /* PITA_NLDYNAMTIMEINTEGRATOR_H */
