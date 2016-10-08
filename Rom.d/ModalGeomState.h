#ifndef _MODAL_GEOM_STATE_H
#define _MODAL_GEOM_STATE_H

#include <Math.d/Vector.h>

namespace Rom {
class PodProjectionNonLinDynamic;
class LumpedPodProjectionNonLinDynamic;
}

class ModalGeomState {
/* class used by PodProjectionNonLinDynamic for storing the geometric state of the structure
   passed as a tmeplate parameter GeomType to class NLDynamSolver defined in
   Driver.d/NLDynamProbType.[Ch]
*/
protected:

  Vector q;          // modal deformation coefficients
  Vector vel;
  Vector acc;

  int numFlex;

public:

  ModalGeomState(int numModes);
  ModalGeomState(const ModalGeomState& mgs);

  void update(const Vector &, int SO3param = 0);
  void get_inc_displacement(Vector &incDsp, ModalGeomState &stepState, bool)
    {
      incDsp = q - stepState.q;
    }
  void midpoint_step_update(Vector &veloc_n, Vector &accel_n, double delta, ModalGeomState &ss,
                            double beta, double gamma, double alphaf, double alpham,
                            bool zeroRot);

  void printState(const char* = "");
  void setVelocity(const Vector &_vel, int SO3param = 0) { vel = _vel; }
  void setAcceleration(const Vector &_acc, int SO3param = 0) { acc = _acc; } 
  void setVelocityAndAcceleration(Vector &_vel, Vector &_acc) { vel = _vel; acc = _acc; }
  void push_forward(Vector &) {}
  void pull_back(Vector &) {}
  void print() {}

  Vector *getModalq() {return &q; }

  friend class Rom::PodProjectionNonLinDynamic;
  friend class Rom::LumpedPodProjectionNonLinDynamic;
};

#endif
