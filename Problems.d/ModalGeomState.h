#ifndef _MODAL_GEOM_STATE_H
#define _MODAL_GEOM_STATE_H

#include <Math.d/Vector.h>

class ModalGeomState{
/* class used by NLModalDescr for stroing the geometric state of the structure
   contains global translaion; global rotaion; and modal deformations
   passed as a tmeplate parameter GeomType to class NLDynamSolver defined in
     Driver.d/NLDynamProbType.[Ch]
*/
private:

  double glT[3];     // rigid body translation vector
  double glR[3][3];  // rigid body rotation tensor : global = glR*local
  Vector q;          // modal deformation coefficients
  Vector lam;        // Lagrange multipliers to enforce constraints if any

  Vector vel;
  Vector acc;

  // const int numRBM, numFlex, numConstr;
  int numRBM, numFlex, numConstr;  // PJSA: doesn't compile on thunderbird with uninitialized consts

  // the state at time step n+1, used to evaluate constraints
  double glTnp1[3];
  double glRnp1[3][3];
  Vector qnp1;

public:

  ModalGeomState(){}
  ModalGeomState(int nRMB, int nFlex, int nConstr);
  ModalGeomState(const ModalGeomState& mgs);

//  void populate_np1();

  void update(Vector &dsp, double dton2);
  void get_inc_displacement(Vector &incDsp, ModalGeomState &stepState, bool=true)
    { /* NLModalDescr does not use incremental displacement so
           this function does nothing */
    }
  void midpoint_step_update(double delta, ModalGeomState &stepSt);

  void printState(const char* = "");
  void printRotation(const char* = "ModalGeomState.glR");
  void print(){}
  void setVelocity(Vector&, Vector&) {}

  friend class NLModalDescr;
};



#endif
