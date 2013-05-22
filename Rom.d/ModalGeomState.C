#include <Rom.d/ModalGeomState.h>
#include <Corotational.d/utilities.h>

ModalGeomState::ModalGeomState(int nFlex) :
  q(nFlex, 0.0), vel(nFlex, 0.0),
  acc(nFlex, 0.0), numFlex(nFlex){
/*PRE: ModalGeomState instantiated with a given number of flexible modes
 POST: parameterized constructor
*/
}

//------------------------------------------------------------------------------

ModalGeomState::ModalGeomState(const ModalGeomState& mgs) : q(mgs.q),
   vel(mgs.vel), acc(mgs.acc), numFlex(mgs.numFlex){
/*PRE: none
 POST: copy constructor
*/
}

//------------------------------------------------------------------------------



//------------------------------------------------------------------------------

void ModalGeomState::update(const Vector &dsp, int){
/*PRE: q is approximation to the solution at t^{n+1-alphaf}
       dsp is the difference between the next approx and the current one
 POST: update q to be the next approximation to soln at n+1-alphaf
*/
  q += dsp;
}

//------------------------------------------------------------------------------

void ModalGeomState::midpoint_step_update(Vector &v_n, Vector &a_n, double delta, ModalGeomState &stepState,
                                          double beta, double gamma, double alphaf, double alpham,
                                          bool zeroRot) {
/*PRE: q is converged solution at t^{n+1-alphaf}
       vel and acc are at t^{n}
       stepState contains the solution at t^{n}
 POST: update *this, and stepState to the solution at n+1
*/

  double dt = 2.0*delta;
  double vdcoef, vvcoef, vacoef, avcoef, aacoef;
  vdcoef = (gamma/(dt*beta))/(1-alphaf);
  vvcoef = ((1-(1-alphaf)*gamma/beta)-alphaf)/(1-alphaf);
  vacoef = dt*(2*beta-gamma)/(2*beta);
  avcoef = 1/(dt*gamma);
  aacoef = -(1-gamma)/gamma;

  vel = vdcoef*(q - stepState.q) + vvcoef*v_n + vacoef*a_n;
  acc = avcoef*(vel - v_n) + aacoef*a_n;
  
  double tcoef = 1/(1-alphaf);
  q = tcoef*(q - alphaf*stepState.q);

  stepState.q = q;
  stepState.vel = vel;
  stepState.acc = acc;
  v_n = vel;
  a_n = acc;
}

//------------------------------------------------------------------------------

void ModalGeomState::printState(const char* text){
/*PRE: none
 POST: print to stderr, the private data members
*/
  if(text){ fprintf(stderr, "%s\n", text); }

  q.print("", "  q");

  fprintf(stderr, "velocity and acceleration:\n");
  for(int i = 0; i < numFlex; ++i){
    fprintf(stderr, "  %16.8e  %16.8e\n", vel[i], acc[i]);
  }

}

