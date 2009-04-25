#include <Problems.d/ModalGeomState.h>
#include <Corotational.d/utilities.h>

ModalGeomState::ModalGeomState(int nRBM, int nFlex, int nConstr) :
  q(nFlex, 0.0), lam(nConstr, 0.0), vel(nFlex+nRBM, 0.0),
  acc(nFlex+nRBM, 0.0), numRBM(nRBM), numFlex(nFlex), numConstr(nConstr), qnp1(nFlex, 0.0){
/*PRE: ModalGeomState instantiated with a given number of rigid body modes,
         flexible modes and constraints
 POST: parameterized constructor
*/
  glT[0] = glTnp1[0] = 0.0;
  glT[1] = glTnp1[1] = 0.0;
  glT[2] = glTnp1[2] = 0.0;

  glR[0][0] = 1.0; glR[0][1] = 0.0; glR[0][2] = 0.0;
  glR[1][0] = 0.0; glR[1][1] = 1.0; glR[1][2] = 0.0;
  glR[2][0] = 0.0; glR[2][1] = 0.0; glR[2][2] = 1.0;

  glRnp1[0][0] = 1.0; glRnp1[0][1] = 0.0; glRnp1[0][2] = 0.0;
  glRnp1[1][0] = 0.0; glRnp1[1][1] = 1.0; glRnp1[1][2] = 0.0;
  glRnp1[2][0] = 0.0; glRnp1[2][1] = 0.0; glRnp1[2][2] = 1.0;
}

//------------------------------------------------------------------------------

ModalGeomState::ModalGeomState(const ModalGeomState& mgs) : q(mgs.q),
   lam(mgs.lam), vel(mgs.vel), acc(mgs.acc), numRBM(mgs.numRBM),
   numFlex(mgs.numFlex), numConstr(mgs.numConstr), qnp1(mgs.qnp1){
/*PRE: none
 POST: copy constructor
*/
  for(int i = 0; i < 3; ++i){
    glT[i] = mgs.glT[i];
    glTnp1[i] = mgs.glTnp1[i];
    for(int j = 0; j < 3; ++j){
      glR[i][j] = mgs.glR[i][j];
      glRnp1[i][j] = mgs.glRnp1[i][j];
    }
  }
}

//------------------------------------------------------------------------------

//void ModalGeomState::populate_np1(){
/*PRE: data members glT, glR and q are populated
 POST: populate glTnp1, glRnp1 and qnp1
*/
//}

//------------------------------------------------------------------------------

void ModalGeomState::update(Vector &dsp, double dton2){
/*PRE: *this, vel and acc are approximations to the solution at n+1/2
       dsp is the difference between the next approx and the current one
 POST: update *this to be the next approximation to soln at n+1/2
*/
  glT[0] += dsp[0];  vel[0] += dsp[0]/dton2;  acc[0] += dsp[0]/(dton2*dton2);
  glT[1] += dsp[1];  vel[1] += dsp[1]/dton2;  acc[1] += dsp[1]/(dton2*dton2);
  glT[2] += dsp[2];  vel[2] += dsp[2]/dton2;  acc[2] += dsp[2]/(dton2*dton2);

  glTnp1[0] += 2*dsp[0];
  glTnp1[1] += 2*dsp[1];
  glTnp1[2] += 2*dsp[2];

  double dtheta[3] = { dsp[3], dsp[4], dsp[5] };
  inc_rottensor(dtheta, glR);
  vel[3] += dsp[3]/dton2; acc[3] += dsp[3]/(dton2*dton2);
  vel[4] += dsp[4]/dton2; acc[4] += dsp[4]/(dton2*dton2);
  vel[5] += dsp[5]/dton2; acc[5] += dsp[5]/(dton2*dton2);

  dtheta[0] = 2*dsp[3];
  dtheta[1] = 2*dsp[4];
  dtheta[2] = 2*dsp[5];
  inc_rottensor(dtheta, glRnp1);

  int b, bRBM;
  for(b = 0; b < q.size(); ++b){
    bRBM = b+numRBM;
    q[b] += dsp[bRBM];
    qnp1[b] += 2*dsp[bRBM];
    vel[bRBM] += dsp[bRBM]/dton2;
    acc[bRBM] += dsp[bRBM]/(dton2*dton2);
  }
  for(b = 0; b < numConstr; ++b){
    lam[b] += dsp[numRBM+q.size()+b];
  }

}

//------------------------------------------------------------------------------

void ModalGeomState::midpoint_step_update(double delta,
  ModalGeomState &stepState){
/*PRE: *this, vel, acc and lam are converged to the solution at time step n+1/2
       stepState contains the solution at time step n
 POST: update *this, and stepState to the solution at n+1
*/
#ifdef VICTORY
  printState(" ^^^^^ from ModalGeomState::midpoint_step_update ^^^^^^^^^^^^^^^");
#endif

  int i, j;
  // update velocities
  for(i = 0; i < vel.size(); ++i)
    vel[i] += delta*acc[i];

  // update displacements
  glT[0] = 2.*glT[0] - stepState.glT[0];
  glT[1] = 2.*glT[1] - stepState.glT[1];
  glT[2] = 2.*glT[2] - stepState.glT[2];
  stepState.glT[0] = glT[0];
  stepState.glT[1] = glT[1];
  stepState.glT[2] = glT[2];

  double dtheta[3];
  double dRot[3][3];
  mat_mult_mat(glR, stepState.glR, dRot, 2);
  mat_to_vec(dRot, dtheta);
  inc_rottensor(dtheta, glR);
  for(i = 0; i < 3; ++i)
    for(j = 0; j < 3; ++j)
      stepState.glR[i][j] = glR[i][j];

  for(i = 0; i < q.size(); ++i){
    q[i] = 2.*q[i] - stepState.q[i];
    stepState.q[i] = q[i];
  }
/*
  for(i = 0; i < numConstr; ++i){
    lam[i] = 2*lam[i] - stepState.lam[i];
    stepState.lam[i] = lam[i];
  }
*/

#ifdef VICTORY
  fflush(stderr);
  fprintf(stderr, "glT: %e %e %e\n", glT[0], glT[1], glT[2]);
  fflush(stderr);
  lam.print("lam from mgs::midpoint_step_update", "lambda");
  fflush(stderr);
  fprintf(stderr, "END midpoint_step_update, hence END TIME STEP ######################\n");
#endif

}

//------------------------------------------------------------------------------

void ModalGeomState::printState(const char* text){
/*PRE: none
 POST: print to stderr, the private data members
*/
  if(text){ fprintf(stderr, "%s\n", text); }

  fprintf(stderr, "glT and glR:\n");
  int i;
  for(i = 0; i < 3; ++i){
    fprintf(stderr, " %16.8e", glT[i]);
    fprintf(stderr, "    %16.8e%16.8e%16.8e\n", glR[i][0], glR[i][1], glR[i][2]);
  }

  q.print("", "  q");
  lam.print("", "lam");

  fprintf(stderr, "velocity and acceleration:\n");
  for(i = 0; i < numRBM+numFlex; ++i){
    fprintf(stderr, "  %16.8e  %16.8e\n", vel[i], acc[i]);
  }

}

//------------------------------------------------------------------------------

void ModalGeomState::printRotation(const char* text){

  fprintf(stderr, "%s\n", text);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      fprintf(stderr, " %16.8e", glR[i][j]);
    }
    fprintf(stderr, "\n");
  }

}




