#include <Utils.d/NodeSpaceArray.h>
#include <cmath>

template<typename Material>
int
PronyViscoElastic<Material>::getNumStates()
{
  // storage order: h1 (6), h2 (6), h3 (6), sigma(tn-1) (6),
  return 24; // store 4 Prony PK2 stresses (24 = 6*4) from previous time step
}

template<typename Material>
void
PronyViscoElastic<Material>::initStates(double *state)
{
  for(int i=0; i<getNumStates(); ++i) state[i] = 0;
}

template<typename Material>
void
PronyViscoElastic<Material>::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  // compute hyperelastic response
  Material::getStress(_stress, _strain, state, temp);

  // add visco-hyperelastic contribution to long-term hyperelastic response
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress); // 2nd Piola-Kirchhoff stress tensor

  for(int i = 0; i < 6; i++) {
    // use stress history stored in state container
    (*stress)[i] = ginf*(*stress)[i] + state[i] + state[i+6] + state[i+12];
  }
}

template<typename Material>
void 
PronyViscoElastic<Material>::integrate(Tensor *_stress, Tensor *_tm, Tensor &en, Tensor &enp,
                                       double *staten, double *statenp, double temp, Tensor *cache, double dt)
{
  // _stress - constainer for 2nd P-K stress tensor that is updated by calling integrate
  // _tm     - elasticity tensor (derivative of 2nd P-K stress tensor) 
  // enp     - Green-Lagrange strain tensor
  // staten  - internal variables before computing material response (get them out of here)
  // statenp - internal variables after computing material response  (put them into here)
  // temp    - material temperature
  // dt      - current time step size

  // compute hyperelastic response
  Material::integrate(_stress, _tm, en, enp, staten, statenp, temp, cache, dt);

  // add viscoelastic contribution to long-term hyperelastic response (both stress and tm)
  // using the the time step dt and history variables at t_n (statenp),
  // and also compute the updated history variables (statenp)
  Tensor_d0s4_Ss12s34 *tm  = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);

  double expTau1 = std::exp(-dt/tau1);
  double expTau2 = std::exp(-dt/tau2);
  double expTau3 = std::exp(-dt/tau3);

  for(int i = 0; i < 6; i++) {
    // update stress history variables
    statenp[i]    = expTau1*staten[i]    + (g1*(1-expTau1)/(dt/tau1))*((*stress)[i] - staten[i+18]);
    statenp[i+6]  = expTau2*staten[i+6]  + (g2*(1-expTau2)/(dt/tau2))*((*stress)[i] - staten[i+18]);
    statenp[i+12] = expTau3*staten[i+12] + (g3*(1-expTau3)/(dt/tau3))*((*stress)[i] - staten[i+18]);
    statenp[i+18] = (*stress)[i];
    // add viscoelastic contribution to stress
    (*stress)[i]  *= ginf; 
    (*stress)[i]  += statenp[i] + statenp[i+6] + statenp[i+12];
    for(int j=0; j<6; j++){
      // add viscoelastic contribution to tm 
      (*tm)[i][j] *= ginf+g1*(1-expTau1)/(dt/tau1)+g2*(1-expTau2)/(dt/tau2)+g3*(1-expTau3)/(dt/tau3);
    }
  }
}

template<typename Material>
void
PronyViscoElastic<Material>::integrate(Tensor *_stress, Tensor &en, Tensor &enp,
                                       double *staten, double *statenp, double temp, Tensor *cache, double dt)
{
  // compute hyperelastic response
  Material::integrate(_stress, en, enp, staten, statenp, temp, cache, dt);

  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress); // second P-K stress tensor

  double expTau1 = std::exp(-dt/tau1);
  double expTau2 = std::exp(-dt/tau2);
  double expTau3 = std::exp(-dt/tau3);

  for(int i = 0; i < 6; i++) {
    // update stress history variables   
    statenp[i]    = expTau1*staten[i]    + (g1*(1-expTau1)/(dt/tau1))*((*stress)[i] - staten[i+18]);
    statenp[i+6]  = expTau2*staten[i+6]  + (g2*(1-expTau2)/(dt/tau2))*((*stress)[i] - staten[i+18]);
    statenp[i+12] = expTau3*staten[i+12] + (g3*(1-expTau3)/(dt/tau3))*((*stress)[i] - staten[i+18]);
    statenp[i+18] = (*stress)[i];
    // add viscoelastic contribution to stress
    (*stress)[i]  *= ginf; 
    (*stress)[i]  += statenp[i] + statenp[i+6] + statenp[i+12];
  }
}

