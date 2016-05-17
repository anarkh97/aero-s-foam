#include <Utils.d/NodeSpaceArray.h>
#include <cmath>
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

template<typename Material>
int
PronyViscoElastic<Material>::getNumStates()
{
  // storage order: h1 (9), h2 (9), h3 (9), sigma(tn-1) (9),
  return 36; // store 4 Prony stresses (36 = 9*4) from previous time step
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
  // Note: this function is called for post-processing.
  // In this case we prefer to output the PK2 stress to be consistent with other
  // materials, and also because the second invariant of the deviatoric PK1 stress
  // can be negative.

  // compute hyperelastic response
  MaterialWrapper<Material>::getStress(_stress, _strain, state, temp);

  // add visco-hyperelastic contribution to long-term hyperelastic response
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress); // symmetric 2nd Piola-Kirchhoff stress tensor
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain); // deformation gradient
  Tensor_d0s2  histry;                                       // container for viscoelastic contribution 

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++) {
       // use stress history stored in state container
       histry[3*i+j] = state[3*i+j] + state[3*i+j+9] + state[3*i+j+18];
    }
  }

#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > F(&strain[0]), P(&histry[0]), s(&(*stress)[0]);
  s *= ginf;
  s += F.inverse()*P; // symmetric 2nd Piola-Kirchhoff stress tensor, S = F^{-1}*P
#else
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = ginf*(*stress)[3*i+j] + histry[3*i+j];
#endif

}

template<typename Material>
void 
PronyViscoElastic<Material>::integrate(Tensor *_stress, Tensor *_tm, Tensor &en, Tensor &_enp,
                                       double *staten, double *statenp, double temp, double dt)
{

  // _stress - constainer for stress tensor that is updated by calling integrate
  // _tm     - first elasticity tensor (time derivative of P-K stress tensor) 
  // en      - doesn't seem to be used 
  // _enp    - deformation gradient (strain)
  // staten  - internal variables before computing material response (get them out of here)
  // statenp - internal variables after computing material response  (put them into here)
  // temp    - material temperature
  // dt      - current time step size

  // compute hyperelastic response
  MaterialWrapper<Material>::integrate(_stress, _tm, en, _enp, staten, statenp, temp, dt);

  // add viscoelastic contribution to long-term hyperelastic response (both stress and tm)
  // using the the time step dt and history variables at t_n (statenp),
  // and also compute the updated history variables (statenp)
  Tensor_d0s4 *tm      = static_cast<Tensor_d0s4 *>(_tm);         // first elasticity tensor
  Tensor_d0s2 *stress  = static_cast<Tensor_d0s2 *>(_stress);     // first P-K stress tensor
  Tensor_d0s2 &enp     = static_cast<Tensor_d0s2 &>(_enp);        // deformation gradient

  double expTau1 = std::exp(-dt/tau1);
  double expTau2 = std::exp(-dt/tau2);
  double expTau3 = std::exp(-dt/tau3);

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      // update stress history variables
      statenp[3*i+j]    = expTau1*staten[3*i+j]    + (g1*(1-expTau1)/(dt/tau1))*((*stress)[3*i+j] - staten[3*i+j+27]);
      statenp[3*i+j+9]  = expTau2*staten[3*i+j+9]  + (g2*(1-expTau2)/(dt/tau2))*((*stress)[3*i+j] - staten[3*i+j+27]);
      statenp[3*i+j+18] = expTau3*staten[3*i+j+18] + (g3*(1-expTau3)/(dt/tau3))*((*stress)[3*i+j] - staten[3*i+j+27]);
      statenp[3*i+j+27] = (*stress)[3*i+j];
      // add viscoelastic contribution to stress
      (*stress)[3*i+j]  *= ginf; 
      (*stress)[3*i+j]  += statenp[3*i+j] + statenp[3*i+j+9] + statenp[3*i+j+18];
      for(int k=0; k<3; k++){
        for(int l=0; l<3; l++){
          // add viscoelastic contribution to tm 
          (*tm)[27*i+9*j+3*k+l] *= ginf+g1*(1-expTau1)/(dt/tau1)+g2*(1-expTau2)/(dt/tau2)+g3*(1-expTau3)/(dt/tau3);
        }
      }
    }
  }

}

template<typename Material>
void
PronyViscoElastic<Material>::integrate(Tensor *_stress, Tensor &en, Tensor &_enp,
                                       double *staten, double *statenp, double temp, double dt)
{

  // compute hyperelastic response
  MaterialWrapper<Material>::integrate(_stress, en, _enp, staten, statenp, temp, dt);

  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress); // first P-K stress tensor
  Tensor_d0s2 &enp    = static_cast<Tensor_d0s2 &>(_enp);    // deformation gradient

  double expTau1 = std::exp(-dt/tau1);
  double expTau2 = std::exp(-dt/tau2);
  double expTau3 = std::exp(-dt/tau3);

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
       // update stress history variables
       statenp[3*i+j]    = expTau1*staten[3*i+j]    + (g1*(1-expTau1)/(dt/tau1))*((*stress)[3*i+j] - staten[3*i+j+27]);
       statenp[3*i+j+9]  = expTau2*staten[3*i+j+9]  + (g2*(1-expTau2)/(dt/tau2))*((*stress)[3*i+j] - staten[3*i+j+27]);
       statenp[3*i+j+18] = expTau3*staten[3*i+j+18] + (g3*(1-expTau3)/(dt/tau3))*((*stress)[3*i+j] - staten[3*i+j+27]);
       statenp[3*i+j+27] = (*stress)[3*i+j];
       // add viscoelastic contribution to stress
       (*stress)[3*i+j]  *= ginf;    
       (*stress)[3*i+j]  += statenp[3*i+j] + statenp[3*i+j+9] + statenp[3*i+j+18];
    }
  }
}

