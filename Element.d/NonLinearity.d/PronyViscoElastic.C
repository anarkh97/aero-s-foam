#include <Utils.d/NodeSpaceArray.h>
#include <cmath>

template<typename Material, class tensor_policy>
PronyViscoElastic<Material, tensor_policy>::PronyViscoElastic(double p1,
                      double _ginf, double _g1, double _tau1, double _g2, double _tau2, double _g3, double _tau3)
 : Material(p1), ginf(_ginf), g1(_g1), tau1(_tau1), g2(_g2), tau2(_tau2), g3(_g3), tau3(_tau3) {}

template<typename Material, class tensor_policy>
PronyViscoElastic<Material, tensor_policy>::PronyViscoElastic(double p1, double p2,
                      double _ginf, double _g1, double _tau1, double _g2, double _tau2, double _g3, double _tau3)
 : Material(p1,p2), ginf(_ginf), g1(_g1), tau1(_tau1), g2(_g2), tau2(_tau2), g3(_g3), tau3(_tau3) {}

template<typename Material, class tensor_policy>
PronyViscoElastic<Material, tensor_policy>::PronyViscoElastic(double p1, double p2, double p3,
                      double _ginf, double _g1, double _tau1, double _g2, double _tau2, double _g3, double _tau3)
 : Material(p1,p2,p3), ginf(_ginf), g1(_g1), tau1(_tau1), g2(_g2), tau2(_tau2), g3(_g3), tau3(_tau3) {}

template<typename Material, class tensor_policy>
PronyViscoElastic<Material, tensor_policy>::PronyViscoElastic(double p1, double p2, double p3, double p4,
                      double _ginf, double _g1, double _tau1, double _g2, double _tau2, double _g3, double _tau3)
 : Material(p1,p2,p3,p4), ginf(_ginf), g1(_g1), tau1(_tau1), g2(_g2), tau2(_tau2), g3(_g3), tau3(_tau3) {}

template<typename Material, class tensor_policy>
PronyViscoElastic<Material, tensor_policy>::PronyViscoElastic(double p1, double p2, double p3, double p4, double p5,
                      double _ginf, double _g1, double _tau1, double _g2, double _tau2, double _g3, double _tau3)
 : Material(p1,p2,p3,p4,p5), ginf(_ginf), g1(_g1), tau1(_tau1), g2(_g2), tau2(_tau2), g3(_g3), tau3(_tau3) {}

template<typename Material, class tensor_policy>
PronyViscoElastic<Material, tensor_policy>::PronyViscoElastic(double p1, double p2, double p3, double p4, double p5, double p6,
                      double _ginf, double _g1, double _tau1, double _g2, double _tau2, double _g3, double _tau3)
 : Material(p1,p2,p3,p4,p5,p6), ginf(_ginf), g1(_g1), tau1(_tau1), g2(_g2), tau2(_tau2), g3(_g3), tau3(_tau3) {}

template<typename Material, class tensor_policy>
PronyViscoElastic<Material, tensor_policy>::PronyViscoElastic(double p1, double p2, double p3, double p4, double p5, double p6, double p7,
                      double _ginf, double _g1, double _tau1, double _g2, double _tau2, double _g3, double _tau3)
 : Material(p1,p2,p3,p4,p5,p6,p7), ginf(_ginf), g1(_g1), tau1(_tau1), g2(_g2), tau2(_tau2), g3(_g3), tau3(_tau3) {}

template<typename Material, class tensor_policy>
PronyViscoElastic<Material, tensor_policy>::PronyViscoElastic(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8,
                      double _ginf, double _g1, double _tau1, double _g2, double _tau2, double _g3, double _tau3)
 : Material(p1,p2,p3,p4,p5,p6,p7,p8), ginf(_ginf), g1(_g1), tau1(_tau1), g2(_g2), tau2(_tau2), g3(_g3), tau3(_tau3) {}

template<typename Material, class tensor_policy>
PronyViscoElastic<Material, tensor_policy>::PronyViscoElastic(double p1, int i1, int i2, int i3, double p2, double p3, double p4, FabricMap::StrainMeasure i4,
                      double _ginf, double _g1, double _tau1, double _g2, double _tau2, double _g3, double _tau3)
 : Material(p1, i1, i2, i3, p2, p3, p4, i4), ginf(_ginf), g1(_g1), tau1(_tau1), g2(_g2), tau2(_tau2), g3(_g3), tau3(_tau3) {}

template<typename Material, class tensor_policy>
PronyViscoElastic<Material, tensor_policy>::PronyViscoElastic(double p1, int i1, int i2, double p2, double p3, double p4, double p5, double p6, double p7,
                      FabricMat::StrainMeasure i3, double _ginf, double _g1, double _tau1, double _g2, double _tau2, double _g3, double _tau3)
 : Material(p1, i1, i2, p2, p3, p4, p5, p6, p7, i3), ginf(_ginf), g1(_g1), tau1(_tau1), g2(_g2), tau2(_tau2), g3(_g3), tau3(_tau3) {}

template<typename Material, class tensor_policy>
NLMaterial *
PronyViscoElastic<Material, tensor_policy>::clone() const
{
  return new PronyViscoElastic<Material, tensor_policy>(*this);
}

template<typename Material, class tensor_policy>
int
PronyViscoElastic<Material, tensor_policy>::getNumStates()
{
  std::size_t s = tensor_policy::stride;

  // storage order: h1 (s), h2 (s), h3 (s), sigma(tn-1) (s)
  return s * 4; // store 4 Prony stress tensors from previous time step
}

template<typename Material, class tensor_policy>
void
PronyViscoElastic<Material, tensor_policy>::initStates(double *state)
{
  for(int i = 0; i < getNumStates(); ++i) state[i] = 0;
}

template<typename Material, class tensor_policy>
void
PronyViscoElastic<Material, tensor_policy>::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  // compute hyperelastic response
  Material::getStress(_stress, _strain, state, temp);

  // add visco-hyperelastic contribution to long-term hyperelastic response
  typename tensor_policy::d0s2_S *stress = static_cast<typename tensor_policy::d0s2_S *>(_stress); // stress tensor

  std::size_t s = tensor_policy::stride;

  for(int i = 0; i < s; i++) {
    // use stress history stored in state container
    (*stress)[i] = ginf*(*stress)[i] + state[i] + state[i+s] + state[i+2*s];
  }
}

template<typename Material, class tensor_policy>
void
PronyViscoElastic<Material, tensor_policy>::integrate(Tensor *_stress, Tensor *_tm, Tensor &en, Tensor &enp,
                                                      double *staten, double *statenp, double temp, Tensor *cache, double dt)
{
  // _stress - constainer for stress tensor that is updated by calling integrate
  // _tm     - elasticity tensor (derivative of stress tensor)
  // enp     - strain tensor
  // staten  - internal variables before computing material response (get them out of here)
  // statenp - internal variables after computing material response  (put them into here)
  // temp    - material temperature
  // dt      - current time step size

  // compute hyperelastic response
  Material::integrate(_stress, _tm, en, enp, staten, statenp, temp, cache, dt);

  typename tensor_policy::d0s4_S *tm  = static_cast<typename tensor_policy::d0s4_S *>(_tm);
  typename tensor_policy::d0s2_S *stress = static_cast<typename tensor_policy::d0s2_S *>(_stress);

  std::size_t s = tensor_policy::stride;

  if(statenp == 0) {
    for(int i = 0; i < s; i++) {
      (*stress)[i] *= ginf;
      for(int j = 0; j < s; j++) (*tm)[i][j] *= ginf;
    }
  }
  else {
    // add viscoelastic contribution to long-term hyperelastic response (both stress and tm)
    // using the the time step dt and history variables at t_n (staten),
    // and also compute the updated history variables (statenp)

    double g[3] = { g1, g2, g3 }, tau[3] = { tau1, tau2, tau3 }, expTau[3], gxxTau[3];
    for(int k = 0; k < 3; ++k) {
      if(tau[k] == 0) {
        expTau[k] = 0;
        gxxTau[k] = 0;
      }
      else {
        double x;
        if((x = dt/tau[k]) <= std::numeric_limits<double>::epsilon()) {
          expTau[k] = 1;
          gxxTau[k] = g[k];
        }
        else {
          expTau[k] = std::exp(-x);
          gxxTau[k] = g[k]*(1 - expTau[k])/x;
        }
      }
    }

    for(int i = 0; i < s; i++) {
      // update stress history variables
      statenp[i]     = expTau[0]*staten[i]     + gxxTau[0]*((*stress)[i] - staten[i+3*s]);
      statenp[i+s]   = expTau[1]*staten[i+s]   + gxxTau[1]*((*stress)[i] - staten[i+3*s]);
      statenp[i+2*s] = expTau[2]*staten[i+2*s] + gxxTau[2]*((*stress)[i] - staten[i+3*s]);
      statenp[i+3*s] = (*stress)[i];
      // add viscoelastic contribution to stress
      (*stress)[i] *= ginf;
      (*stress)[i] += statenp[i] + statenp[i+s] + statenp[i+2*s];
      for(int j = 0; j < s; j++) {
        // add viscoelastic contribution to tm
        (*tm)[i][j] *= ginf + gxxTau[0] + gxxTau[1] + gxxTau[2];
      }
    }
  }
}

template<typename Material, class tensor_policy>
void
PronyViscoElastic<Material, tensor_policy>::integrate(Tensor *_stress, Tensor &en, Tensor &enp,
                                                      double *staten, double *statenp, double temp, Tensor *cache, double dt)
{
  // compute hyperelastic response
  Material::integrate(_stress, en, enp, staten, statenp, temp, cache, dt);

  typename tensor_policy::d0s2_S *stress = static_cast<typename tensor_policy::d0s2_S *>(_stress); // stress tensor

  std::size_t s = tensor_policy::stride;

  if(statenp == 0) {
    for(int i = 0; i < s; i++) {
      (*stress)[i] *= ginf;
    }
  }
  else {
    // add viscoelastic contribution to long-term hyperelastic response
    // using the the time step dt and history variables at t_n (staten),
    // and also compute the updated history variables (statenp)

    double g[3] = { g1, g2, g3 }, tau[3] = { tau1, tau2, tau3 }, expTau[3], gxxTau[3];
    for(int k = 0; k < 3; ++k) {
      if(tau[k] == 0) {
        expTau[k] = 0;
        gxxTau[k] = 0;
      }
      else {
        double x;
        if((x = dt/tau[k]) <= std::numeric_limits<double>::epsilon()) {
          expTau[k] = 1;
          gxxTau[k] = g[k];
        }
        else {
          expTau[k] = std::exp(-x);
          gxxTau[k] = g[k]*(1 - expTau[k])/x;
        }
      }
    }

    for(int i = 0; i < s; i++) {
      // update stress history variables
      statenp[i]     = expTau[0]*staten[i]     + gxxTau[0]*((*stress)[i] - staten[i+3*s]);
      statenp[i+s]   = expTau[1]*staten[i+2*s] + gxxTau[1]*((*stress)[i] - staten[i+3*s]);
      statenp[i+2*s] = expTau[2]*staten[i+3*s] + gxxTau[2]*((*stress)[i] - staten[i+3*s]);
      statenp[i+3*s] = (*stress)[i];
      // add viscoelastic contribution to stress
      (*stress)[i] *= ginf;
      (*stress)[i] += statenp[i] + statenp[i+s] + statenp[i+2*s];
    }
  }
}

template<typename Material, class tensor_policy>
void
PronyViscoElastic<Material, tensor_policy>::print(std::ostream &out) const
{
  out << "Visco";
  Material::print(out);
  out << " " << g1 << " " << tau1 << " " << g2 << " " << tau2 << " " << g3 << " " << tau3;
}

