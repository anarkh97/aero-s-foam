#ifndef _PRONYVISCOELASTIC_H_
#define _PRONYVISCOELASTIC_H_

#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NeoHookeanMat.h>
#include <Element.d/NonLinearity.d/MooneyRivlinMat.h>

class Tensor;

template<typename Material>
class PronyViscoElastic : public Material
{
  public:
    PronyViscoElastic(double* params);

    int getNumStates();

    void initStates(double *);

    void getStress(Tensor *stress, Tensor &strain, double*, double temp);

    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, Tensor *cache, double dt=0);

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, Tensor *cache, double dt=0);

  private:
    double ginf; 
    double g1, tau1;
    double g2, tau2;
    double g3, tau3;
};


// set Prony amplitudes and times scales for Linear-elastic, Mooney Rivlin, and NeoHookean

template<>
inline
PronyViscoElastic<ElaLinIsoMat>::PronyViscoElastic(double *params) : ElaLinIsoMat(params[0], params[1], params[2], 0, 0)
{
 ginf = params[3];
 g1   = params[4];
 tau1 = params[5];
 g2   = params[6];
 tau2 = params[7];
 g3   = params[8];
 tau3 = params[9];
}

template<>
inline
PronyViscoElastic<NeoHookeanMat>::PronyViscoElastic(double *params) : NeoHookeanMat(params[0], params[1], params[2]) 
{
 ginf = params[3];
 g1   = params[4];
 tau1 = params[5];
 g2   = params[6];
 tau2 = params[7];
 g3   = params[8];
 tau3 = params[9];
}

template<>
inline
PronyViscoElastic<MooneyRivlinMat>::PronyViscoElastic(double *params) : MooneyRivlinMat(params[0], params[1], params[2], params[3])
{
 ginf = params[4];
 g1   = params[5];
 tau1 = params[6];
 g2   = params[7];
 tau2 = params[8];
 g3   = params[9];
 tau3 = params[10];
}

#ifdef _TEMPLATE_FIX_
  #include <Element.d/NonLinearity.d/PronyViscoElastic.C>
#endif

#endif
