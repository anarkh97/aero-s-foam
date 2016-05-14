#ifndef _PRONYVISCOELASTIC_H_
#define _PRONYVISCOELASTIC_H_

#include <Element.d/NonLinearity.d/MaterialWrapper.h>

class Tensor;

template<typename Material>
class PronyViscoElastic : public MaterialWrapper<Material>
{
  public:
    PronyViscoElastic(double* params) : MaterialWrapper<Material>(params){};

    int getNumStates();

    void initStates(double *);

    void getStress(Tensor *stress, Tensor &strain, double*, double temp);

    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, double dt=0);

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, double dt=0);

  private:
    double ginf; 
    double g1, tau1;
    double g2, tau2;
    double g3, tau3;
    
};



// set Prony amplitudes and times scales for Linear-elastic, Mooney Rivlin, and NeoHookean

template<>
inline
PronyViscoElastic<IsotropicLinearElastic>::PronyViscoElastic(double *params) : MaterialWrapper<IsotropicLinearElastic>(params)
{
 ginf = params[3];
 g1   = params[4];
 tau1 = params[5];
 g2   = params[6];
 tau2 = params[7];
 g3   = params[8];
 tau3 = params[9];
 std::cout << " ginf = " << ginf << " g1 = " << g1 << " tau1 = " << tau1 << " g2 = " << g2 << " tau2 = " << tau2 << " g3 = " << g3 << " tau3 = " << tau3 << std::endl;
}

template<>
inline
PronyViscoElastic<NeoHookean>::PronyViscoElastic(double *params) : MaterialWrapper<NeoHookean>(params) 
{
 ginf = params[4];
 g1   = params[5];
 tau1 = params[6]; 
 g2   = params[7]; 
 tau2 = params[8];
 g3   = params[9];
 tau3 = params[10];
}

template<>
inline
PronyViscoElastic<MooneyRivlin>::PronyViscoElastic(double *params) : MaterialWrapper<MooneyRivlin>(params)
{
 ginf = params[5];
 g1   = params[6];
 tau1 = params[7];
 g2   = params[8];
 tau2 = params[9];
 g3   = params[10];
 tau3 = params[11];
}

#ifdef _TEMPLATE_FIX_
  #include <Element.d/NonLinearity.d/PronyViscoElastic.C>
#endif

#endif
