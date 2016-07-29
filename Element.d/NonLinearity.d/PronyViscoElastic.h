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
    PronyViscoElastic(double p1,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3, double p4,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3, double p4, double p5,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3, double p4, double p5, double p6,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3, double p4, double p5, double p6, double p7,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);
    PronyViscoElastic(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8,
                      double ginf, double g1, double tau1, double g2, double tau2, double g3, double tau3);

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

#ifdef _TEMPLATE_FIX_
  #include <Element.d/NonLinearity.d/PronyViscoElastic.C>
#endif

#endif
