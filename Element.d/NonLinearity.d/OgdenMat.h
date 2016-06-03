#ifndef _OGDENMAT_H_
#define _OGDENMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class OgdenMat : public NLMaterial
{
  protected:
    // isotropic material properties
    double rho; // density
    double mu[9], alpha[9]; // material properties characterizing distortional response
    int m, n; // number of terms in the Ogden series
    double K[9]; // material properties characterizing volumetric response

  public:
    template<int _m, int _n>
    OgdenMat(double _rho, double (&_mu)[_m], double (&_alpha)[_m], double (&_K)[_n]);

    int getNumStates() { return 0; }

    void getStress(Tensor *stress, Tensor &strain, double*, double temp);

    void getTangentMaterial(Tensor *tm, Tensor &strain, double*, double temp);

    void getElasticity(Tensor *tm) {};

    void updateStates(Tensor &en, Tensor &enp, double *state, double temp) {};

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double temp);
     
    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, double dt=0);

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, double dt=0);

    void initStates(double *) {};

    double getDensity() { return rho; }

    StrainEvaluator * getStrainEvaluator();

    double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp);

    void print(std::ostream &out) const;

    NLMaterial * clone() const;

    void getMaterialConstants(std::vector<double> &c);
};

#endif
