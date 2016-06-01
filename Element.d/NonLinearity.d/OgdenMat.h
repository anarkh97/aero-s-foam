#ifndef _OGDENMAT_H_
#define _OGDENMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

// This material and those derived from it can now be either isotropic or anisotropic
class OgdenMat : public NLMaterial
{
  protected:
    // isotropic material properties
    double rho; // density
    double mu[9], alpha[9]; // material properties characterizing distortional response
    int N, M; // number of terms in the Ogden series
    double K[9]; // material properties characterizing volumetric response

  public:
    template<int _N, int _M>
    OgdenMat(double _rho, double (&_mu)[_N], double (&_alpha)[_N], double (&_K)[_M]);

    int getNumStates() { return 0; }

    void getStress(Tensor *stress, Tensor &strain, double*, double temp);

    void transformStress(Tensor &stress, Tensor &gradU, Tensor_d0s2_Ss12 &S);

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
