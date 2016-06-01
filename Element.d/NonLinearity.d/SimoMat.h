#ifndef _SIMOMAT_H_
#define _SIMOMAT_H_

// Isotropic elastic material charactarized by an uncoupled free energy function, quadratic in principle logarithmic stretches
// Reference: Simo, J. C. "Algorithms for static and dynamic multiplicative plasticity that preserve the classical return mapping
//            schemes of the infinitesimal theory." Computer Methods in Applied Mechanics and Engineering 99.1 (1992): 61-112.
//            (section 5)

#include <Element.d/NonLinearity.d/NLMaterial.h>

class SimoMat : public NLMaterial
{
  protected:
    // isotropic material properties
    double rho; // density
    double E, nu; // Young's modulus and Poisson's ratio

  public:
    SimoMat(double _rho, double _E, double _nu);

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
