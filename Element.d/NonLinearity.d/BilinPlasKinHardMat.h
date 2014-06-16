#ifndef _BILINPLASKINHARDMAT_H_
#define _BILINPLASKINHARDMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class StructProp;

template<int e>
class ElasPlasKinHardMat : public NLMaterial
{
  protected:
    // Ep is the tangent modulus from the uniaxial strain-stress curve (Et for DYNA3D Material type 3)
    // theta is hardening parameter which specifies and arbitrary combination of isotropic and kinematic hardening (beta for DYNA3D Material type 3)
    // theta = 0 (default) is purely kinematic hardening, while theta = 1 is purely isotropic hardening
    double rho, E, nu, Ep, sigE, theta;

  public:
    ElasPlasKinHardMat(StructProp *p);
    ElasPlasKinHardMat(double _rho, double _E, double _nu, double _Ep, double _sigE, double _theta = 0)
       {rho = _rho; E = _E; nu = _nu; Ep = _Ep; sigE = _sigE; theta = _theta; }

    void getStress(Tensor *stress, Tensor &strain, double *, double temp);

    void getTangentMaterial(Tensor *tm, Tensor &strain, double *);

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double *, double temp);

    void getElasticity(Tensor *tm);

    void updateStates(Tensor &en, Tensor &enp, double *state, double temp);

    int getNumStates() { return 13; } // the internal variables are : the plastic strain (6 doubles),
                                      // the center of the yield surface in sigma space (6 doubles),
                                      // and the equivalent plastic strain (1 double)

    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp);

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp);

    void initStates(double *);

    double getDensity() { return rho; }

    StrainEvaluator * getStrainEvaluator();

    double getEquivPlasticStrain(double *statenp) { return statenp[12]; }

    bool getBackStress(double *statenp, Tensor *backstress);

    bool getPlasticStrain(double *statenp, Tensor *plasticstrain);

    double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp = 0.0);

    double getDissipatedEnergy(double *statenp);

    void print(std::ostream &out) const;
};

typedef ElasPlasKinHardMat<0> BilinPlasKinHardMat;
typedef ElasPlasKinHardMat<1> FiniteStrainPlasKinHardMat;
typedef ElasPlasKinHardMat<2> LogStrainPlasKinHardMat;

#ifdef _TEMPLATE_FIX_
#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.C>
#endif

#endif
