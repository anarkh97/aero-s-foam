#ifndef _SIMOPLASTICMAT_H_
#define _SIMOPLASTICMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Utils.d/MFTT.h>

// Isotropic elasto-plastic material whose elastic response is charactarized by an uncoupled free energy function, quadratic in
// principle logarithmic stretches, and featuring a multiplicative decomposition of the elastic and plastic deformation gradients.
// Reference: Simo, J. C. "Algorithms for static and dynamic multiplicative plasticity that preserve the classical return mapping
//            schemes of the infinitesimal theory." Computer Methods in Applied Mechanics and Engineering 99.1 (1992): 61-112.

class SimoPlasticMat : public NLMaterial
{
  protected:
    // Ep is the tangent modulus from the uniaxial strain-stress curve (Et for DYNA3D Material type 3)
    // theta is hardening parameter which specifies and arbitrary combination of isotropic and kinematic hardening (beta for DYNA3D Material type 3)
    // theta = 0 (default) is purely kinematic hardening, while theta = 1 is purely isotropic hardening
    double rho, E, nu, Ep, sigE, theta;
    // strain dependent material properties
    MFTTData *ysst;
    // tolerance for convergence of nonlinear solve
    double tol;
    // rate dependent material properties
    int yssrtid;
    MFTTData *yssrt;

  public:
    SimoPlasticMat(double rho, double E, double nu, double _Ep, double _sigE, double _theta = 0,
                   double _tol = 1e-6, int _yssrtid = 0);

    void getStress(Tensor *stress, Tensor &strain, double *, double temp);

    void getTangentMaterial(Tensor *tm, Tensor &strain, double *, double temp);

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double *, double temp);

    void getElasticity(Tensor *tm);

    void updateStates(Tensor &en, Tensor &enp, double *state, double temp);

    int getNumStates();

    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp,
                   Tensor *cache, double dt=0);

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp,
                   Tensor *cache, double dt=0);

    void initStates(double *);

    StrainEvaluator * getStrainEvaluator();

    double getEquivPlasticStrain(double *statenp);

    bool getBackStress(double *statenp, Tensor *backstress);

    bool getPlasticStrain(double *statenp, Tensor *plasticstrain);

    double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp = 0.0);

    double getDissipatedEnergy(double *statenp);

    void print(std::ostream &out) const;

    void setSDProps(MFTTData *_ysst) { if(sigE < 0 && _ysst && _ysst->getID() == -int(sigE)) ysst = _ysst; }

    void setSRDProps(MFTTData *_yssrt) { if(yssrtid > 0 && _yssrt && _yssrt->getID() == yssrtid) yssrt = _yssrt; }
};

#endif
