#ifndef _BILINPLASKINHARDMAT_H_
#define _BILINPLASKINHARDMAT_H_
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>

class StructProp;

//Declaration of the material properties
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

    void getStress(Tensor *stress, Tensor &strain, double *);

    void getTangentMaterial(Tensor *tm, Tensor &strain, double *);

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double *);

    void getElasticity(Tensor *tm);

    void updateStates(Tensor en, Tensor enp, double *state);

    int getNumStates() { return 13; } // the internal variables are : the plastic strain (6 doubles), the center of the yield surface in sigma space (6 doubles) and the equivalent plastic strain (1 double)

    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double);

    void initStates(double *);

    double getDensity() { return rho; }

    StrainEvaluator * getStrainEvaluator();

    double getEquivPlasticStrain(double *statenp) { return statenp[12]; }
};

typedef ElasPlasKinHardMat<0> BilinPlasKinHardMat;
typedef ElasPlasKinHardMat<1> FiniteStrainPlasKinHardMat;
typedef ElasPlasKinHardMat<2> LogStrainPlasKinHardMat;

//template <class T> inline double delta(T a, T b) { return (a==b) ? 1.0 : 0.0; }

#ifdef _TEMPLATE_FIX_
#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.C>
#endif

#endif
