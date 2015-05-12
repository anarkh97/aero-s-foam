#ifndef _MATERIALWRAPPER_H_
#define _MATERIALWRAPPER_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Material.d/Material.h>
#include <Material.d/IsotropicLinearElasticJ2PlasticMaterial.h>
#include <Material.d/IsotropicLinearElasticJ2PlasticPlaneStressMaterial.h>
#include <Material.d/MooneyRivlin.h>

class Tensor;

// This is a wrapper class to interface Adrian Lew's material library with FEM
template<typename Material>
class MaterialWrapper : public NLMaterial
{
  protected:
    Material *mat;
    union {
      double lambda;
      double mu1;
    };
    union {
      double mu;
      double mu2;
    };
    union {
      double kappa;
      int yssrtid;
    };
    double posdefifyTol;

  public:
    MaterialWrapper(Material*);
    MaterialWrapper(double*);
    ~MaterialWrapper() { delete mat; }

    int getNumStates();

    void getStress(Tensor *stress, Tensor &strain, double*, double temp);

    void getTangentMaterial(Tensor *tm, Tensor &strain, double*, double temp);

    void getElasticity(Tensor *tm) {}

    void updateStates(Tensor &en, Tensor &enp, double *state, double temp) {}

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double temp);
     
    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, double dt=0);

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, double dt=0);

    void initStates(double *);

    double getDensity();

    StrainEvaluator * getStrainEvaluator();

    double getEquivPlasticStrain(double *statenp);

    double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp);

    double getPosdefifyTol() { return posdefifyTol; }

    Material* getMaterial() { return mat; }

    void print(std::ostream &out) const;

    void setSDProps(MFTTData *ysst);
    void setSRDProps(MFTTData *yssrt);

    void getMaterialConstants(std::vector<double> &c);
};

template<>
inline
MaterialWrapper<IsotropicLinearElastic>::MaterialWrapper(double *params)
{
  double rho    = params[0];
  double E      = params[1];
  double nu     = params[2];
  lambda = E*nu/((1.+nu)*(1.-2.*nu));
  mu     = E/(2.*(1.+nu));
  mat = new IsotropicLinearElastic(lambda,mu,rho);
  posdefifyTol = -1;
}

template<>
inline
MaterialWrapper<NeoHookean>::MaterialWrapper(double *params) 
{
  double rho    = params[0];
  double E      = params[1]; 
  double nu     = params[2];
  lambda = E*nu/((1.+nu)*(1.-2.*nu));
  mu     = E/(2.*(1.+nu));
  mat = new NeoHookean(lambda,mu,rho);
  posdefifyTol = params[3];
}

template<>
inline
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::MaterialWrapper(double *params)
{
  double rho    = params[0];
  double E      = params[1];
  double nu     = params[2];
  double sigmaY = params[3];
  double K      = params[4];
  double H      = params[5];
  double Tol    = params[6];
  double epsF   = (params[7] <= 0) ? std::numeric_limits<double>::infinity() : params[7];
  yssrtid   = int(params[8]);
  lambda = E*nu/((1.+nu)*(1.-2.*nu));
  mu     = E/(2.*(1.+nu));
  mat = new IsotropicLinearElasticJ2PlasticMaterial(lambda,mu,sigmaY,K,H,Tol,epsF);
  posdefifyTol = -1;
}

template<>
inline
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::MaterialWrapper(double *params)
{
  double rho    = params[0];
  double E      = params[1];
  double nu     = params[2];
  double sigmaY = params[3];
  double K      = params[4];
  double H      = params[5];
  double Tol    = params[6];
  double epsF   = (params[7] <= 0) ? std::numeric_limits<double>::infinity() : params[7];
  yssrtid   = int(params[8]);
  lambda = E*nu/((1.+nu)*(1.-2.*nu));
  mu     = E/(2.*(1.+nu));
  mat = new IsotropicLinearElasticJ2PlasticPlaneStressMaterial(lambda,mu,sigmaY,K,H,Tol,epsF);
  posdefifyTol = -1;
}

template<>
inline
MaterialWrapper<MooneyRivlin>::MaterialWrapper(double *params)
{
  double rho    = params[0];
  mu1    = params[1];
  mu2    = params[2];
  kappa  = params[3];
  mat = new MooneyRivlin(mu1, mu2, kappa, rho);
  posdefifyTol  = params[4];
}

#ifdef _TEMPLATE_FIX_
  #include <Element.d/NonLinearity.d/MaterialWrapper.C>
#endif

#endif
