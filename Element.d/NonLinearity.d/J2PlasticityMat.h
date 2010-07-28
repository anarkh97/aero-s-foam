#ifndef _J2PLASTICITYMAT_H_
#define _J2PLASTICITYMAT_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>

class J2PlasticityMat : public NLMaterial
{
  public:
    double ematpro[5]; // Young's modulus, Poisson's ratio, mass density, etc.

    J2PlasticityMat(double rho, double E, double nu, double Ep, double sigE)
      { ematpro[0] = E; ematpro[1] = nu; ematpro[2] = rho; ematpro[3] = sigE; ematpro[4] = Ep; }

    int getNumStates() { return 0; }

    void getStress(Tensor *stress, Tensor &strain, double*)
      { std::cerr << "HypoelasticMat::getStress is not implemented\n"; }

    void getTangentMaterial(Tensor *tm, Tensor &strain, double*)
      { std::cerr << "HypoelasticMat::getTangentMaterial is not implemented\n"; }

    void getElasticity(Tensor *tm)
      { std::cerr << "HypoelasticMat::getElasticity is not implemented\n"; }

    void updateStates(Tensor en, Tensor enp,double *state) {}

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*)
      { std::cerr << "HypoelasticMat::getStressAndTangentMaterial is not implemented\n"; }
     
    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double)
      { std::cerr << "HypoelasticMat::integrate is not implemented\n"; }

    void initStates(double *) {}

    double getDensity() { return ematpro[2]; } 
};

#endif
