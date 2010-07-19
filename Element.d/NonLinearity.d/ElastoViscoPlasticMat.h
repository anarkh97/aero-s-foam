#ifndef _ELASTOVISCOPLASTICMAT_H_
#define _ELASTOVISCOPLASTICMAT_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>

class ElastoViscoPlasticMat : public NLMaterial
{
  public:
    double ematpro[8]; // Young's modulus, Poisson's ratio, mass density, etc

    ElastoViscoPlasticMat(double rho, double E, double nu, double p4, double p5, double p6, double p7, double p8)
      { ematpro[0] = E; ematpro[1] = nu; ematpro[2] = rho; ematpro[3] = p4; 
        ematpro[4] = p5; ematpro[5] = p6; ematpro[6] = p7; ematpro[7] = p8; }

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
