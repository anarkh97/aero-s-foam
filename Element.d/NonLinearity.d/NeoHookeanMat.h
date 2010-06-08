#ifndef _NEOHOOKEANMAT_H_
#define _NEOHOOKEANMAT_H_
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Utils.d/NodeSpaceArray.h>
#include <vector>

class StructProp;

// This class is based on Material.d/NeoHookean.cpp (from Adrian Lew)
// modified to take the green-lagrange strain as input and
// output the PK2 stress tensor and the material elasticity tensor
class NeoHookeanMat : public NLMaterial
{
  protected:
    double rho, lambda, mu;

  public:
    NeoHookeanMat(StructProp *p);
    NeoHookeanMat(double _rho, double _E, double _nu);

    int getNumStates() { return 0; }

    void getStress(Tensor *stress, Tensor &strain, double*);

    void getTangentMaterial(Tensor *tm, Tensor &strain, double*);

    void getElasticity(Tensor *tm) {};

    void updateStates(Tensor en,Tensor enp,double *state) {};

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*);
     
    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double);

    void initStates(double *){};

    double getDensity() { return rho; } // PJSA

  private:
    bool GetConstitutiveResponse(const std::vector<double> * strain,
                                       std::vector<double> * stress,
                                       std::vector<double> * tangents) const;
};

#endif
 
