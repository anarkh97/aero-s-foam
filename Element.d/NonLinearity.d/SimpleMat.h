#ifndef _SIMPLEMAT_H_
#define _SIMPLEMAT_H_
#include <Element.d/Element.h>
#include <Utils.d/NodeSpaceArray.h>
#include <Material.d/Material.h>

class StructProp;
class SimpleMaterial;

// This is a wrapper class to interface Adrian Lew's SimpleMaterial class with FEM
class SimpleMat : public NLMaterial
{
  protected:
    SimpleMaterial *mat;

  public:
    SimpleMat(int type, double rho, double E, double nu);

    int getNumStates() { return 0; }

    void getStress(Tensor *stress, Tensor &strain, double*);

    void getTangentMaterial(Tensor *tm, Tensor &strain, double*);

    void getElasticity(Tensor *tm) {}

    void updateStates(Tensor en, Tensor enp, double *state) {}

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*);
     
    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double);

    void initStates(double *) {}

    double getDensity();
};

#endif
 
