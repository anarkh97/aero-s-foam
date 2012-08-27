#ifndef _NLMATERIAL_H_
#define _NLMATERIAL_H_
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Utils.d/NodeSpaceArray.h>

//Declaration of the material properties
class StrainEvaluator;
template <typename Tensor> class GenStrainEvaluator;
template <int n> class TwoDTensorTypes;

class NLMaterial
{
   public:
     NLMaterial() {}  

     virtual int getNumStates() = 0;  

     virtual void getTangentMaterial(Tensor *tm, Tensor &strain, double *state) = 0;

     virtual void getElasticity(Tensor *tm) = 0;

     virtual void getStress(Tensor *stress, Tensor &strain, double *state) = 0;

     virtual void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double *state) = 0;

     virtual void updateStates(Tensor en, Tensor enp, double *state) = 0;

     virtual void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                            double *staten, double *statenp, double dt = 0.0) = 0;

     virtual void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                            double *staten, double *statenp, double dt = 0.0) = 0;

     virtual void initStates(double *) = 0;

     virtual double getDensity() { return 0; } // PJSA

     virtual StrainEvaluator * getStrainEvaluator() { return NULL; } // return default strain evaluator 

     virtual GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator() { return NULL; }

     virtual double getEquivPlasticStrain(double *statenp) { return 0; }

     virtual double getThickness() { return 0; }

     virtual double getPosdefifyTol() { return -1; }
};

#endif
