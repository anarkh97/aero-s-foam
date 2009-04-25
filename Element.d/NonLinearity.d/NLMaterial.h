#ifndef _NLMATERIAL_H_
#define _NLMATERIAL_H_
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Utils.d/NodeSpaceArray.h>

//Declaration of the material properties
class Material {
   public:
     Material() {}  
     virtual int getNumStates() = 0;  

     virtual void getTangentMaterial(Tensor *tm,Tensor &strain, double *state) = 0;

     virtual void getElasticity(Tensor *tm) = 0;

     virtual void getStress(Tensor *stress,Tensor &strain, double *state)=0;

     virtual void getStressAndTangentMaterial(Tensor *stress,Tensor *tm, Tensor &strain, double *state)= 0;

     virtual void updateStates(Tensor en, Tensor enp, double *state)=0;
     virtual void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                      double *staten, double *statenp, double dt = 0.0) = 0;
     virtual void initStates(double *) = 0;

     
};

#endif
