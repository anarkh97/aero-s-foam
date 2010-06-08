#ifndef _2DMAT_H_
#define _2DMAT_H_
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Utils.d/NodeSpaceArray.h>

class StructProp;

//Declaration of the material properties
class ElaLinIsoMat2D : public NLMaterial
{
   protected:
     double rho, E, nu, t;

   public:
     ElaLinIsoMat2D(StructProp *p);
     ElaLinIsoMat2D(double _rho, double _E, double _nu, double _t);

     int getNumStates() { return 0; }

     void getStress(Tensor *stress, Tensor &strain, double*);

     void getTangentMaterial(Tensor *tm, Tensor &strain, double*);

     void getElasticity(Tensor *tm) {};

     void updateStates(Tensor en, Tensor enp, double *state) {};

     void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*);
     
     void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                    double *staten, double *statenp, double);

     void initStates(double *) {};
};

#endif
