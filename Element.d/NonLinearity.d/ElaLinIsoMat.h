#ifndef _ELALINISOMAT_H_
#define _ELALINISOMAT_H_
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Utils.d/NodeSpaceArray.h>

class StructProp;

//Declaration of the material properties
class ElaLinIsoMat : public Material {
   
   protected:

     double rho, E, nu;

   public:
     ElaLinIsoMat(StructProp *p);
     ElaLinIsoMat(double _rho, double _E, double _nu) {
      rho = _rho; nu = _nu; E = _E;
     }

     int getNumStates(){return 0;}

     void getStress(Tensor *stress,Tensor &strain,double*);

     void getTangentMaterial(Tensor *tm,Tensor &strain,double*);

     void getElasticity(Tensor *tm){};

     void updateStates(Tensor en,Tensor enp,double *state){};

     void getStressAndTangentMaterial(Tensor *stress,Tensor *tm,Tensor &strain,double*);
     
     void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                      double *staten, double *statenp, double);
     void initStates(double *){};
};

#endif
 
