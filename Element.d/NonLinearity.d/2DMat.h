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

     GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator();
};

// same equation as ElaLinIsoMat2D but with different Green-Lagrange strain evaluator
// (also known as St. Venant-Kirchhoff hyperelastic material
class StVenantKirchhoffMat2D : public ElaLinIsoMat2D
{
  public:
    StVenantKirchhoffMat2D(StructProp *p) : ElaLinIsoMat2D(p) {}
    StVenantKirchhoffMat2D(double rho, double E, double nu, double t) : ElaLinIsoMat2D(rho, E, nu, t) {}

    GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator();
};


#endif
