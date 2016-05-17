#ifndef _2DMAT_H_
#define _2DMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class StructProp;

class ElaLinIsoMat2D : public NLMaterial
{
   protected:
     double rho, E, nu, t, Tref, alpha;

   public:
     ElaLinIsoMat2D(StructProp *p);
     ElaLinIsoMat2D(double _rho, double _E, double _nu, double _t, double _Tref, double _alpha);

     int getNumStates() { return 0; }

     void getStress(Tensor *stress, Tensor &strain, double*, double temp);

     void transformStress(Tensor &stress, Tensor &gradU, Tensor_d0s2_Ss12 &S);

     void getTangentMaterial(Tensor *tm, Tensor &strain, double*, double temp);

     void getElasticity(Tensor *tm) {};

     void updateStates(Tensor &en, Tensor &enp, double *state, double temp) {};

     void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double temp);
     
     void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                    double *staten, double *statenp, double temp, double dt=0);

     void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                    double *staten, double *statenp, double temp, double dt=0);

     void initStates(double *) {};

     GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator();

     double getThickness() { return t; }
};

// same equation as ElaLinIsoMat2D but with different Green-Lagrange strain evaluator
// (also known as St. Venant-Kirchhoff hyperelastic material
class StVenantKirchhoffMat2D : public ElaLinIsoMat2D
{
  public:
    StVenantKirchhoffMat2D(StructProp *p) : ElaLinIsoMat2D(p) {}
    StVenantKirchhoffMat2D(double rho, double E, double nu, double t, double Tref, double alpha) : ElaLinIsoMat2D(rho, E, nu, t, Tref, alpha) {}

    GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator();
};


#endif
