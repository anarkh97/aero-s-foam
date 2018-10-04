#ifndef _FABRICMAT_H_
#define _FABRICMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class StructProp;

class FabricMat : public NLMaterial
{
   public:
     enum StrainMeasure { INFINTESIMAL = 0, GREEN_LAGRANGE = 1 };

   protected:
     int x_ymst_id, y_ymst_id;
     double rho, Gxy, nuxy, nuyx, t, Tref, alpha;
     enum StrainMeasure strain_measure;
     MFTTData *x_ymst, *y_ymst;

   public:
     FabricMat(StructProp *p);
     FabricMat(double _rho, int _x_ymst_id, int _y_ymst_id, double _Gxy, double _nuxy, double _nyx, double _t,
               double _Tref, double _alpha, StrainMeasure _strain_measure);

     int getNumStates() { return 0; }

     void getStress(Tensor *stress, Tensor &strain, double*, double temp);

     void getTangentMaterial(Tensor *tm, Tensor &strain, double*, double temp);

     void getElasticity(Tensor *tm) {};

     void updateStates(Tensor &en, Tensor &enp, double *state, double temp) {};

     void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double temp);
     
     void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                    double *staten, double *statenp, double temp,
                    Tensor *cache, double dt=0);

     void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                    double *staten, double *statenp, double temp,
                    Tensor *cache, double dt=0);

     void initStates(double *) {};

     GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator();

     double getDensity() { return rho; }

     double getThickness() { return t; }

     double getReferenceTemperature() { return Tref; }

     void setEDProps(MFTTData *ymst);
};

#endif
