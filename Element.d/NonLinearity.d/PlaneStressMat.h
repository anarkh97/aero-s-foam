#ifndef _PLANESTRESSMAT_H_
#define _PLANESTRESSMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class StructProp;

template<class BaseMaterial>
class PlaneStressMat : public BaseMaterial
{
     double t;

   public:
     PlaneStressMat(double p1, double p2, double p3, double _t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double _t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double _t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double _t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double _t);
     PlaneStressMat(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double _t);
     PlaneStressMat(double *p, double _t);

     int getNumStates();

     void initStates(double *);

     void getStress(Tensor *stress, Tensor &strain, double*, double temp);

     void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                    double *staten, double *statenp, double temp,
                    Tensor *cache, double dt=0);

     void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                    double *staten, double *statenp, double temp,
                    Tensor *cache, double dt=0);

     void print(std::ostream &out) const;

     GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator();

     double getThickness() { return t; }
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/NonLinearity.d/PlaneStressMat.C>
#endif

#endif
