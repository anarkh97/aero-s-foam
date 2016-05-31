#ifndef _BRITTLEFRACTURETB_H_
#define _BRITTLEFRACTURETB_H_

class Tensor;

template<typename BaseMaterial>
class BrittleFractureTB : public BaseMaterial
{
    double maxprs;
    double exponent; 
    double Kf; //stress impulse

  public:
    BrittleFractureTB(StructProp *p) : BaseMaterial(p) {}
    BrittleFractureTB(double rho, double E, double nu, double Tref, double alpha, 
                      double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(rho, E, nu, Tref, alpha), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}
    BrittleFractureTB(double rho, double C[6][6], double Tref, double alphas[6], 
                      double _maxprs, double _exponent, double _Kf)
     : BaseMaterial(rho, C, Tref, alphas), maxprs(_maxprs), exponent(_exponent), Kf(_Kf) {}

    int getNumStates();

    void initStates(double *);

    void getStress(Tensor *stress, Tensor &strain, double*, double temp);

    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, double dt=0);

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, double dt=0);

    double getDamage(double *statenp);

    void print(std::ostream &out) const;
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/NonLinearity.d/BrittleFractureTB.C>
#endif

#endif
