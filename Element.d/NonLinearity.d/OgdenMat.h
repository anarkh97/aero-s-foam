#ifndef _OGDENMAT_H_
#define _OGDENMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

// This material and those derived from it can now be either isotropic or anisotropic
class OgdenMat : public NLMaterial
{
  protected:
    // isotropic material properties
    double rho; // density
    double mu[3], alpha[3]; // material properties characterizing distortional response
    int N; // number of terms in the Ogden series
    double c, d; // material properties characterizing volumetric response
    int vol; // identifies the particular form of the volumetric stored energy function

  public:
    OgdenMat(double _rho, double _mu1, double _alpha1, double _c, double _d, int _vol=1);
    OgdenMat(double _rho, double _mu1, double _mu2, double _alpha1, double _alpha2, double _c, double _d, int _vol=1);
    OgdenMat(double _rho, double _mu1, double _mu2, double _mu3, double _alpha1, double _alpha2, double _alpha3, double _c, double _d, int _vol=1);

    int getNumStates() { return 0; }

    void getStress(Tensor *stress, Tensor &strain, double*, double temp);

    void getTangentMaterial(Tensor *tm, Tensor &strain, double*, double temp);

    void getElasticity(Tensor *tm) {};

    void updateStates(Tensor &en, Tensor &enp, double *state, double temp) {};

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double temp);
     
    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, double dt=0);

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp, double dt=0);

    void initStates(double *) {};

    double getDensity() { return rho; }

    StrainEvaluator * getStrainEvaluator();

    double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp);

    void print(std::ostream &out) const {
      switch(N) {
        case 1:
          out << "Ogden " << rho << " " << mu[0] << " " << alpha[0] << " " << c << " " << d;
          break;
        case 2:
          out << "Ogden " << rho << " " << mu[0] << " " << mu[1] << " " << alpha[0] << " " << alpha[1] << " " << c << " " << d;
          break;
        case 3:
          out << "Ogden " << rho << " " << mu[0] << " " << mu[1] << " " << mu[2] << " " << alpha[0] << " " << alpha[1] << " " << alpha[2] << " " 
              << c << " " << d;
          break;
      }
    }

    NLMaterial * clone() const;

    void getMaterialConstants(std::vector<double> &c);
};

#endif
