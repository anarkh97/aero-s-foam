#ifndef _NLMATERIAL_H_
#define _NLMATERIAL_H_

#include <stdexcept>
#include <iostream>

class StrainEvaluator;
class Tensor;
template <typename Tensor> class GenStrainEvaluator;
template <int n> class TwoDTensorTypes;
class MFTTData;

class NLMaterial
{
   public:
     NLMaterial() {}  

     virtual ~NLMaterial() {}

     virtual int getNumStates() = 0;  

     virtual void getTangentMaterial(Tensor *tm, Tensor &strain, double *state, double temp) = 0;

     virtual void getElasticity(Tensor *tm) = 0;

     virtual void getStress(Tensor *stress, Tensor &strain, double *state, double temp) = 0;

     virtual void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double *state, double temp) = 0;

     virtual void updateStates(Tensor& en, Tensor& enp, double *state, double temp) = 0;

     virtual void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                            double *staten, double *statenp, double temp) = 0;

     virtual void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                            double *staten, double *statenp, double temp) = 0;

     virtual void initStates(double *) = 0;

     virtual double getDensity() { return 0; }

     virtual StrainEvaluator * getStrainEvaluator() { return NULL; } // return default strain evaluator 

     virtual GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator() { return NULL; }

     virtual double getEquivPlasticStrain(double *statenp) { return 0; }

     virtual bool getBackStress(double *statenp, Tensor *backstress) { return false; }

     virtual bool getPlasticStrain(double *statenp, Tensor *plasticstrain) { return false; }

     virtual double getDamage(double *statenp) { return 0; }

     virtual double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp) {
       std::cerr << "material law does not implement getStrainEnergyDensity function\n";
       return 0;
     }

     virtual double getDissipatedEnergy(double *statenp) { return 0; }

     virtual double getThickness() { return 0; }

     virtual double getPosdefifyTol() { return -1; }

     virtual void print(std::ostream &out) const {
       throw std::range_error("material law does not implement print function");
     }

     virtual NLMaterial * clone() const {
       std::cerr << "material law does not implement clone function\n";
       return 0;
     }

     virtual void setTangentMaterial(double C[6][6]) {
       std::cerr << "material law does not implement setTangentMaterial function\n";
     }

     virtual void setThermalExpansionCoef(double alpha[6]) {
       std::cerr << "material law does not implement setThermalExpansionCoef function\n";
     }

     virtual void setTDProps(MFTTData *ymtt, MFTTData *ctett) {};
};

#endif
