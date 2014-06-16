#ifndef _EXPMAT_H_
#define _EXPMAT_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>

#include <iostream>
#include <iterator>

class ExpMat : public NLMaterial
{
  public:
    int optctv; // constitutive law (1 for hypoelastic, 3 for elasto viscoplastic, 5 for j2 explicit)
    double ematpro[20]; // Young's modulus, Poisson's ratio, mass density, etc.

    ExpMat(int _optctv, double e1, double e2, double e3, double e4, double e5, double e6,
           double e7, double e8, double e9, double e10, double e11, double e12, double e13,
           double e14, double e15, double e16, double e17, double e18, double e19, double e20)
      { optctv = _optctv;
        ematpro[0] = e1; ematpro[1] = e2; ematpro[2] = e3; ematpro[3] = e4; ematpro[4] = e5;
        ematpro[5] = e6; ematpro[6] = e7; ematpro[7] = e8; ematpro[8] = e9; ematpro[9] = e10;
        ematpro[10] = e11; ematpro[11] = e12; ematpro[12] = e13; ematpro[13] = e14; ematpro[14] = e15;
        ematpro[15] = e16; ematpro[16] = e17; ematpro[17] = e18; ematpro[18] = e19; ematpro[19] = e20; }

    int getNumStates() { return 0; }

    void getStress(Tensor *stress, Tensor &strain, double*, double)
      { std::cerr << "ExpMat::getStress is not implemented\n"; }

    void getTangentMaterial(Tensor *tm, Tensor &strain, double*)
      { std::cerr << "ExpMat::getTangentMaterial is not implemented\n"; }

    void getElasticity(Tensor *tm)
      { std::cerr << "ExpMat::getElasticity is not implemented\n"; }

    void updateStates(Tensor &en, Tensor &enp, double *state, double) {}

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double)
      { std::cerr << "ExpMat::getStressAndTangentMaterial is not implemented\n"; }
     
    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double)
      { std::cerr << "ExpMat::integrate is not implemented\n"; }

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double)
      { std::cerr << "ExpMat::integrate is not implemented\n"; }

    void initStates(double *) {}

    double getDensity() { return ematpro[2]; } 

    StrainEvaluator * getStrainEvaluator()
      { std::cerr << "ExpMat::getStrainEvaluator is not implemented\n"; return NULL; }

    void print(std::ostream &out) const {
      std::string type;
      switch (optctv) {
        case 1: type = "HypoElastic"; break;
        case 5: type = "J2Plasticity"; break;
        case 6: type = "KK1"; break;
        case 7: type = "KK2"; break;
        default: throw std::range_error("Unknown material law type");
      }
      out << type << " ";
      std::copy(&ematpro[0], &ematpro[20], std::ostream_iterator<double>(out, " "));
    }
};

#endif
