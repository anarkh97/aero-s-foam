#ifndef _ELALINISOMAT_H_
#define _ELALINISOMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class StructProp;
class Tensor_d0s4_Ss12s34;

// This material and those derived from it can now be either isotropic or anisotropic
class ElaLinIsoMat : public NLMaterial
{
  protected:
    // isotropic material properties
    double rho, E, nu, alpha;
    // anisotropic material properties
    Tensor_d0s4_Ss12s34 *m_tm;
    double alphas[6];
    // reference temperature
    double Tref;

  public:
    ElaLinIsoMat(StructProp *p);
    ElaLinIsoMat(double _rho, double _E, double _nu, double _Tref, double _alpha);
    ElaLinIsoMat(double _rho, double C[6][6], double _Tref, double _alpha);
    ElaLinIsoMat(double _rho, double C[6][6], double _Tref, double _alphas[6]);
    ~ElaLinIsoMat();

    int getNumStates() { return 0; }

    void getStress(Tensor *stress, Tensor &strain, double*, double temp);

    void getTangentMaterial(Tensor *tm, Tensor &strain, double*);

    void getElasticity(Tensor *tm) {};

    void updateStates(Tensor &en, Tensor &enp, double *state, double temp) {};

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double temp);
     
    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp);

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double temp);

    void initStates(double *){};

    double getDensity() { return rho; }

    StrainEvaluator * getStrainEvaluator();

    double getStrainEnergyDensity(Tensor &enp, double *statenp, double temp);

    void print(std::ostream &out) const {
      out << "Linear "; if(!m_tm) out << rho << " " << E << " " << nu;
    }

    NLMaterial * clone() const;

    void setTangentMaterial(double C[6][6]);
};

// same equation as ElaLinIsoMat but with different Green-Lagrange strain evaluator
// (also known as St. Venant-Kirchhoff hyperelastic material
class StVenantKirchhoffMat : public ElaLinIsoMat
{
  public:
    StVenantKirchhoffMat(StructProp *p) : ElaLinIsoMat(p) {}
    StVenantKirchhoffMat(double rho, double E, double nu, double Tref, double alpha) : ElaLinIsoMat(rho, E, nu, Tref, alpha) {}
    StVenantKirchhoffMat(double rho, double C[6][6], double Tref, double alpha) : ElaLinIsoMat(rho, C, Tref, alpha) {}
    StVenantKirchhoffMat(double rho, double C[6][6], double Tref, double alphas[6]) : ElaLinIsoMat(rho, C, Tref, alphas) {}

    StrainEvaluator * getStrainEvaluator();
    void print(std::ostream &out) const {
      out << "StVenantKirchhoff "; if(!m_tm) out << rho << " " << E << " " << nu;
    }
    NLMaterial * clone() const;
};

class HenckyMat : public ElaLinIsoMat
{
  public:
    HenckyMat(StructProp *p) : ElaLinIsoMat(p) {}
    HenckyMat(double rho, double E, double nu, double Tref, double alpha) : ElaLinIsoMat(rho, E, nu, Tref, alpha) {}
    HenckyMat(double rho, double C[6][6], double Tref, double alpha) : ElaLinIsoMat(rho, C, Tref, alpha) {}
    HenckyMat(double rho, double C[6][6], double Tref, double alphas[6]) : ElaLinIsoMat(rho, C, Tref, alphas) {}

    StrainEvaluator * getStrainEvaluator();
    void print(std::ostream &out) const {
      out << "HenckyElastic "; if(!m_tm) out << rho << " " << E << " " << nu;
    }
    NLMaterial * clone() const;
};

#endif
