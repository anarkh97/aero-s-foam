#ifndef _ELALINISOMAT_H_
#define _ELALINISOMAT_H_
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Utils.d/NodeSpaceArray.h>

class StructProp;

// This material and those derived from it can now be either isotropic or anisotropic
class ElaLinIsoMat : public NLMaterial
{
  protected:
    double rho, E, nu;
    Tensor_d0s4_Ss12s34 *m_tm;

  public:
    ElaLinIsoMat(StructProp *p);
    ElaLinIsoMat(double _rho, double _E, double _nu);
    ElaLinIsoMat(double _rho, double C[6][6]);
    ~ElaLinIsoMat();

    int getNumStates() { return 0; }

    void getStress(Tensor *stress, Tensor &strain, double*);

    void getTangentMaterial(Tensor *tm, Tensor &strain, double*);

    void getElasticity(Tensor *tm) {};

    void updateStates(Tensor en,Tensor enp,double *state) {};

    void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*);
     
    void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double);

    void integrate(Tensor *stress, Tensor &en, Tensor &enp,
                   double *staten, double *statenp, double);

    void initStates(double *){};

    double getDensity() { return rho; }

    StrainEvaluator * getStrainEvaluator();

    double getStrainEnergyDensity(Tensor &enp, double *statenp);

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
    StVenantKirchhoffMat(double rho, double E, double nu) : ElaLinIsoMat(rho, E, nu) {}
    StVenantKirchhoffMat(double rho, double C[6][6]) : ElaLinIsoMat(rho, C) {}

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
    HenckyMat(double rho, double E, double nu) : ElaLinIsoMat(rho, E, nu) {}
    HenckyMat(double rho, double C[6][6]) : ElaLinIsoMat(rho, C) {}

    StrainEvaluator * getStrainEvaluator();
    void print(std::ostream &out) const {
      out << "HenckyElastic "; if(!m_tm) out << rho << " " << E << " " << nu;
    }
    NLMaterial * clone() const;
};


#endif
 
