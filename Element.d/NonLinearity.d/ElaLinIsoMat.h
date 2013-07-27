#ifndef _ELALINISOMAT_H_
#define _ELALINISOMAT_H_
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Utils.d/NodeSpaceArray.h>

class StructProp;

//Declaration of the material properties
class ElaLinIsoMat : public NLMaterial
{
  protected:
    double rho, E, nu;

  public:
    ElaLinIsoMat(StructProp *p);
    ElaLinIsoMat(double _rho, double _E, double _nu);

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

    void print(std::ostream &out) const {
      out << "Linear " << rho << " " << E << " " << nu;
    }
};

// same equation as ElaLinIsoMat but with different Green-Lagrange strain evaluator
// (also known as St. Venant-Kirchhoff hyperelastic material
class StVenantKirchhoffMat : public ElaLinIsoMat
{
  public:
    StVenantKirchhoffMat(StructProp *p) : ElaLinIsoMat(p) {}
    StVenantKirchhoffMat(double rho, double E, double nu) : ElaLinIsoMat(rho, E, nu) {}

    StrainEvaluator * getStrainEvaluator();
    void print(std::ostream &out) const {
      out << "StVenantKirchhoff " << rho << " " << E << " " << nu;
    }
};

class HenckyMat : public ElaLinIsoMat
{
  public:
    HenckyMat(StructProp *p) : ElaLinIsoMat(p) {}
    HenckyMat(double rho, double E, double nu) : ElaLinIsoMat(rho, E, nu) {}

    StrainEvaluator * getStrainEvaluator();
    void print(std::ostream &out) const {
      out << "HenckyElastic " << rho << " " << E << " " << nu;
    }
};


#endif
 
