#ifndef _NLHEXAHEDRAL_H_
#define _NLHEXAHEDRAL_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/GaussIntgElem.h>
#include <Element.d/NonLinearity.d/3DShapeFunction.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>

class HexahedralShapeFunction : public ShapeFunction
{
  public:
    HexahedralShapeFunction() : ShapeFunction(24) {}
    void getLocalDerivatives(Tensor *localDerivatives, double xi[3]);
    void getValues(Tensor *gradU, Tensor *dgradUdqk, double *jac, double xi[3]) {}
    Tensor *getValInstance() { return 0; }
};

class NLHexahedral : public GaussIntgElement
{
    int n[8];
    NLMaterial *material;
    int strainMeasure;

  protected:
    int getNumGaussPoints();
    void getGaussPointAndWeight(int i, double *point, double &weight);
    ShapeFunction *getShapeFunction();
    StrainEvaluator *getStrainEvaluator();
    NLMaterial *getMaterial();

  public:
    NLHexahedral(int *nd, int = -1);
    int numNodes() { return 8; }
    int numDofs() { return 24; }
    PrioInfo examine(int sub, MultiFront *); // dec
    void renum(int *);
    void   markDofs(DofSetArray &);
    int*   dofs(DofSetArray &, int *p=0);
    int*   nodes(int * = 0);
    //void updateStates(Node *nodes, double *states, double *un, double *unp) {}
    void setProp(StructProp *);
    void setMaterial(NLMaterial *);
    int getTopNumber();
};

#endif
