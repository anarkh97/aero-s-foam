#ifndef _SOLIDELEMENTTEMPLATE_H_
#define _SOLIDELEMENTTEMPLATE_H_

#include <Element.d/NonLinearity.d/ShapeFunction.h>
#include <Element.d/NonLinearity.d/GaussIntgElem.h>

class NLMaterial;
class ShapeFunction;
class StrainEvaluator;

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes>
class AutoShapeFunction : public ShapeFunction
{
  public:
    AutoShapeFunction() : ShapeFunction(3*NumberOfNodes) {}
    void getLocalDerivatives(Tensor *localDerivatives, double xi[3]);
    void getValues(Tensor *val, double xi[3]) {}
    Tensor *getValInstance() { return 0; }
};

template<template <typename S> class ShapeFunctionTemplate, int NumberOfNodes, int NumIntgPts>
class SolidElementTemplate : public GaussIntgElement
{
  private:
    static const double nodeRefCoords[NumberOfNodes][3];
    static const AutoShapeFunction<ShapeFunctionTemplate,NumberOfNodes> shapeFunction;

  protected:
    int n[NumberOfNodes];
    NLMaterial *material;

    int getNumGaussPoints();
    void getGaussPointAndWeight(int i, double *point, double &weight);
    void getLocalNodalCoords(int i, double *coords);
    ShapeFunction *getShapeFunction();
    StrainEvaluator *getStrainEvaluator();
    NLMaterial *getMaterial();
    void getNodeRefCoords(double (*nodeRefCoords)[3]);

  public:
    SolidElementTemplate(int *nd);
    int numNodes();
    int numDofs();
    void renum(int *);
    void renum(EleRenumMap&);
    void markDofs(DofSetArray &);
    int* dofs(DofSetArray &, int *p=0);
    int* nodes(int * = 0);
    void setMaterial(NLMaterial *);
};

#ifdef _TEMPLATE_FIX_
#include <Element.d/NonLinearity.d/SolidElementTemplate.C>
#endif

#endif
