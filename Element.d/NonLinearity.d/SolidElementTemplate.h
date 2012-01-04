#ifndef _SOLIDELEMENTTEMPLATE_H_
#define _SOLIDELEMENTTEMPLATE_H_

#include <Element.d/NonLinearity.d/ShapeFunction.h>
#include <Element.d/NonLinearity.d/GaussIntgElem.h>
#include <Element.d/Meta.d/ShapeFunctions.h>

class NLMaterial;
class ShapeFunction;
class StrainEvaluator;

//notes:
// the integration points and weights should be static members of an integration template class
// there could be two different integration schemes with the same number of integration points
// there could be two different shape functions with the same number of nodes and the same region so we need some addition mechanism for classifying shape functions and elements
// one possibility is to use template template variables:
//template<typename template<RegionType Region, int NumberOfNodes, class Scalar> ShapeFunctions> class ShapeFunctionTemplate
// another possibility is to have another template parameter
//template<RegionType Region, int NumberOfNodes, ShapeFunctionAttribute Attribute> class ShapeFunctionTemplate, where attribute is for example Serendipity


template<RegionType Region, int NumberOfNodes>
class ShapeFunctionTemplate : public ShapeFunction
{
  public:
    ShapeFunctionTemplate() : ShapeFunction(3*NumberOfNodes) {}
    void getLocalDerivatives(Tensor *localDerivatives, double xi[3]);
    void getValues(Tensor *val, double xi[3]) {}
    Tensor *getValInstance() { return 0; }
};


template<RegionType Region, int NumberOfNodes, int NumIntgPts>
class SolidElementTemplate : public GaussIntgElement
{
  private:
    static const double nodeRefCoords[NumberOfNodes][3];
    static const ShapeFunctionTemplate<Region,NumberOfNodes> shapeFunction;

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
    void markDofs(DofSetArray &);
    int* dofs(DofSetArray &, int *p=0);
    int* nodes(int * = 0);
    void setMaterial(NLMaterial *);
};

#ifdef _TEMPLATE_FIX_
#include <Element.d/NonLinearity.d/SolidElementTemplate.C>
#endif

#endif
