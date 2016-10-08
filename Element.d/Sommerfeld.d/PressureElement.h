#ifndef _PRESSUREELEMENT_H_
#define _PRESSUREELEMENT_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class DofSet;
class GeomState;

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree>
class PressureElement : public SommerElement
{
  protected:
    int nNodes;               // number of nodes
    int *nn;                  // node numbers
    std::vector<BCond> terms;
    int nterms;
    PressureBCond *pbc;

    void addTerms(DofSet);

  public:
    PressureElement(int*, PressureBCond*);
    ~PressureElement();

    void renum(int*);
    void renum(EleRenumMap&);

    int numNodes();
    int* nodes(int* = 0); 
    int getNode(int nd) { return nn[nd]; }
    int* getNodes() { return nn; }

    int numDofs();
    int* dofs(DofSetArray&, int* = 0); 
    void markDofs(DofSetArray&);

    int findAndSetEle(CoordSet& cs, Elemset &eset, Connectivity *nodeToEle, int *eleTouch, int *eleCount, int myNum,
                      int it = 0);

    PressureBCond* getPressure() { return pbc; }
    void neumVector(CoordSet&, Vector&, int pflag = 0, GeomState* = 0, double t = 0);
    void neumVectorJacobian(CoordSet&, FullSquareMatrix&, int pflag = 0, GeomState* = 0, double t = 0);
    FullSquareMatrix sommerMatrix(CoordSet&, double *);
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/Sommerfeld.d/PressureElement.C>
#endif

#endif
