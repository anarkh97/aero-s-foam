#ifndef _PRESSUREELEMENT_H_
#define _PRESSUREELEMENT_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class DofSet;
class GeomState;

template<template <typename S> class VectorValuedFunctionTemplate>
class PressureElement : public SommerElement
{
  protected:
    int nNodes;               // number of nodes
    int *nn;                  // node numbers
    std::vector<BCond> terms;
    int nterms;

    void addTerms(DofSet);

  public:
    PressureElement(int, DofSet, int*);
    PressureElement(int, DofSet*, int*);
    PressureElement(int, DofSet*, DofSet*, int*);
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

    void neumVector(CoordSet&, Vector&, int pflag=0, GeomState* = 0);
    void neumVectorJacobian(CoordSet&, FullSquareMatrix&, int pflag=0, GeomState* = 0);
    FullSquareMatrix sommerMatrix(CoordSet&, double *);

  protected:
    virtual void getConstants(CoordSet&,
                              Eigen::Array<typename VectorValuedFunctionTemplate<double>::ScalarConstantType,
                                           VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> &sconst,
                              Eigen::Array<int,
                                           VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> &iconst) = 0;

   virtual void getInputs(Eigen::Matrix<double,VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q, 
                          CoordSet& c0, GeomState *c1);
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/Sommerfeld.d/PressureElement.C>
#endif

#endif
