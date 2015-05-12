#ifndef _BOUNDARYELEMENT_H_
#define _BOUNDARYELEMENT_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <vector>

class DofSet;

class BoundaryElement : public Element, public Corotator
{
  protected:
    int nNodes;               // number of nodes
    int *nn;                  // node numbers
    std::vector<BCond> terms;
    int nterms;

    void addTerms(DofSet);
    void addTerms(DofSet*);
    void addTerms(DofSet*, DofSet*);

    std::vector<int> inputs;
    std::vector<int> outputs;

  public:
    BoundaryElement(int, DofSet, int*);
    BoundaryElement(int, DofSet*, int*);
    BoundaryElement(int, DofSet*, DofSet*, int*);
   ~BoundaryElement();

    void renum(int*);
    void renum(EleRenumMap&);

    int numNodes();
    int* nodes(int* = 0);

    int numDofs();
    int* dofs(DofSetArray&, int* = 0);
    void markDofs(DofSetArray&);

    Corotator* getCorotator(CoordSet&, double*, int, int);
    double getElementEnergy(GeomState&, CoordSet&) { return 0; }

    bool isSafe() { return false; }

    int getTopNumber() { return 101; }

    void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                             int &nlflag) { nlflag = 2; }

    void getNLVonMises(Vector&, Vector& weight,
                       GeomState &, CoordSet &, int);

    void getGravityForce(CoordSet&, double*, Vector&, int, GeomState* =0);
};
#endif
