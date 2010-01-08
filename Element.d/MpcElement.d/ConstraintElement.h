#ifndef _CONSTRAINTELEMENT_H_
#define _CONSTRAINTELEMENT_H_

#include <Driver.d/Mpc.h>
#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>

class DofSet;

class ConstraintElement : public Element, public Corotator, public LMPCons
{
  protected:
    int nNodes;              // number of nodes (not including "internal node")
    int *nn;                 // node numbers

  public:
    ConstraintElement(int, DofSet, int*);
    ConstraintElement(int, DofSet*, int*);
    ConstraintElement(LMPCons *mpc);
   ~ConstraintElement();

    int getNumMPCs();
    LMPCons** getMPCs();

    void renum(int*);

    int numNodes();
    int* nodes(int* = 0);

    int numInternalNodes();
    void setInternalNodes(int*);

    int numDofs();
    int* dofs(DofSetArray&, int* = 0);
    void markDofs(DofSetArray&);

    FullSquareMatrix stiffness(CoordSet&, double*, int = 1);

    void getGravityForce(CoordSet&, double*, Vector& f, int, GeomState* = 0) { f.zero(); }

    bool isConstraintElement() { return true; }

    Corotator* getCorotator(CoordSet&, double*, int, int);
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*);

    virtual void update(GeomState&, CoordSet&);
    virtual void getHessian(GeomState&, CoordSet&, FullSquareMatrix&);

};
#endif
