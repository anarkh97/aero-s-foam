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
    int nDofs;               // number of dofs (not including lagrange multiplier)
    int *nDofsPerNode;       // number of degrees of freedom per node (not including "internal node")
    DofSet *nodalDofs;       // active dofs on each node
    int *nn;                 // node numbers
    int mode;                // mode 0: constraints extracted by getMPCs, mode 1: otherwise

  public:
    ConstraintElement(int, DofSet, int*);
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

    Corotator* getCorotator(CoordSet&, double*, int, int);
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*);

    virtual void updateLMPCs(GeomState&, CoordSet&);
    virtual void getHessian(GeomState&, CoordSet&, FullSquareMatrix&);

};
#endif
