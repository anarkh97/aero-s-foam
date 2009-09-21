#ifndef _RIGIDMPCELEMENT_H_
#define _RIGIDMPCELEMENT_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>

class RigidMpcElement : public Element, public Corotator
{
  protected:
    int nNodes;              // number of nodes (not including "internal node")
    int nDofsPerNode;        // number of degrees of freedom per node (not including "internal node")
    int nodalDofs;           // active dofs on each node
    int nMpcs;               // number of constraints
    int *nn, *nn_orig;       // node numbers (renumbered and original)
    LMPCons** mpcs;          // constraints
    int mode;                // mode 0: constraints extracted by getMPCs, mode 1: otherwise

  public:
    RigidMpcElement(int, int, int, int, int*);
   ~RigidMpcElement();

    void setProp(StructProp*);

    virtual void computeMPCs(CoordSet&) = 0;
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

    Corotator* getCorotator(CoordSet&, double*, int, int);
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*);

    bool isRigidMpcElement();

};
#endif
