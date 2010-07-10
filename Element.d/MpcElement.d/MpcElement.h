#ifndef _MPCELEMENT_H_
#define _MPCELEMENT_H_

#include <Driver.d/Mpc.h>
#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <map>

class DofSet;

class MpcElement : public Element, public Corotator, public LMPCons
{
  protected:
    int nNodes;              // number of nodes (not including "internal node")
    int *nn;                 // node numbers
    // for direct elimination / coordinate split / state space  set prop->lagrangeMult to false and prop->penalty to 0
    // for lagrange multipliers method set prop->lagrangeMult to true and prop->penalty to 0
    // for penalty method set prop->lagrangeMult to false and prop->penalty to some large number
    // for augmented lagrangian method set prop->lagrangeMult to true and prop->penalty to some large number

    void addTerms(DofSet);
    void addTerms(DofSet*);

  public:
    MpcElement(int, DofSet, int*);
    MpcElement(int, DofSet*, int*);
    MpcElement(LMPCons *mpc);
   ~MpcElement();

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

    bool isMpcElement() { return true; }

    Corotator* getCorotator(CoordSet&, double*, int, int);
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*);
    double getElementEnergy(GeomState&, CoordSet&) { }

    virtual void update(GeomState&, CoordSet&);
    virtual void getHessian(GeomState&, CoordSet&, FullSquareMatrix&);

    PrioInfo examine(int sub, MultiFront *mf);
    bool isSafe() { return false; }

    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs, int cflg);
    int getTopNumber() { return 101; }
};
#endif
