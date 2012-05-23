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
    std::map<int,std::vector<int> > rotation_indices;
    std::map<int,std::vector<double> > rotation_coefs;

    void addTerms(DofSet);
    void addTerms(DofSet*);

  public:
    MpcElement(int, DofSet, int*);
    MpcElement(int, DofSet*, int*);
    MpcElement(LMPCons *mpc, bool nlflag);
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
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double);
    void getResidualCorrection(GeomState& c1, double* r);
    double getElementEnergy(GeomState&, CoordSet&) { return 0; }

    virtual void update(GeomState&, CoordSet&, double);
    virtual void getHessian(GeomState&, CoordSet&, FullSquareMatrix&, double t);
    virtual double getVelocityConstraintRhs(GeomState *gState, CoordSet& cs, double t);
    virtual double getAccelerationConstraintRhs(GeomState *gState, CoordSet& cs, double t);

    PrioInfo examine(int sub, MultiFront *mf);
    bool isSafe() { return false; }

    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0);
    int getTopNumber() { return 101; }

    void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                             int &nlflag) { nlflag = 2; }

    void getNLVonMises(Vector&, Vector& weight,
                       GeomState &, CoordSet &, int);
};
#endif
