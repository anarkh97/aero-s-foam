#ifndef _GENRIGIDMPCBEAM_H_
#define _GENRIGIDMPCBEAM_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>
#include <Utils.d/dofset.h>

class GenRigidMpcBeam : public Element, public Corotator
{
    // general rigid element defining between 1 and 8 constraints between 2 nodes
    // xyz translation, xyz rotation, temperature and pressure
    int nn[2];
    double len0;
    int numcdofs;
    char cdofs[9];
    LMPCons* mpc[8];
    double lambda[8];
    int glMpcNb[8]; 
    DofSet activeDofs;

  public:
    GenRigidMpcBeam(int *);
    GenRigidMpcBeam(int* _nn, int _numcdofs, char* _cdofs); // used by RBE2 superelement
    ~GenRigidMpcBeam() { /* nothing to delete */ };

    virtual void computeMPCs(CoordSet &cs, int &lmpcnum);
    virtual void updateMPCs(GeomState &gState);
    bool isRigidMpcElement() { return true; }
    void setMpcForces(double *mpcForces) { for(int i=0; i<numcdofs; ++i) lambda[i] = mpcForces[glMpcNb[i]]; }

    Corotator *getCorotator(CoordSet &, double *, int, int) { return this; }
    void setProp(StructProp *p) { };

    void renum(int *);

    FullSquareMatrix stiffness(CoordSet &, double *, int = 1);
    FullSquareMatrix massMatrix(CoordSet &, double *, int = 1);

    void markDofs(DofSetArray &);
    int* dofs(DofSetArray &, int * = 0);
    int numDofs();
    int numNodes();
    int* nodes(int * = 0);

    void getStiffAndForce(GeomState &, CoordSet &, FullSquareMatrix &, double *);
    int  getTopNumber() { return 101; }
    bool isSafe() { return true; }

};
#endif
