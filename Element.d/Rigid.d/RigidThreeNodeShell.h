#ifndef _RIGIDTHREENODESHELL_H_
#define _RIGIDTHREENODESHELL_H_

#include <Element.d/SuperElement.h>

class RigidThreeNodeShell : public SuperElement
{
  public:
    RigidThreeNodeShell(int*);
    int getTopNumber() { return 108; }
    bool isRigidElement() { return true; }
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront *mf);

    int getMassType() { return 0; }
    FullSquareMatrix massMatrix(CoordSet&, double* mel, int cmflg = 1);
    double           getMass(CoordSet& cs);
    void             getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
                                     GeomState *gs);
};

#endif

