#ifndef _RIGIDTWONODETRUSS_H_
#define _RIGIDTWONODETRUSS_H_

#include <Element.d/Joint.d/BuildingBlocks.d/ConstantDistanceConstraint.h>

class RigidTwoNodeTruss : public ConstantDistanceConstraint
{
  public:
    RigidTwoNodeTruss(int*);
    int getTopNumber() { return 101; }
    bool isRigidElement() { return true; }
    bool isSafe() { return false; }
    PrioInfo examine(int sub, MultiFront*);

    int getMassType() { return 2; } // both consistent and lumped
    FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
    double getMass(CoordSet&);
    void getGravityForce(CoordSet&, double *g, Vector& f, int gravflg,
                         GeomState *gs);
};

#endif
