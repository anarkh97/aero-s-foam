#ifndef _RIGIDEIGHTNODEBRICK_H_
#define _RIGIDEIGHTNODEBRICK_H_

#include <Element.d/SuperElement.h>

class RigidEightNodeBrick : public SuperElement
{
  public:
    RigidEightNodeBrick(int*);
    int getTopNumber() override { return 117; }
    bool isRigidElement() { return true; }
    bool isSafe() { return true; }
    PrioInfo examine(int sub, MultiFront*);

    FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);
    double getMass(CoordSet& cs);
    void getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
                         GeomState *gs);
};

#endif

