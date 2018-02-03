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

    FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg=1) const;
    double getMass(const CoordSet& cs) const;
    void getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
                         GeomState *gs);
};

#endif

