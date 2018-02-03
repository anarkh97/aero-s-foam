#ifndef _RIGIDTHREENODESHELL_H_
#define _RIGIDTHREENODESHELL_H_

#include <Element.d/SuperElement.h>

class RigidThreeNodeShell : public SuperElement
{
    PressureBCond *pbc;

  public:
    RigidThreeNodeShell(int*);
    int getTopNumber() override { return 108; }
    bool isRigidElement() { return true; }
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront *mf) override;

    int getMassType() { return 0; }
    FullSquareMatrix massMatrix(CoordSet&, double* mel, int cmflg = 1);
    double           getMass(CoordSet& cs);
    void             getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
                                     GeomState *gs);

    void             computeDisp(CoordSet&, State &, const InterpPoint &,
                                 double*, GeomState *gs);
    void             getFlLoad(CoordSet &, const InterpPoint &,
                               double *flF, double *resF, GeomState *gs=0);

    void setPressure(PressureBCond *_pbc) { pbc = _pbc; }
    PressureBCond* getPressure() { return pbc; }
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0);
};

#endif
