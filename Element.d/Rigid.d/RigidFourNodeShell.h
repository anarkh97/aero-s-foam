#ifndef _RIGIDFOURNODESHELL_H_
#define _RIGIDFOURNODESHELL_H_

#include <Element.d/SuperElement.h>

class GeomState;
class MultiFront;
class NLMaterial;
class ExpMat;

class RigidFourNodeShell : public SuperElement
{
    ExpMat *expmat;

  public:
    RigidFourNodeShell(int*);
    int getTopNumber() { return 188; }
    bool isRigidElement() { return true; }
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront *mf);

    // EXPERIMENTAL
    void setPressure(double, MFTTData* = 0, bool = false);
    double getPressure();
    FullSquareMatrix massMatrix(CoordSet&, double* mel, int cmflg = 1);
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState* gs = 0, int cflg = 0, double t = 0);

    void computeDisp(CoordSet&, State&, const InterpPoint&, double*,
                     GeomState*);
    void getFlLoad(CoordSet&, const InterpPoint&, double*, double *,
                   GeomState* = 0);
};

#endif

