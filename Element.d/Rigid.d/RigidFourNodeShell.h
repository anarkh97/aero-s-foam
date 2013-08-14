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
    BlastLoading::BlastData *conwep;

  public:
    RigidFourNodeShell(int*);
    int getTopNumber() { return 188; }
    bool isRigidElement() { return true; }
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront *mf);

    void setPressure(double, MFTTData* = 0, BlastLoading::BlastData* = 0);

    int getMassType() { return 0; }
    FullSquareMatrix massMatrix(CoordSet&, double* mel, int cmflg = 1);
    double getMass(CoordSet& cs);
    void getGravityForce(CoordSet&, double* gravity, Vector&, int gravflg,
                         GeomState *gs);

    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState* gs = 0, int cflg = 0, double t = 0);

    void computeDisp(CoordSet&, State&, const InterpPoint&, double*,
                     GeomState*);
    void getFlLoad(CoordSet&, const InterpPoint&, double*, double *,
                   GeomState* = 0);
};

#endif

