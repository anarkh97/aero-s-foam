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
    PressureBCond *pbc;

  public:
    RigidFourNodeShell(int*);
    int getTopNumber() override { return 188; }
    bool isRigidElement() const override { return true; }
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront *mf) override;

    int getMassType() const override { return 0; }
    FullSquareMatrix massMatrix(const CoordSet&, double* mel, int cmflg = 1) const;
    double getMass(const CoordSet& cs) const override;
    void getGravityForce(CoordSet&, double* gravity, Vector&, int gravflg,
                         GeomState *gs);

    void setPressure(PressureBCond *_pbc) { pbc = _pbc; }
    PressureBCond* getPressure() { return pbc; }
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState* gs = 0, int cflg = 0, double t = 0);

    void computeDisp(CoordSet&, State&, const InterpPoint&, double*,
                     GeomState*);
    void getFlLoad(CoordSet&, const InterpPoint&, double*, double *,
                   GeomState* = 0);
};

#endif

