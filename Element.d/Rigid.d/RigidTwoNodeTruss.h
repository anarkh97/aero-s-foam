#ifndef _RIGIDTWONODETRUSS_H_
#define _RIGIDTWONODETRUSS_H_

#include <Element.d/Joint.d/BuildingBlocks.d/ConstantDistanceConstraint.h>

class RigidTwoNodeTruss : public ConstantDistanceConstraint
{
  public:
    RigidTwoNodeTruss(int*);
    int getTopNumber() override { return 101; }
    bool isRigidElement() const override { return true; }
    bool isSafe() const override { return false; }
    PrioInfo examine(int sub, MultiFront*) override;
};

class RigidTwoNodeTrussWithMass : public ConstantDistanceConstraint
{
  public:
    RigidTwoNodeTrussWithMass(int*);
    int getTopNumber() override { return 101; }
    bool isRigidElement() const override { return true; }
    bool isSafe() const override { return false; }
    PrioInfo examine(int sub, MultiFront*) override;

    int getMassType() const override { return 2; } // both consistent and lumped
    FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg=1) const;
    double getMass(const CoordSet&) const;
    void getGravityForce(CoordSet&, double *g, Vector& f, int gravflg,
                         GeomState *gs);
};

#endif
