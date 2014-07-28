#ifndef _RIGIDBEAM_H_
#define _RIGIDBEAM_H_

#include <Element.d/SuperElement.h>

class RigidBeam : public SuperElement
{
    EFrame *elemframe;
    double c0[3][3];
    int variant;
    double length;
  public:
    RigidBeam(int*, int=0);
    int getTopNumber() { return 106; }
    bool isRigidElement() { return true; }
    bool hasRot() { return true; }
    bool isSafe() { return true; }
    PrioInfo examine(int sub, MultiFront*);

    void buildFrame(CoordSet&);

  private:
    void getLength(CoordSet&, double &length);
};

class RigidBeamWithMass : public SuperElement
{
    EFrame *elemframe;
    double c0[3][3];
    int variant;
    double length;
  public:
    RigidBeamWithMass(int*, int=0);
    int getTopNumber() { return 106; }
    bool isRigidElement() { return true; }
    bool hasRot() { return true; }
    bool isSafe() { return true; }
    PrioInfo examine(int sub, MultiFront*);

    void buildFrame(CoordSet&);
    void setProp(StructProp *p, bool _myProp = false);
    int getMassType() { return 0; } // lumped
    double computeStabilityTimeStep(FullSquareMatrix &K, FullSquareMatrix &M, CoordSet &cs, GeomState *gs,
                                    double stable_tol, int stable_maxit);

  private:
    void getLength(CoordSet&, double &length);
};

#endif
