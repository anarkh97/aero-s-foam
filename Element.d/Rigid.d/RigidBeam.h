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
    //LMPCons** getMPCs();

    void buildFrame(CoordSet&);
    int getMassType() { return 2; } // both consistent and lumped
    FullSquareMatrix massMatrix(CoordSet &cs, double *mel, int cmflg = 1);
    double getMass(CoordSet& cs);
    void getGravityForce(CoordSet&, double *g, Vector &f, int gravflg, GeomState *gs);
    double computeStabilityTimeStep(FullSquareMatrix &K, FullSquareMatrix &M, CoordSet &cs, GeomState *gs,
                                    double stable_tol, int stable_maxit);

  private:
    void getLength(CoordSet&, double &length);
    void updTransMatrix(CoordSet&, GeomState *gs, double t[3][3], double &len);
};

#endif
