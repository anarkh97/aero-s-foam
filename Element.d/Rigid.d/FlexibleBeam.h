#ifndef _FLEXIBLEBEAM_H_
#define _FLEXIBLEBEAM_H_

#include <Element.d/SuperElement.h>

class FlexibleBeam : public SuperElement
{
    EFrame *elemframe;
    double c0[3][3];
    double l0;
  public:
    FlexibleBeam(int*);
    int getTopNumber() { return 106; }
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);

    void setProp(StructProp *p, bool myProp);
    void buildFrame(CoordSet&);
    int getMassType() { return 2; }
    FullSquareMatrix massMatrix(CoordSet &cs, double *mel, int cmflg = 1);
    double getMass(CoordSet& cs);
    void getGravityForce(CoordSet&, double *g, Vector &f, int gravflg, GeomState *gs);

  private:
    void getLength(CoordSet&, double &length);
    void updTransMatrix(CoordSet&, GeomState *gs, double t[3][3], double &len);
};

#endif
