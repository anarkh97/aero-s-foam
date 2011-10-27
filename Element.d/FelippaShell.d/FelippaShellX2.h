#ifndef _FELIPPASHELLX2_H_
#define _FELIPPASHELLX2_H_

#include <Element.d/SuperElement.h>

class FelippaShellX2 : public SuperElement
{
  public:
    FelippaShellX2(int *nodenums);

    Element* clone();
    int  getTopNumber();
    bool isShell() { return true; }

    // aero functions
    void computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs=0);
    void getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, double *res, GeomState *gs=0);
    bool hasRot() { return true; }

    // Routines for the decomposer
    PrioInfo examine(int sub, MultiFront *);

};

#endif
