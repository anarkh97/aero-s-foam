#ifndef _FOURNODESHELL_H_
#define _FOURNODESHELL_H_

#include <Element.d/SuperElement.h>

class FourNodeShell : public SuperElement 
{
  public:
    FourNodeShell(int *nodenums);

    Element* clone();
    int  getTopNumber();
    bool isShell() { return true; }

    // aero functions
    void computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs=0);
    void getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, double *res, GeomState *gs=0);
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront *mf);

};

#endif
