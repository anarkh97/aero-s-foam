#ifndef _COMPO4NODESHELL_H_
#define _COMPO4NODESHELL_H_

#include <Element.d/SuperElement.h>

class Compo4NodeShell : public SuperElement
{
  public:
    Compo4NodeShell(int *nodenums);

    Element* clone() override;
    int getTopNumber() override;
    bool isShell() const override { return true; }

    // aero functions
    void computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs=0);
    void getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, double *res, GeomState *gs=0);
    bool hasRot() { return true; }

    // Routines for the decomposer
    PrioInfo examine(int sub, MultiFront *) override;
    int nDecFaces() const override { return 1; }
    int getDecFace(int iFace, int *fn) { for(int i=0; i<4; i++) fn[i] = nn[i]; return 4; }

    int getFace(int iFace, int *fn) { return getDecFace(iFace,fn); }
};

#endif
