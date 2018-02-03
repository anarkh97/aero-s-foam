#ifndef _FELIPPASHELLX2_H_
#define _FELIPPASHELLX2_H_

#ifdef USE_EIGEN3
#include <Element.d/SuperElement.h>

class FelippaShellX2 : public SuperElement
{
  public:
    FelippaShellX2(int *nodenums);

    Element* clone();
    int getTopNumber() override;
    bool isShell() { return true; }

    // aero functions
    void computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs=0);
    void getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, double *res, GeomState *gs=0);
    bool hasRot() { return true; }

    // Routines for the decomposer
    PrioInfo examine(int sub, MultiFront *) override;
    int nDecFaces() { return 1; }
    int getDecFace(int iFace, int *fn) { for(int i=0; i<4; i++) fn[i] = nn[i]; return 4; }

    int getFace(int iFace, int *fn) { return getDecFace(iFace,fn); }
};

#endif
#endif
