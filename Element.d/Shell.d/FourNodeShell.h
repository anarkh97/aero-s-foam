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

    int nDecFaces() { return 1;}
    int getDecFace(int iFace, int *fn) {
          for(int i=0;i<4;i++) fn[i] = nn[i]; return 4; }

    // aero functions
    void computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs=0);
    void getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, double *res, GeomState *gs=0);
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront *mf);

};

#endif
