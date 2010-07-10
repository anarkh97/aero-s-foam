#ifndef _RIGIDBEAM_H_
#define _RIGIDBEAM_H_

#include <Element.d/SuperElement.h>

class RigidBeam : public SuperElement
{
  public:
    RigidBeam(int*);
    int getTopNumber() { return 106; }
    bool isRigidElement() { return true; }
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
    LMPCons** getMPCs();
/*
    // experimental 
    int getMassType();
    FullSquareMatrix massMatrix(CoordSet &cs, double *mel, int cmflg = 1);
    double getMass(CoordSet& cs);
*/
};

#endif
