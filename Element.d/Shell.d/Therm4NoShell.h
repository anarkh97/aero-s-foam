#ifndef _THERM4NOSHELL_H_
#define _THERM4NOSHELL_H_

#include <Element.d/SuperElement.h>
#include <Element.d/Element.h>

class Therm4NoShell : public SuperElement 
{
  public:
    Therm4NoShell(int *nodenums);

    Element* clone();
    int  getTopNumber();
    //bool isShell() { return true; }
    //int numNodes(void) { return 4; };
    //int numDofs(void)  { return 4; };
    PrioInfo examine(int sub, MultiFront *);
    bool hasRot(){return true;}
    // aero functions
    void computeTemp(CoordSet &cs, State &state, double[2], double *res);
    void getFlFlux(double[2], double *flF, double *res);
};

#endif
