#ifndef _THERM4NOSHELL_H_
#define _THERM4NOSHELL_H_

#include <Element.d/SuperElement.h>
#include <Element.d/Element.h>

class Therm4NoShell : public SuperElement 
{
  public:
    Therm4NoShell(int *nodenums);

    Element* clone() override;
    int getTopNumber() override;
    PrioInfo examine(int sub, MultiFront *) override;
    bool hasRot() const{return true;}
    // aero functions
    void computeTemp(CoordSet &cs, State &state, double[2], double *res);
    void getFlFlux(double[2], double *flF, double *res);
};

#endif
