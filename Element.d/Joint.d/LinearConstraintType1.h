#ifndef _LINEARCONSTRAINTTYPE1_H_
#define _LINEARCONSTRAINTTYPE1_H_

#include <Element.d/MpcElement.d/MpcElement.h>

class LinearConstraintType1 : public MpcElement
{
  public:
    LinearConstraintType1(int*, DofSet*);
    void buildFrame(CoordSet&);
    void update(GeomState& gState, CoordSet& cs, double);
};

#endif
