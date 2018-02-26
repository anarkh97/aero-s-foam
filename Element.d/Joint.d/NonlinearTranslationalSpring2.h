#ifndef _NONLINEARTRANSLATIONALSPRING2_H_
#define _NONLINEARTRANSLATIONALSPRING2_H_

#include <Element.d/MpcElement.d/DotType3ConstraintElement.h>

class NonlinearTranslationalSpring2 : public DotType3ConstraintElement
{
    double sp0;

  public:
    NonlinearTranslationalSpring2(int*, int);
    void setProp(StructProp *p, bool _myProp) override;
    void buildFrame(CoordSet&);
    void update(GeomState *refState, GeomState& gState, CoordSet& cs, double);

    bool isSpring() { return true; }
    bool hasRot() const override { return true; }
};

#endif
