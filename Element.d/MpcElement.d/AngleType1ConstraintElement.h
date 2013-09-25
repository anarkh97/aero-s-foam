#ifndef _ANGLETYPE1CONSTRAINTELEMENT_H_
#define _ANGLETYPE1CONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/AngleType1ConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class AngleType1ConstraintElement : public ConstraintFunctionElement<Simo::AngleType1ConstraintFunction>
{
  protected:
    double (*C0)[3]; // initial frame (axes stored row-wise)
    int axis1, axis2;
    double offset;

  public:
    AngleType1ConstraintElement(int*, int, int, double = M_PI/2); 
    ~AngleType1ConstraintElement();
    void buildFrame(CoordSet&);
    void setFrame(EFrame *);

  protected:
    void getConstants(CoordSet&, Eigen::Array<double,7,1>& sconst, Eigen::Array<int,0,1>&, GeomState *gs = NULL);
};

#endif
