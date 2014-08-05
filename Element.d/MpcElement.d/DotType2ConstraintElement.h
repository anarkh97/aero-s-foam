#ifndef _DOTTYPE2CONSTRAINTELEMENT_H_
#define _DOTTYPE2CONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/DotType2ConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class DotType2ConstraintElement : public ConstraintFunctionElement<Simo::DotType2ConstraintFunction>
{
  protected:
    double (*C0)[3]; // initial frame (axes stored row-wise)
    int axis;
    double d0;

  public:
    DotType2ConstraintElement(int*, int);
    ~DotType2ConstraintElement();
    void setFrame(EFrame *);
    void buildFrame(CoordSet&);
    static const DofSet NODALDOFS[2];
    void setConstantTerm(double _d0) { d0 = _d0; }
    double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double);
    double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&, double);

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,7,1>& sconst, Eigen::Array<int,0,1>&, GeomState *gs = NULL);
};

#endif
