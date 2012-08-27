#ifndef _DOTTYPE2CFE_H_
#define _DOTTYPE2CFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/DotType2ConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class DotType2ConstraintElement : public ConstraintFunctionElement<DotType2ConstraintFunction>
{
  protected:
    double (*C0)[3]; // initial frame (axes stored row-wise)
    int axis;

  public:
    DotType2ConstraintElement(int*, int, int=2);
    ~DotType2ConstraintElement();
    void setFrame(EFrame *);
    void buildFrame(CoordSet&);
    static const DofSet NODALDOFS[2];

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,6,1>& sconst, Eigen::Array<int,0,1>&, GeomState *gs = NULL);
};

#endif
