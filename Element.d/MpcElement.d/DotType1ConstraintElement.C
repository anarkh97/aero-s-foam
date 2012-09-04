#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/DotType1ConstraintElement.h>
#include <Corotational.d/GeomState.h>

DotType1ConstraintElement::DotType1ConstraintElement(int* _nn, int _axis1, int _axis2, double _d0, int _rotdescr)
 : ConstraintFunctionElement<DotType1ConstraintFunction>(2, DofSet::XYZrot, _nn, 0, _rotdescr)
{
  C0 = 0;
  axis1 = _axis1;
  axis2 = _axis2;
  d0 = _d0;
}

DotType1ConstraintElement::~DotType1ConstraintElement()
{
  if(C0) delete [] C0;
}

void
DotType1ConstraintElement::getConstants(CoordSet&, Eigen::Array<double,7,1>& sconst, Eigen::Array<int,0,1>&, GeomState *gs)
{
  // note #1: gs is interpreted as the reference configuration for updated lagrangian and the current configuration for eulerian
  // note #2: gs = NULL means that the reference/updated configuration is identical to the undeformed configuration
  if(rotdescr == 0 || gs == NULL) {
    Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&DotType1ConstraintElement::C0[0][0]);
    sconst << C0(axis1,0), C0(axis1,1), C0(axis1,2), C0(axis2,0), C0(axis2,1), C0(axis2,2), d0;
  }
  else {
    Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&DotType1ConstraintElement::C0[0][0]),
                                                           R1(&(*gs)[nn[0]].R[0][0]),
                                                           R2(&(*gs)[nn[1]].R[0][0]);
    // rotate specified axes to their position in the specified (either reference or current) configuration
    Eigen::Matrix<double,1,3> a0 = C0.row(axis1)*R1.transpose(), b0 = C0.row(axis2)*R2.transpose();
    sconst << a0[0], a0[1], a0[2], b0[0], b0[1], b0[2], d0;
  }
}

void
DotType1ConstraintElement::setFrame(EFrame *elemframe) 
{ 
  C0 = new double[3][3];
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) 
      C0[i][j] = (*elemframe)[i][j]; 
}

void 
DotType1ConstraintElement::buildFrame(CoordSet& cs)
{
  if(!C0) {
    cerr << " *** ERROR: element frame is not defined for DotType1 constraint function element\n";
    exit(-1);
  }

  ConstraintFunctionElement<DotType1ConstraintFunction>::buildFrame(cs);
}
#endif
