#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/DotType2v2ConstraintElement.h>
#include <Corotational.d/GeomState.h>

DotType2v2ConstraintElement::DotType2v2ConstraintElement(int* _nn, int _axis, int _rotdescr)
 : ConstraintFunctionElement<DotType2v2ConstraintFunction>(2, DofSet::XYZdisp | DofSet::XYZrot, _nn, 0, _rotdescr)
{
  axis = _axis;
  C0 = 0;
}

DotType2v2ConstraintElement::~DotType2v2ConstraintElement()
{
  if(C0) delete [] C0;
}

void
DotType2v2ConstraintElement::getConstants(CoordSet& cs, Eigen::Array<double,9,1>& sconst, Eigen::Array<int,0,1>&, GeomState *gs)
{
  // if rotdesc = 0, or rotdesc = 1/2 and gs = NULL which means that the reference/updated configuration
  // is identical to the undeformed configuration
  if(rotdescr == 0 || gs == NULL) {
    Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&DotType2v2ConstraintElement::C0[0][0]);
    sconst << cs[nn[1]]->x - cs[nn[0]]->x, cs[nn[1]]->y - cs[nn[0]]->y, cs[nn[1]]->z - cs[nn[0]]->z,
              C0(axis,0), C0(axis,1), C0(axis,2),
              C0(axis,0), C0(axis,1), C0(axis,2);
  }
  else { // if rotdesc = 1/2 and gs != NULL
         // gs is interpreted as the reference configuration for updated lagrangian and the current configuration for eulerian
    Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&DotType2v2ConstraintElement::C0[0][0]),
                                                           R1(&(*gs)[nn[0]].R[0][0]),
                                                           R2(&(*gs)[nn[1]].R[0][0]);
    Eigen::Matrix<double,1,3> b0 = C0.row(axis)*R1.transpose(),
                              c0 = C0.row(axis)*R2.transpose();
    sconst << cs[nn[1]]->x - cs[nn[0]]->x, cs[nn[1]]->y - cs[nn[0]]->y, cs[nn[1]]->z - cs[nn[0]]->z,
              b0[0], b0[1], b0[2],
              c0[0], c0[1], c0[2];
  }
}

void
DotType2v2ConstraintElement::setFrame(EFrame *elemframe)
{
  C0 = new double[3][3];
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) 
      C0[i][j] = (*elemframe)[i][j];
}

void 
DotType2v2ConstraintElement::buildFrame(CoordSet& cs)
{
  if(!C0) {
    cerr << " *** ERROR: element frame is not defined for DotType2v2 constraint function element\n";
    exit(-1);
  }

  ConstraintFunctionElement<DotType2v2ConstraintFunction>::buildFrame(cs);
}
#endif
