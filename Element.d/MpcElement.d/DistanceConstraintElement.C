#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/DistanceConstraintElement.h>

DistanceConstraintElement::DistanceConstraintElement(int* _nn, double _f0, int _type)
 : ConstraintFunctionElement<DistanceConstraintFunction>(2, DofSet::XYZdisp, _nn, _type)
{
  f0 = _f0;
}

void
DistanceConstraintElement::getConstants(CoordSet& cs, Eigen::Array<double,4,1>& sconst, Eigen::Array<int,0,1>&, GeomState*)
{
  sconst << cs[nn[1]]->x - cs[nn[0]]->x, cs[nn[1]]->y - cs[nn[0]]->y, cs[nn[1]]->z - cs[nn[0]]->z,
            f0;
}
#endif
