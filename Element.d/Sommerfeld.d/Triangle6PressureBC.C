#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Sommerfeld.d/Triangle6PressureBC.h>

Triangle6PressureBC::Triangle6PressureBC(int* _nn, double _pressure)
 : PressureElement<Tri6LagrangePolynomialSurfacePressureForceFunction>(6, DofSet::XYZdisp, _nn),
   pressure(_pressure)
{
}

void
Triangle6PressureBC::getConstants(CoordSet& cs, Eigen::Array<double,26,1> &sconst, Eigen::Array<int,2,1> &iconst)
{
  if(!conwep) {
    sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z,
              cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
              cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z,
              cs[nn[3]]->x, cs[nn[3]]->y, cs[nn[3]]->z,
              cs[nn[4]]->x, cs[nn[4]]->y, cs[nn[4]]->z,
              cs[nn[5]]->x, cs[nn[5]]->y, cs[nn[5]]->z,
              pressure, 0, 0, 0, 0, 0, 0, 0;
    iconst << 1, // quadrature rule degree
              0;
  }
  else {
    sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z,
              cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
              cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z,
              cs[nn[3]]->x, cs[nn[3]]->y, cs[nn[3]]->z,
              cs[nn[4]]->x, cs[nn[4]]->y, cs[nn[4]]->z,
              cs[nn[5]]->x, cs[nn[5]]->y, cs[nn[5]]->z,
              conwep->ExplosivePosition[0],
              conwep->ExplosivePosition[1],
              conwep->ExplosivePosition[2],
              conwep->ExplosiveDetonationTime,
              conwep->ExplosiveWeight,
              conwep->ScaleLength,
              conwep->ScaleTime,
              conwep->ScaleMass;
     iconst << 1, // quadrature rule degree
               (conwep->BlastType == BlastLoading::BlastData::SurfaceBurst ? 1 : 2);
  }
}
#endif
