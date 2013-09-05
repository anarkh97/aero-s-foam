#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Sommerfeld.d/Quad9PressureBC.h>

Quad9PressureBC::Quad9PressureBC(int* _nn, PressureBCond* _pbc)
 : PressureElement<Quad9LagrangePolynomialSurfacePressureForceFunction>(9, DofSet::XYZdisp, _nn, _pbc)
{}

void
Quad9PressureBC::getConstants(CoordSet& cs, Eigen::Array<double,36,1> &sconst, Eigen::Array<int,2,1> &iconst)
{
  if(!(pbc->conwep && pbc->conwepswitch)) {
    sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z,
              cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
              cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z,
              cs[nn[3]]->x, cs[nn[3]]->y, cs[nn[3]]->z,
              cs[nn[4]]->x, cs[nn[4]]->y, cs[nn[4]]->z,
              cs[nn[5]]->x, cs[nn[5]]->y, cs[nn[5]]->z,
              cs[nn[6]]->x, cs[nn[6]]->y, cs[nn[6]]->z,
              cs[nn[7]]->x, cs[nn[7]]->y, cs[nn[7]]->z,
              cs[nn[8]]->x, cs[nn[8]]->y, cs[nn[8]]->z,
              pbc->val, 0, 0, 0, 0, 0, 0, 0, 0;
    iconst << 2, // quadrature rule degree
              0;
  }
  else {
    sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z,
              cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
              cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z,
              cs[nn[3]]->x, cs[nn[3]]->y, cs[nn[3]]->z,
              cs[nn[4]]->x, cs[nn[4]]->y, cs[nn[4]]->z,
              cs[nn[5]]->x, cs[nn[5]]->y, cs[nn[5]]->z,
              cs[nn[6]]->x, cs[nn[6]]->y, cs[nn[6]]->z,
              cs[nn[7]]->x, cs[nn[7]]->y, cs[nn[7]]->z,
              cs[nn[8]]->x, cs[nn[8]]->y, cs[nn[8]]->z,
              pbc->val,
              pbc->conwep->ExplosivePosition[0],
              pbc->conwep->ExplosivePosition[1],
              pbc->conwep->ExplosivePosition[2],
              pbc->conwep->ExplosiveDetonationTime,
              pbc->conwep->ExplosiveWeight,
              pbc->conwep->ScaleLength,
              pbc->conwep->ScaleTime,
              pbc->conwep->ScaleMass;
     iconst << 3, // quadrature rule degree
               (pbc->conwep->BlastType == BlastLoading::BlastData::SurfaceBurst ? 1 : 2);
  }
}
#endif
