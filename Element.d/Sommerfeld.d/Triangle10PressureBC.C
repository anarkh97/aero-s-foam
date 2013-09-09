#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Sommerfeld.d/Triangle10PressureBC.h>

Triangle10PressureBC::Triangle10PressureBC(int* _nn, PressureBCond* _pbc)
 : PressureElement<Tri10LagrangePolynomialSurfacePressureForceFunction>(10, DofSet::XYZdisp, _nn, _pbc)
{}

void
Triangle10PressureBC::getConstants(CoordSet& cs, Eigen::Array<double,39,1> &sconst, Eigen::Array<int,2,1> &iconst)
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
              cs[nn[9]]->x, cs[nn[9]]->y, cs[nn[9]]->z,
              pbc->val, 0, 0, 0, 0, 0, 0, 0, 0;
    iconst << 1, // XXX quadrature rule degree
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
              cs[nn[9]]->x, cs[nn[9]]->y, cs[nn[9]]->z,
              pbc->val,
              pbc->conwep->ExplosivePosition[0],
              pbc->conwep->ExplosivePosition[1],
              pbc->conwep->ExplosivePosition[2],
              pbc->conwep->ExplosiveDetonationTime,
              pbc->conwep->ExplosiveWeight,
              pbc->conwep->ScaleLength,
              pbc->conwep->ScaleTime,
              pbc->conwep->ScaleMass;
     iconst << 1, // XXX quadrature rule degree
               (pbc->conwep->BlastType == BlastLoading::BlastData::SurfaceBurst ? 1 : 2);
  }
}
#endif
