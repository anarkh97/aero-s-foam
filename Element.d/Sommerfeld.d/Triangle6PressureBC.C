#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Sommerfeld.d/Triangle6PressureBC.h>

Triangle6PressureBC::Triangle6PressureBC(int* _nn, PressureBCond* _pbc)
 : PressureElement<Tri6LagrangePolynomialSurfacePressureForceFunction>(6, DofSet::XYZdisp, _nn, _pbc)
{}

void
Triangle6PressureBC::getConstants(CoordSet& cs, Eigen::Array<double,27,1> &sconst, Eigen::Array<int,2,1> &iconst)
{
  if(!(pbc->conwep && pbc->conwepswitch)) {
    sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z,
              cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
              cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z,
              cs[nn[3]]->x, cs[nn[3]]->y, cs[nn[3]]->z,
              cs[nn[4]]->x, cs[nn[4]]->y, cs[nn[4]]->z,
              cs[nn[5]]->x, cs[nn[5]]->y, cs[nn[5]]->z,
              pbc->val, 0, 0, 0, 0, 0, 0, 0, 0;
    iconst << 3, // number of gauss points
              0;
  }
  else {
    sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z,
              cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
              cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z,
              cs[nn[3]]->x, cs[nn[3]]->y, cs[nn[3]]->z,
              cs[nn[4]]->x, cs[nn[4]]->y, cs[nn[4]]->z,
              cs[nn[5]]->x, cs[nn[5]]->y, cs[nn[5]]->z,
              pbc->val,
              pbc->conwep->ExplosivePosition[0],
              pbc->conwep->ExplosivePosition[1],
              pbc->conwep->ExplosivePosition[2],
              pbc->conwep->ExplosiveDetonationTime,
              pbc->conwep->ExplosiveWeight,
              pbc->conwep->ScaleLength,
              pbc->conwep->ScaleTime,
              pbc->conwep->ScaleMass;
     iconst << 6, // number of gauss points
               (pbc->conwep->BlastType == BlastLoading::BlastData::SurfaceBurst ? 1 : 2);
  }
}
#endif
