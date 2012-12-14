#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Sommerfeld.d/Quad8PressureBC.h>

Quad8PressureBC::Quad8PressureBC(int* _nn, double _pressure)
 : PressureElement<Quad8LagrangePolynomialSurfacePressureForceFunction>(8, DofSet::XYZdisp, _nn),
   pressure(_pressure)
{
}

void
Quad8PressureBC::getConstants(CoordSet& cs, Eigen::Array<double,25,1> &sconst, Eigen::Array<int,1,1> &iconst)
{
  sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z,
            cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
            cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z,
            cs[nn[3]]->x, cs[nn[3]]->y, cs[nn[3]]->z,
            cs[nn[4]]->x, cs[nn[4]]->y, cs[nn[4]]->z,
            cs[nn[5]]->x, cs[nn[5]]->y, cs[nn[5]]->z,
            cs[nn[6]]->x, cs[nn[6]]->y, cs[nn[6]]->z,
            cs[nn[7]]->x, cs[nn[7]]->y, cs[nn[7]]->z,
            pressure;
  iconst << 2; // quadrature rule degree
}
#endif
