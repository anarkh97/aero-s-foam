#include <Utils.d/dbg_alloca.h>
#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Sommerfeld.d/TrianglePressureBC.h>

TrianglePressureBC::TrianglePressureBC(int* _nn, double _pressure)
 : PressureElement<Tri3LagrangePolynomialSurfacePressureForceFunction>(3, DofSet::XYZdisp, _nn),
   pressure(_pressure)
{
}

void
TrianglePressureBC::getConstants(CoordSet& cs, Eigen::Array<double,17,1> &sconst, Eigen::Array<int,2,1> &iconst)
{
  if(!conwep) {
    sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z,
              cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
              cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z,
              pressure, 0, 0, 0, 0, 0, 0, 0;
    iconst << 1, // quadrature rule degree
              0;
  }
  else {
    sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z,
              cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
              cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z,
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

#else
#include <Element.d/Sommerfeld.d/TrianglePressureBC.h>
#include <Utils.d/dbg_alloca.h>
#include <Corotational.d/GeomState.h>

extern "C" {
  void _FORTRAN(elefbc3dtri)(int&, int&, double*, double*, double*, double*);
};


TrianglePressureBC::TrianglePressureBC(int *_nn, double _pressure)
{
  nnode = 3;
  nndof = 3;
  ndime = 3;
  optele = 3;
  nn[0] = _nn[0];
  nn[1] = _nn[1]; 
  nn[2] = _nn[2]; 
  pressure = _pressure;
  conwep = 0;
  dom = 0;
}

FullSquareMatrix
TrianglePressureBC::sommerMatrix(CoordSet &cs, double *d)
{
  FullSquareMatrix sommerM(9,d);
  sommerM.zero();

  return sommerM;
}

void
TrianglePressureBC::neumVector(CoordSet &cs, Vector &f, int, GeomState *geomState, double t)
{
  // Check if Conwep is being used. If so, use the pressure from the blast loading function.
  if (conwep) {
    double* CurrentElementNodePositions = (double*) dbg_alloca(sizeof(double)*3*4);
    int Offset;
    for(int i = 0; i < 4; ++i) {
      Offset = i*3;
      if (i==3) {
        CurrentElementNodePositions[Offset+0] = cs[nn[2]]->x;
        CurrentElementNodePositions[Offset+1] = cs[nn[2]]->y;
        CurrentElementNodePositions[Offset+2] = cs[nn[2]]->z;
      }
      else {
        CurrentElementNodePositions[Offset+0] = cs[nn[i]]->x;
        CurrentElementNodePositions[Offset+1] = cs[nn[i]]->y;
        CurrentElementNodePositions[Offset+2] = cs[nn[i]]->z;
      }
    }
    pressure = BlastLoading::ComputeShellPressureLoad(CurrentElementNodePositions, t, *conwep);
  }
  int opttrc = 0; // 0 : pressure
                  // 1 : traction
  double* ecord = (double*) dbg_alloca(sizeof(double)*nnode*ndime);
  double* edisp = (double*) dbg_alloca(sizeof(double)*nnode*ndime); // translations only
  for(int i = 0; i < nnode; ++i) {
    ecord[i*ndime+0] = cs[nn[i]]->x;
    ecord[i*ndime+1] = cs[nn[i]]->y;
    ecord[i*ndime+2] = cs[nn[i]]->z;
    edisp[i*ndime+0] = (geomState != NULL) ? (*geomState)[nn[i]].x - cs[nn[i]]->x : 0;
    edisp[i*ndime+1] = (geomState != NULL) ? (*geomState)[nn[i]].y - cs[nn[i]]->y : 0;
    edisp[i*ndime+2] = (geomState != NULL) ? (*geomState)[nn[i]].z - cs[nn[i]]->z : 0;
  }
  double trac[3] = { -pressure, 0, 0 };

  _FORTRAN(elefbc3dtri)(opttrc, optele, ecord, edisp, trac, f.data());
}

#include <Driver.d/Domain.h>

int*
TrianglePressureBC::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[nndof*nnode];

  for(int i = 0; i < nnode; ++i)
    dsa.number(nn[i], DofSet::XYZdisp, p+i*nndof);

  return p;
}

int
TrianglePressureBC::numDofs()
{
  return nndof*nnode;
}

void
TrianglePressureBC::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, nnode,  DofSet::XYZdisp);
}

void
TrianglePressureBC::getNormal(CoordSet &cs, double normal[3])
{
  Node nd1 = cs.getNode(nn[0]);
  Node nd2 = cs.getNode(nn[1]);
  Node nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];
  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  double u1, u2, u3, v1, v2, v3, w1, w2, w3;
  u1 = x[1]-x[0];
  u2 = y[1]-y[0];
  u3 = z[1]-z[0];

  v1 = x[2]-x[0];
  v2 = y[2]-y[0];
  v3 = z[2]-z[0];

  w1 = u2*v3-u3*v2;
  w2 = u3*v1-u1*v3;
  w3 = u1*v2-u2*v1;

  double l = sqrt(w1*w1+w2*w2+w3*w3);

  if (l == 0) l = 1.0;

  normal[0] = w1/l;
  normal[1] = w2/l;
  normal[2] = w3/l;
}

#endif
