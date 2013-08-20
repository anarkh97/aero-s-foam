#include <Utils.d/dbg_alloca.h>
#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Sommerfeld.d/QuadPressureBC.h>

QuadPressureBC::QuadPressureBC(int* _nn, PressureBCond* _pbc)
 : PressureElement<Quad4LagrangePolynomialSurfacePressureForceFunction>(4, DofSet::XYZdisp, _nn, _pbc)
{}

void
QuadPressureBC::getConstants(CoordSet& cs, Eigen::Array<double,21,1> &sconst, Eigen::Array<int,2,1> &iconst)
{
  if(!(pbc->conwep && pbc->conwepswitch)) {
    sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z,
              cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
              cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z,
              cs[nn[3]]->x, cs[nn[3]]->y, cs[nn[3]]->z,
              pbc->val, 0, 0, 0, 0, 0, 0, 0, 0;
    iconst << 2, // quadrature rule degree
              0;
  }
  else {
    sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z,
              cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
              cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z,
              cs[nn[3]]->x, cs[nn[3]]->y, cs[nn[3]]->z,
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

#else
#include <Element.d/Sommerfeld.d/QuadPressureBC.h>
#include <Utils.d/dbg_alloca.h>
#include <Corotational.d/GeomState.h>

extern "C" {
  void _FORTRAN(elefbc3dbrkshl2)(int&, int&, double*, double*, double*, double*);
};

QuadPressureBC::QuadPressureBC(int *_nn, PressureBCond *_pbc)
{
  nnode = 4;
  nndof = 3;
  ndime = 3;
  optele = 3;
  nn[0] = _nn[0];
  nn[1] = _nn[1]; 
  nn[2] = _nn[2]; 
  nn[3] = _nn[3]; 
  pbc = _pbc;
  dom = 0;
}

FullSquareMatrix
QuadPressureBC::sommerMatrix(CoordSet &cs, double *d)
{
  FullSquareMatrix sommerM(12,d);
  sommerM.zero();

  return sommerM;
}

void
QuadPressureBC::neumVector(CoordSet &cs, Vector &f, int, GeomState *geomState, double t)
{
  double pressure = pbc->val;
  // Check if Conwep is being used. If so, add the pressure from the blast loading function.
  if (pbc->conwep && pbc->conwepswitch) {
    double* CurrentElementNodePositions = (double*) dbg_alloca(sizeof(double)*3*4);
    int Offset;
    for(int i = 0; i < 4; ++i) {
      Offset = i*3;
      CurrentElementNodePositions[Offset+0] = cs[nn[i]]->x;
      CurrentElementNodePositions[Offset+1] = cs[nn[i]]->y;
      CurrentElementNodePositions[Offset+2] = cs[nn[i]]->z;
    }
    pressure += BlastLoading::ComputeShellPressureLoad(CurrentElementNodePositions, t, *(pbc->conwep));
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

  _FORTRAN(elefbc3dbrkshl2)(opttrc, optele, ecord, edisp, trac, f.data());
}

#include <Driver.d/Domain.h>

int*
QuadPressureBC::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[nndof*nnode];

  for(int i = 0; i < nnode; ++i)
    dsa.number(nn[i], DofSet::XYZdisp, p+i*nndof);

  return p;
}

int
QuadPressureBC::numDofs()
{
  return nndof*nnode;
}

void
QuadPressureBC::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, nnode, DofSet::XYZdisp);
}

void
QuadPressureBC::getNormal(CoordSet &cs, double normal[3])
{
  double x[4], y[4], z[4];

  Node nd1 = cs.getNode(nn[0]);
  Node nd2 = cs.getNode(nn[1]);
  Node nd3 = cs.getNode(nn[2]);
  Node nd4 = cs.getNode(nn[3]);

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

  // Define normal vector as vector product of diagonals

  double nx, ny, nz;

  nx = (y[2]-y[0])*(z[3]-z[1]) - (z[2]-z[0])*(y[3]-y[1]);
  ny = (z[2]-z[0])*(x[3]-x[1]) - (x[2]-x[0])*(z[3]-z[1]);
  nz = (x[2]-x[0])*(y[3]-y[1]) - (y[2]-y[0])*(x[3]-x[1]);

  double l = sqrt(nx*nx+ny*ny+nz*nz);

  normal[0] = nx/l;
  normal[1] = ny/l;
  normal[2] = nz/l;
}

#endif
