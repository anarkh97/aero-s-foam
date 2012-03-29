#include <Element.d/Triangle3.d/ThreeNodeTri3D.h>
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Element.d/State.h>
#include <Hetero.d/InterpPoint.h>
#include <Utils.d/dbg_alloca.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/PhantomCorotator.h>

extern "C" {
void _FORTRAN(elefbc3dtri)(int&, int&, double*, double*, double*, double*);
}

ThreeNodeTri3D::ThreeNodeTri3D(int* nodenums)
{
  for(int i = 0; i < numNodes(); ++i)
    nn[i] = nodenums[i];
}

void
ThreeNodeTri3D::renum(int *table)
{
  for(int i = 0; i < numNodes(); ++i)
    nn[i] = table[nn[i]];
}

int
ThreeNodeTri3D::numNodes()
{
  return 3;
}

int*
ThreeNodeTri3D::nodes(int *p)
{
  if(p == 0) p = new int[numNodes()];

  for(int i = 0; i < numNodes(); ++i) 
    p[i] = nn[i];

  return p;
}

int
ThreeNodeTri3D::numDofs()
{
  return 3*numNodes();
}

int*
ThreeNodeTri3D::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[numDofs()];

  for(int i = 0; i < numNodes(); ++ i)
    dsa.number(nn[i],DofSet::XYZdisp, p + 3*i);

  return p;
}

void
ThreeNodeTri3D::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

FullSquareMatrix
ThreeNodeTri3D::stiffness(CoordSet &cs, double *d, int flg)
{
  FullSquareMatrix result(1, d);
  result.setSize(numDofs());
  result.zero();
  return result;
}

FullSquareMatrix
ThreeNodeTri3D::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  FullSquareMatrix result(1, mel);
  result.setSize(numDofs());
  result.zero();
  return result;
}

double
ThreeNodeTri3D::getMass(CoordSet& cs)
{
  return 0.0;
}

void
ThreeNodeTri3D::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                                Vector& gravityForce, int gravflg, GeomState *geomState)
{
  gravityForce.zero();
}

void
ThreeNodeTri3D::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
                            Vector& elDisp, int strInd, int, double *ndTemps,
                            double ylayer, double zlayer, int avgnum)
{
  stress.zero();
  weight.zero();
}

void
ThreeNodeTri3D::getAllStress(FullM& stress, Vector& weight, CoordSet &cs,
                             Vector& elDisp, int strInd, int, double *ndTemps)
{
  stress.zero();
  weight.zero();
}

void
ThreeNodeTri3D::computeDisp(CoordSet &cs, State &state, const InterpPoint &ip,
                            double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double xyz[3][6];
  state.getDV(nn[0], xyz[0], xyz[0]+3);
  state.getDV(nn[1], xyz[1], xyz[1]+3);
  state.getDV(nn[2], xyz[2], xyz[2]+3);

  for(int j = 0; j < 6; ++j)
    res[j] = (1.0-gp[0]-gp[1]) * xyz[0][j] + gp[0]*xyz[1][j] + gp[1]*xyz[2][j];
}

void
ThreeNodeTri3D::getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, 
                          double *resF, GeomState *gs)
{
  const double *gp = ip.xy;
  for(int i = 0; i < 3; ++i) {
    resF[i]   = (1.0-gp[0]-gp[1]) * flF[i];
    resF[3+i] = gp[0] * flF[i];
    resF[6+i] = gp[1] * flF[i];
  }
}

void
ThreeNodeTri3D::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                     GeomState *geomState, int cflg, double)
{
  int opttrc = 0; // 0 : pressure
                  // 1 : traction
  double* ecord = (double*) dbg_alloca(sizeof(double)*numNodes()*3);
  double* edisp = (double*) dbg_alloca(sizeof(double)*numNodes()*3);
  for(int i = 0; i < numNodes(); ++i) {
    ecord[i*3+0] = cs[nn[i]]->x;
    ecord[i*3+1] = cs[nn[i]]->y;
    ecord[i*3+2] = cs[nn[i]]->z;
    for(int j = 0; j < 3; ++j) edisp[i*3+j] = (geomState) ? (*geomState)[nn[i]].d[j] : 0;
  }
  double trac[3] = { -pressure, 0, 0 };
  double *efbc = (double*) dbg_alloca(sizeof(double)*numNodes()*3);

  int optele = 3;
  _FORTRAN(elefbc3dtri)(opttrc, optele, ecord, edisp, trac, elPressureForce.data());
}

int
ThreeNodeTri3D::getTopNumber()
{
  return 104;
}

Corotator *
ThreeNodeTri3D::getCorotator(CoordSet &, double *, int , int)
{
  return new PhantomCorotator();
}

