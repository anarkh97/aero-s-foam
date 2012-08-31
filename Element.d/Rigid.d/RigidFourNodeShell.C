#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/RigidFourNodeShell.h>
#include <Element.d/Rigid.d/RigidBeam.h>

RigidFourNodeShell::RigidFourNodeShell(int *_nn)
 : SuperElement(true)
{
  nnodes = 4;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  for(int i = 0; i < nSubElems; ++i) {
    int indices[2] = { i+1, 0 };
    subElems[i] = new RigidBeam(indices);
  }
}

//EXPERIMENTAL: equip this element with mass matrix and pressure load vector from BelytschkoTsayShell
#include <Utils.d/dbg_alloca.h>
#include <Corotational.d/GeomState.h>
#include <Element.d/NonLinearity.d/ExpMat.h>

extern "C" {
  void _FORTRAN(elemaslbt)(int&, double*, double*, double*, double*);
  void _FORTRAN(elefbc3dbrkshl2)(int&, int&, double*, double*, double*, double*);
}

void
RigidFourNodeShell::setPressure(double pres, MFTTData *)
{
  pressure = pres;
}

double
RigidFourNodeShell::getPressure()
{ 
  return pressure;
}

FullSquareMatrix
RigidFourNodeShell::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  int nndof = 6, ndime = 3;
  FullSquareMatrix ret(numDofs(), mel);
  ret.zero();

  // Check for element which has no mass
  if(prop && prop->rho != 0 && prop->eh != 0) {
    double* ecord = new double[nnodes*ndime];
    double* edisp = new double[nnodes*nndof];
    for(int i = 0; i < nnodes; ++i) {
      ecord[i*ndime+0] = cs[nn[i]]->x;
      ecord[i*ndime+1] = cs[nn[i]]->y;
      ecord[i*ndime+2] = cs[nn[i]]->z;
      for(int j = 0; j < nndof; ++j) edisp[i*nndof+j] = 0; // constant mass matrix
    }

    double* emasl = (double*) dbg_alloca(sizeof(double)*nnodes*nndof);
    double ematpro[20];
    ematpro[2] = prop->rho; ematpro[19] = prop->eh;
    // get bt shell element lumped mass
    _FORTRAN(elemaslbt)(nndof, ematpro, ecord, edisp, emasl);
       // input : nndof,ematpro,ecord,edisp
       // output : emasl
    delete [] ecord;
    delete [] edisp;

    for(int i = 0; i < nnodes*nndof; ++i) ret[i][i] = emasl[i];
  }

  return ret;
}

void
RigidFourNodeShell::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                         GeomState *geomState, int cflg, double)
{
  int opttrc = 0; // 0 : pressure
                  // 1 : traction
  int optele = 3, ndime = 3;
  double* ecord = (double*) dbg_alloca(sizeof(double)*nnodes*ndime);
  double* edisp = (double*) dbg_alloca(sizeof(double)*nnodes*ndime); // translations only
  for(int i = 0; i < nnodes; ++i) {
    ecord[i*ndime+0] = cs[nn[i]]->x;
    ecord[i*ndime+1] = cs[nn[i]]->y;
    ecord[i*ndime+2] = cs[nn[i]]->z;
    for(int j = 0; j < 3; ++j) edisp[i*ndime+j] = (*geomState)[nn[i]].d[j];
  }
  double trac[3] = { -pressure, 0, 0 };
  double *efbc = (double*) dbg_alloca(sizeof(double)*nnodes*ndime); // translations only

  _FORTRAN(elefbc3dbrkshl2)(opttrc, optele, ecord, edisp, trac, efbc);

  for(int i = 0; i < nnodes; ++i)
    for(int j = 0; j < ndime; ++j)
      elPressureForce[6*i+j] = efbc[3*i+j];

  for(int i=nnodes*ndime; i<numDofs(); ++i) elPressureForce[i] = 0; // lagrange multiplier dofs, if any
}

#include <Element.d/State.h>
#include <Hetero.d/InterpPoint.h>

void
RigidFourNodeShell::computeDisp(CoordSet& cs, State& state, const InterpPoint& ip,
                                double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double xyz[4][6];
  state.getDV(nn[0], xyz[0], xyz[0]+3);
  state.getDV(nn[1], xyz[1], xyz[1]+3);
  state.getDV(nn[2], xyz[2], xyz[2]+3);
  state.getDV(nn[3], xyz[3], xyz[3]+3);

  int j;
  for(j=0; j<6; ++j)
    res[j] = (1-gp[0])*(1-gp[1])* xyz[0][j] +
             gp[0]*(1-gp[1])* xyz[1][j] +
             (1-gp[0])*gp[1]* xyz[3][j] +
             gp[0]*gp[1]*xyz[2][j];
}

void
RigidFourNodeShell::getFlLoad(CoordSet& cs, const InterpPoint& ip, double *flF,
                              double *resF, GeomState *gs)
{
  // PJSA 9/10/2010 reversed resF[12+i] and resF[18+i] to match 
  // FlExchanger::getQuadFlLoad in Xfem/Hetero.d/FlExchange.C
  const double *gp = ip.xy;
  for(int i = 0; i < 3; ++i) {
    resF[i]    = (1-gp[0])*(1-gp[1])* flF[i];
    resF[6+i]  = gp[0]*(1-gp[1])* flF[i];
    resF[12+i] = gp[0]*gp[1]* flF[i];
    resF[18+i] = (1-gp[0])*gp[1]* flF[i];
    resF[i+3]  = resF[i+9] = resF[i+15] = resF[i+21] = 0.0;
  }
}
#endif
