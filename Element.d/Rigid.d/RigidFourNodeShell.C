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
RigidFourNodeShell::setPressure(double pres)
{
  pressure = pres;
}

double
RigidFourNodeShell::getPressure()
{ 
  return pressure;
}

void
RigidFourNodeShell::setMaterial(NLMaterial *m)
{
  expmat = dynamic_cast<ExpMat *>(m);
}

FullSquareMatrix
RigidFourNodeShell::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  int nndof = 6, ndime = 3;
  FullSquareMatrix ret(nnodes*nndof, mel);
  ret.zero();

  // Check for phantom element, which has no mass
  if(prop) {
    double* ecord = new double[nnodes*ndime];
    double* edisp = new double[nnodes*nndof];
    for(int i = 0; i < nnodes; ++i) {
      ecord[i*ndime+0] = cs[nn[i]]->x;
      ecord[i*ndime+1] = cs[nn[i]]->y;
      ecord[i*ndime+2] = cs[nn[i]]->z;
      for(int j = 0; j < nndof; ++j) edisp[i*nndof+j] = 0; // constant mass matrix
    }

    double* emasl = (double*) dbg_alloca(sizeof(double)*nnodes*nndof);
    // get bt shell element lumped mass
    _FORTRAN(elemaslbt)(nndof, expmat->ematpro, ecord, edisp, emasl);
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
                                         GeomState *geomState, int cflg)
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
}

