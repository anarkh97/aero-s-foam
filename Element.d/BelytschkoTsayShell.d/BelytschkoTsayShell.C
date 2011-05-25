#include <Element.d/BelytschkoTsayShell.d/BelytschkoTsayShell.h>
#include <Utils.d/dbg_alloca.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Element.d/NonLinearity.d/ExpMat.h>
#include <Element.d/State.h>
#include <Hetero.d/InterpPoint.h>
#include <Material.d/IsotropicLinearElasticJ2PlasticPlaneStressMaterial.h>
#include <Material.d/KorkolisKyriakidesPlaneStressMaterial.h>
#include <Material.d/KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding.h>

#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <vector>

using namespace Eigen;
using std::vector;
#endif

extern "C" {
  void _FORTRAN(getgqsize)(int&, int&, int*, int*, int*);
  void _FORTRAN(getgq1d)(int&, double*, double*);
  void _FORTRAN(elemaslbt)(int&, double*, double*, double*, double*);

  void _FORTRAN(getlocbvecbt)(double*, double*);
  void _FORTRAN(getbmat1pt)(double*, double&, double*);
  void _FORTRAN(getgamma4nod)(double*, double&, double*, double&);
  void _FORTRAN(getbcmat1pt)(double*, double&, double*, double&, double*);
  void _FORTRAN(getbsmat1pt)(double*, double*);
  void _FORTRAN(updstrnbt)(int*, double&, double*, double&, double*,
                           double*, double*, double*,
                           double*, double*, double*);
  void _FORTRAN(updstrsbt)(int&, int&, double&, double*, double&, double*,
                           double*, double&, double&, double*,
                           double*);
  void _FORTRAN(gqfintbt)(int*, double&, double*, double&, double&,
                          double&, double*, double*, double*,
                          double*, double*);
  void _FORTRAN(updhgcstrsbt)(double*, double&, double*, double*, double&,
                              double*, double*, double&, double*);
  void _FORTRAN(gqfhgcbt)(double*, double&, double*, double*);
  void _FORTRAN(elefbc3dbrkshl2opt)(double&, double*, double&, double*);
}

double BelytschkoTsayShell::t1 = 0;
double BelytschkoTsayShell::t2 = 0;
double BelytschkoTsayShell::t3 = 0;
double BelytschkoTsayShell::t4 = 0;
double BelytschkoTsayShell::t5 = 0;
double BelytschkoTsayShell::t6 = 0;
double BelytschkoTsayShell::t7 = 0;

BelytschkoTsayShell::BelytschkoTsayShell(int* nodenums)
{
  // TODO most of these should come from StructProp via MATERIAL
  optdmg = 0; // no damage
  opthgc = 1; // perturbation type hourglass control
  opttrc = -1; // no pressure or traction
  optdmp = 0; // no damping
  optcor[0] = 1; // was 1 // warping correction on
  optcor[1] = 0; // was 0 // shear correction off
  nndof  = 6; // number of dofs per node
  ndime  = 3;
  nnode  = 4;
  ngqpt[0] = 1; 
  ngqpt[1] = 1; 
  ngqpt[2] = 3;
  ngqpt4 = 3;
  // hourglass control parameters
  prmhgc[0] = 2.5e-3;
  prmhgc[1] = 2.5e-3;
  prmhgc[2] = 2.5e-3;
  for(int i = 3; i < 10; ++i) prmhgc[i] = 0;
  // damping parameters
  prmdmp[0] = 4e-1;
  prmdmp[1] = 1e6;
  for(int i = 2; i < 10; ++i) prmdmp[i] = 0;

  for(int i = 0; i < nnode; ++i)
    nn[i] = nodenums[i];

  // ---------------------------------------------------------------
  // set gq size
  // -----------
  int optmhd = 0, optele = 3;
  _FORTRAN(getgqsize)(optmhd, optele, ngqpt, mgaus, mgqpt);
     // input : optmhd,optele,ngqpt
     // output : mgaus,mgqpt

  // ---------------------------------------------------------------
  // get 1d gq rule: through thickness integration
  // --------------
  gqpoin3 = new double[mgaus[2]];
  gqweigt3 = new double[mgaus[2]];
  _FORTRAN(getgq1d)(mgaus[2], gqpoin3, gqweigt3);
     // input : mgaus[2]
     // output : gqpoin3,gqweigt3

  // ---------------------------------------------------------------
  // allocate and initialize history variables
  // ---------------------
  // TODO depends on optctv, see inihistvar3d
  evar1 = new double[5*mgqpt[0]];
  evar2 = new double[5*mgqpt[0]];
  evoit1 = new double[6*mgqpt[0]];
  evoit2 = new double[6*mgqpt[0]];
  evoit3 = new double[6*mgqpt[0]];
  for(int i = 0; i < 5*mgqpt[0]; ++i) evar1[i] = evar2[i] = 0;
  for(int i = 0; i < 6*mgqpt[0]; ++i) evoit1[i] = evoit2[i] = evoit3[i] = 0;

  expmat = 0;
  mftt = 0;
}

BelytschkoTsayShell::~BelytschkoTsayShell()
{
  delete [] gqpoin3;
  delete [] gqweigt3;
  delete [] evar1;
  delete [] evar2;
  delete [] evoit1;
  delete [] evoit2;
  delete [] evoit3;
}

void
BelytschkoTsayShell::setMaterial(NLMaterial *m)
{
  expmat = dynamic_cast<ExpMat *>(m);
  if(expmat->optctv != 1) {
    double E = expmat->ematpro[0], nu = expmat->ematpro[1];
    double lambda = E*nu/((1+nu)*(1-2*nu)), mu = E/(2*(1+nu));
    mat = new ElastoPlasticPlaneStressMaterial * [mgaus[2]];
    for(int i=0; i<mgaus[2]; ++i) {
      switch(expmat->optctv) {
      case 5 :
        mat[i] = new IsotropicLinearElasticJ2PlasticPlaneStressMaterial(lambda, mu, expmat->ematpro[3], expmat->ematpro[4], expmat->ematpro[5]);
        break;
      case 6 :
        mat[i] = new KorkolisKyriakidesPlaneStressMaterial(lambda, mu, expmat->ematpro[3], expmat->ematpro[4], expmat->ematpro[5]);
        break;
      case 7 :
        mat[i] = new KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding(lambda, mu);
      }
    }
  }
}

void
BelytschkoTsayShell::setPressure(double _pressure, MFTTData *_mftt)
{
  pressure = _pressure;
  mftt = _mftt;
  opttrc = 0;
}

double
BelytschkoTsayShell::getPressure()
{
  // the return value of this function is used to determine whether or not
  // computePressureForce should be called. Since the pressure for this element
  // in computed along with the internal force inside elefintbt1, it is
  // not necessary to call that function, so we return 0 here
  return double(0);
}

Element *
BelytschkoTsayShell::clone()
{
  return new BelytschkoTsayShell(*this);
}

void
BelytschkoTsayShell::renum(int *table)
{
  for(int i = 0; i < nnode; ++i)
    nn[i] = table[nn[i]];
}

void
BelytschkoTsayShell::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
                                 Vector& elDisp, int strInd, int surface,
                                 double *ndTemps, double ylayer, double zlayer, int avgnum)
{ 
#ifdef USE_EIGEN3
  // voight rule in xed3d code: [xx,yy,zz,yz,xz,xy]
  weight = 1.0;
  int j;
  switch(surface) {
    case 1 : j = 0; break;              // upper
    case 2 : j = (mgaus[2]-1)/2; break; // median (mgaus[2] is the number of gauss points through the thickness)
    case 3 : j = mgaus[2]-1; break;     // lower
  }
  int k[6] = { 0, 1, 2, 5, 3, 4 };
  for(int i = 0; i < nnode; ++i) {
    switch(strInd) {
      case 0 : case 1 : case 2 : case 3 : case 4 : case 5 : // sxx, syy, szz, sxy, syz, sxz
        stress[i] = evoit2[6*j+k[strInd]];
        break;
      case 6 : { // effective stress
        Matrix3d M; M << evoit2[6*j+k[0]], evoit2[6*j+k[3]], evoit2[6*j+k[5]],
                         evoit2[6*j+k[3]], evoit2[6*j+k[1]], evoit2[6*j+k[4]],
                         evoit2[6*j+k[5]], evoit2[6*j+k[4]], evoit2[6*j+k[2]];
        // compute the deviatoric stress/strain tensor and it's second invariant
        Matrix3d dev = M - (M.trace()/3)*Matrix3d::Identity();
        double J2 = 0.5*(dev*dev).trace();
        stress[i] = sqrt(3*J2);
      } break;
      case 7 : case 8 : case 9 : case 10: case 11: case 12: // exx, eyy, ezz, exy, eyz, exz
        stress[i] = evoit1[6*j+k[strInd-7]];
        break;
      case 13 : { // effective strain
        Matrix3d M; M << evoit3[6*j+k[0]], evoit3[6*j+k[3]], evoit3[6*j+k[5]],
                         evoit3[6*j+k[3]], evoit3[6*j+k[1]], evoit3[6*j+k[4]],
                         evoit3[6*j+k[5]], evoit3[6*j+k[4]], evoit3[6*j+k[2]];
        // compute the deviatoric stress/strain tensor and it's second invariant
        Matrix3d dev = M - (M.trace()/3)*Matrix3d::Identity();
        double J2 = 0.5*(dev*dev).trace();
        stress[i] = sqrt(3*J2);

      } break;
      case 17 : // damage for hypoelas
        if(expmat->optctv == 1)
          stress[i] = evar1[5*j+1];
        break;
      case 18 : // effective plastic strain for elasto plastic materials
        if(expmat->optctv != 1)
          stress[i] = mat[j]->GetMaterialEquivalentPlasticStrain();
        break;
    }
  }
#else
  cerr << "USE_EIGEN3 is not defined here in BelytschkoTsayShell::getVonMises\n";
  exit(-1);
#endif
}

void
BelytschkoTsayShell::getAllStress(FullM& stress, Vector& weight, CoordSet &cs,
                                  Vector& elDisp, int strInd, int surface,
                                  double *ndTemps)
{
  cerr << "BelytschkoTsayShell::getAllStress not implemented\n";
}

double
BelytschkoTsayShell::getMass(CoordSet& cs)
{
  if (prop == NULL) return 0.0;

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);

  Vector r1(3), r2(3), r3(3), r4(3);

  r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = 0.0;
  r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = 0.0;
  r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = 0.0;
  r4[0] = nd4.x; r4[1] = nd4.y; r4[2] = 0.0;

  Vector v1(3), v2(3), v3(3), v4(3), v5(3);

  v1 = r2 - r1;
  v2 = r3 - r1;
  v3 = r4 - r1;

  v4 = v1.cross(v2);
  v5 = v2.cross(v3);

  double area = 0.5*(v4.magnitude() + v5.magnitude());
  double rho = expmat->ematpro[2];
  double eh = expmat->ematpro[19];
  double mass = area*rho*eh;

  return mass;
}

void
BelytschkoTsayShell::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                                     Vector& gravityForce, int gravflg, GeomState *geomState)
{
  cerr << "BelytschkoTsayShell::getGravityForce not implemented\n";
  gravityForce.zero();
}

FullSquareMatrix
BelytschkoTsayShell::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  FullSquareMatrix ret(nnode*nndof, mel);
  ret.zero();

  // Check for phantom element, which has no mass
  if(prop) {
    double* ecord = new double[nnode*ndime];
    double* edisp = new double[nnode*nndof];
    for(int i = 0; i < nnode; ++i) {
      ecord[i*ndime+0] = cs[nn[i]]->x;
      ecord[i*ndime+1] = cs[nn[i]]->y;
      ecord[i*ndime+2] = cs[nn[i]]->z;
      for(int j = 0; j < nndof; ++j) edisp[i*nndof+j] = 0; // constant mass matrix
    }

    double* emasl = (double*) dbg_alloca(sizeof(double)*nnode*nndof);
    // get bt shell element lumped mass
    _FORTRAN(elemaslbt)(nndof, expmat->ematpro, ecord, edisp, emasl);
       // input : nndof,ematpro,ecord,edisp
       // output : emasl
    delete [] ecord;
    delete [] edisp;
  
    for(int i = 0; i < nnode*nndof; ++i) ret[i][i] = emasl[i];
  }

  return ret;
}

FullSquareMatrix
BelytschkoTsayShell::stiffness(CoordSet &cs, double *d, int flg)
{
  if(prop) cerr << "BelytschkoTsayShell::stiffness not implemented. This element only works for explicit dynamics\n";
  FullSquareMatrix ret(nnode*nndof, d);
  ret.zero();
  return ret;
}

int
BelytschkoTsayShell::numNodes()
{
  return nnode;
}

int*
BelytschkoTsayShell::nodes(int *p)
{
   if(p == 0) p = new int[nnode];

   for(int i = 0; i < nnode; ++i)
     p[i] = nn[i];

  return p;
}

int
BelytschkoTsayShell::numDofs()
{
  return nndof*nnode;
}

int*
BelytschkoTsayShell::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[nndof*nnode];

  for(int i = 0; i < nnode; ++i) 
    dsa.number(nn[i], DofSet::XYZdisp | DofSet::XYZrot, p+i*nndof);

  return p;
}

void
BelytschkoTsayShell::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, nnode,  DofSet::XYZdisp | DofSet::XYZrot);
}

Corotator *
BelytschkoTsayShell::getCorotator(CoordSet &, double *, int, int)
{
  return this;
}

void
BelytschkoTsayShell::getStiffAndForce(GeomState& geomState, CoordSet& cs, FullSquareMatrix& k, double* efint, double delt, double time)
{
  //=======================================================================
  //  compute internal force vector including hourglass control, pressure
  //  and damping for Belytschko Tsay shell element
  // 
  //  output:
  //  ------
  //  efint(24) : internal force
  //
  // ======================================================================

  if(prop) { // check for phantom
    // ---------------------------------------------------------------
    // get element nodal displacement and velocity
    // -------------------------------------------
    double* ecord = (double*) dbg_alloca(sizeof(double)*nnode*ndime);
    double* edisp = (double*) dbg_alloca(sizeof(double)*nnode*nndof); // d^{n+1}
    double* evelo = (double*) dbg_alloca(sizeof(double)*nnode*nndof); // v^{n+0.5}
    int iloc,jloc;
    for(int i = 0; i < nnode; ++i) {
      iloc = i*ndime;
      ecord[iloc+0] = cs[nn[i]]->x;
      ecord[iloc+1] = cs[nn[i]]->y;
      ecord[iloc+2] = cs[nn[i]]->z;
      for(int j = 0; j < nndof; ++j) {
        jloc = i*nndof+j;
        edisp[jloc] = geomState[nn[i]].d[j];
        evelo[jloc] = geomState[nn[i]].v[j] + delt*0.5*geomState[nn[i]].a[j];
      }
    }
    double trac[3] = { 0, 0, pressure };
    double tmftval = (mftt) ? mftt->getVal(std::max(time,0.0)) : 1.0;

    // ---------------------------------------------------------------
    // internal force, hourglass control and pressure
    // ------------------
    Elefintbt1(delt, ecord, edisp, evelo, trac, tmftval, efint);
     // input : optctv,optdmb,opttrc,prmhgc,delt,ematpro,nndof,mgaus(3),mgqpt[0],gqpoin3,gqweigt3,ecord,edisp,evelo,trac,tmftval
     // inoutput : evar1,evoit2,evoit3
     // output : evar2,efint
 
    // ---------------------------------------------------------------
    // damping force
    // -----------------------
    if(optdmp) {
      double dmpfctr = prmdmp[0]; // damping factor: 0.0-1.0
      double freq = prmdmp[1]; // frequency: 1.0d0 / cycle(time)
      double cnst = 2.0 * freq * dmpfctr;
      double* emasl = (double*) dbg_alloca(sizeof(double)*nnode*nndof);
      _FORTRAN(elemaslbt)(nndof, expmat->ematpro, ecord, edisp, emasl);
      for(int i = 0; i < nnode*nndof; ++i) efint[i] += cnst*emasl[i]*evelo[i];
    }
  }
}

void
BelytschkoTsayShell::computeDisp(CoordSet& cs, State& state, const InterpPoint& ip,
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
BelytschkoTsayShell::getFlLoad(CoordSet& cs, const InterpPoint& ip, double *flF,
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

int
BelytschkoTsayShell::getTopNumber()
{
  return 188;
}

void
BelytschkoTsayShell::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                          GeomState *geomState, int cflg)
{
/* now the pressure force is added in the same routine as the internal force
  int opttrc = 0; // 0 : pressure
                  // 1 : traction
  double* ecord = (double*) dbg_alloca(sizeof(double)*nnode*ndime);
  double* edisp = (double*) dbg_alloca(sizeof(double)*nnode*ndime); // translations only
  for(int i = 0; i < nnode; ++i) {
    ecord[i*ndime+0] = cs[nn[i]]->x;
    ecord[i*ndime+1] = cs[nn[i]]->y;
    ecord[i*ndime+2] = cs[nn[i]]->z;
    for(int j = 0; j < 3; ++j) edisp[i*ndime+j] = (*geomState)[nn[i]].d[j];
  }
  double trac[3] = { -pressure, 0, 0 };
  double *efbc = (double*) dbg_alloca(sizeof(double)*nnode*ndime); // translations only

  _FORTRAN(elefbc3dbrkshl2)(opttrc, optele, ecord, edisp, trac, efbc);

  for(int i = 0; i < nnode; ++i)
    for(int j = 0; j < ndime; ++j)
      elPressureForce[6*i+j] = efbc[3*i+j];
*/
}

void
BelytschkoTsayShell::getThermalForce(CoordSet& cs, Vector& ndTemps,
                                     Vector &elThermalForce, int glflag, 
                                     GeomState *geomState)
{
  cerr << "BelytschkoTsayShell::getThermalForce not implemented\n";
  elThermalForce.zero();
}

// New include files for Restart file
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

void
BelytschkoTsayShell::writeHistory(int fn)
{
  // ---------------------------------------------------------------
  // write history variables to file for restart
  // ---------------------
  int writeSize;

  writeSize = write(fn, evar1, 5*mgqpt[0]*sizeof(double));
  if(writeSize != 5*mgqpt[0]*sizeof(double))
    fprintf(stderr," *** ERROR: Inconsistent restart file 5.1\n");

  writeSize = write(fn, evar2, 5*mgqpt[0]*sizeof(double));
  if(writeSize != 5*mgqpt[0]*sizeof(double))
    fprintf(stderr," *** ERROR: Inconsistent restart file 5.2\n");

  writeSize = write(fn, evoit1, 6*mgqpt[0]*sizeof(double));
  if(writeSize != 6*mgqpt[0]*sizeof(double))
    fprintf(stderr," *** ERROR: Inconsistent restart file 5.3\n");

  writeSize = write(fn, evoit2, 6*mgqpt[0]*sizeof(double));
  if(writeSize != 6*mgqpt[0]*sizeof(double))
    fprintf(stderr," *** ERROR: Inconsistent restart file 5.4\n");

  writeSize = write(fn, evoit3, 6*mgqpt[0]*sizeof(double));
  if(writeSize != 6*mgqpt[0]*sizeof(double))
    fprintf(stderr," *** ERROR: Inconsistent restart file 5.5\n");
}

void
BelytschkoTsayShell::readHistory(int fn)
{
  // ---------------------------------------------------------------
  // read history variables from file for restart
  // ---------------------
  int readSize;

  readSize = read(fn, evar1, 5*mgqpt[0]*sizeof(double));
  if(readSize != 5*mgqpt[0]*sizeof(double))
    fprintf(stderr," *** ERROR: Inconsistent restart file 5.1\n");

  readSize = read(fn, evar2, 5*mgqpt[0]*sizeof(double));
  if(readSize != 5*mgqpt[0]*sizeof(double))
    fprintf(stderr," *** ERROR: Inconsistent restart file 5.2\n");

  readSize = read(fn, evoit1, 6*mgqpt[0]*sizeof(double));
  if(readSize != 6*mgqpt[0]*sizeof(double))
    fprintf(stderr," *** ERROR: Inconsistent restart file 5.3\n");

  readSize = read(fn, evoit2, 6*mgqpt[0]*sizeof(double));
  if(readSize != 6*mgqpt[0]*sizeof(double))
    fprintf(stderr," *** ERROR: Inconsistent restart file 5.4\n");

  readSize = read(fn, evoit3, 6*mgqpt[0]*sizeof(double));
  if(readSize != 6*mgqpt[0]*sizeof(double))
    fprintf(stderr," *** ERROR: Inconsistent restart file 5.5\n");
}

void
BelytschkoTsayShell::Elefintbt1(double delt, double *_ecord, double *_edisp, double *_evelo,
                                double trac[3], double tmftval, double *_efint)
{
#ifdef USE_EIGEN3
  //=======================================================================
  //  elefintbt1 = compute internal force matrix for bt shell, including hourglass force
  //
  //               note:
  //               ----
  //               rotation projection for the drilling dof is considered
  //
  //  arguments description
  //  ---------------------
  //  input:
  //  -----
  //  delt : time increment
  // 
  //  ecord(3,4) : element nodal coordinate
  //
  //  edisp(nndof,4) : element nodal displacement data
  //
  //  evelo(nndof,4) : element nodal velocity: v_x, v_y, v_z, theta_x, theta_y, theta_z
  //
  //  output:
  //  ------
  //  efint(nndof*4,1) : element nodal internal force
  //
  // ======================================================================

  Map<Matrix<double,3,4,ColMajor> > ecord(_ecord);
  Map<Matrix<double,6,4,ColMajor> > edisp(_edisp);
  Map<Matrix<double,6,4,ColMajor> > evelo(_evelo); // TODO use stride here to skip theta_z
  Map<Matrix<double,3,8,ColMajor> > efint(_efint);

  // ====================================
  // local variable
  // ==============
  Matrix<double,3,3,ColMajor> locbvec; 
  Matrix<double,3,4,ColMajor> ecurn;
  Matrix<double,3,8,ColMajor> efintloc;
  Matrix<double,3,4,ColMajor> ecurnloc;
  Matrix<double,5,4,ColMajor> eveloloc;
  double bmat1pt[8];  // dimension(2,4)
  double bcmat1pt[8]; // dimension(2,4)
  double bsmat1pt[24]; // dimension(2,3,4)
  double area;
  double gamma[4];
  double zgamma;
  double ipstrndot[3];
  double tsstrndot[2];
  // ====================================

  // initialize
  efintloc = Matrix<double,3,8,ColMajor>::Zero();

  // get current nodal coordinate
  ecurn = ecord + edisp.block<3,4>(0,0);

  // compute co rotational local base vector
  // ---------------------------------------
  _FORTRAN(getlocbvecbt)(ecurn.data(),
                         locbvec.data());
        // input : ecurn
        // output : locbvec

  // get local nodal coordinates and velocity
  // ----------------------------------------
  ecurnloc = locbvec.transpose()*ecurn;
  eveloloc.block<3,4>(0,0) = locbvec.transpose()*evelo.block<3,4>(0,0);
  eveloloc.block<2,4>(3,0) = locbvec.block<3,2>(0,0).transpose()*evelo.block<3,4>(3,0);

  // compute area
  // ------------
  area = 0.50*( (ecurnloc(0,2)-ecurnloc(0,0))*(ecurnloc(1,3)-ecurnloc(1,1)) 
               +(ecurnloc(0,1)-ecurnloc(0,3))*(ecurnloc(1,2)-ecurnloc(1,0)) );

  // check current element configuration
  // -----------------------------------
  if(area <= 0) {
     cerr << " *** ERROR: in BelytschkoTsayShell::Elefintbt1 current element has negative or zero area: " << area << endl;
     exit(-1);
  }

  //compute b matrix: b matrix, b^c matrix, and b^s matrix
  // ----------------
  // compute b matrix: one point integration
  _FORTRAN(getbmat1pt)(ecurnloc.data(), area, bmat1pt);
        // input : ecurnloc,area
        // output : bmat1pt

  // compute b^c matrix: warping correction
  _FORTRAN(getgamma4nod)(ecurnloc.data(), area, gamma, zgamma);
  if(optcor[0] > 0) {
     _FORTRAN(getbcmat1pt)(ecurnloc.data(), area, gamma, zgamma, bcmat1pt);
           // input : ecurnloc,area,gamma,zgamma
           // output : bcmat1pt
  }

  // compute b^s matrix: transverse shear projection
  if(optcor[1] > 0) {
    _FORTRAN(getbsmat1pt)(ecurnloc.data(), bsmat1pt);
      // input : ecurnloc
      // output : bsmat1pt
  }

  // loop over gauss quadrature
  // note: this loop is for through thickness integration
  //       in plane, we use 1 point rule
  for(int igaus = 0; igaus < mgaus[2]; ++igaus) {

    // compute rates of deformation and update strain at gq
    // ----------------------------------------------------
    _FORTRAN(updstrnbt)(optcor, delt, expmat->ematpro, gqpoin3[igaus], eveloloc.data(), bmat1pt, 
                        bcmat1pt, bsmat1pt, evoit3+6*igaus, ipstrndot, tsstrndot);
          // input : optcor,delt,ematpro,gqpoin3[igaus],eveloloc,bmat1pt,bcmat1pt,bsmat1pt
          // inoutput : strnvoitloc
          // output : ipstrndot, tsstrndot

    // update hypo stresses at gq
    // --------------------------
    _FORTRAN(updstrsbt)(expmat->optctv, optdmg, delt, expmat->ematpro, area, ipstrndot, tsstrndot,
                        evar1[5*igaus+0], evar1[5*igaus+1], evoit2+6*igaus, evoit3+6*igaus);
          // input : optctv,optdmg,delt,ematpro,area,ipstrndot,tsstrndot
          // inoutput : effpstrn,hardvar,sigvoitloc,strnvoitloc

    // update elasto plastic stresses at gq
    // ------------------------
    if(expmat->optctv != 1) {
      vector<double> F(9), CauchyStress(9);
      // get Fnp1 from strnvoitloc, i.e. evoit3[6*igaus+0]
      // note: voight rule in xfem code: [xx,yy,zz,yz,xz,xy]
      F[0] = 1+evoit3[6*igaus+0]; // xx
      F[1] = 0.5*evoit3[6*igaus+5]; // xy
      F[2] = 0.5*evoit3[6*igaus+4]; // xz
      F[3] = 0.5*evoit3[6*igaus+5]; // yx
      F[4] = 1+evoit3[6*igaus+1]; // yy
      F[5] = 0.5*evoit3[6*igaus+3]; // yz
      F[6] = 0.5*evoit3[6*igaus+4]; // zx
      F[7] = 0.5*evoit3[6*igaus+3]; // zy
      F[8] = 1+evoit3[6*igaus+2]; // zz
      if(!mat[igaus]->ComputeElastoPlasticConstitutiveResponse(F, &CauchyStress)) {
        cerr << " *** ERROR: ComputeElastoPlasticConstitutiveResponse failed\n";
        exit(-1);
      }
      // copy CauchyStress into sigvoitloc, i.e. evoit2[6*igaus+0]
      evoit2[6*igaus+0] = CauchyStress[0]; // xx
      evoit2[6*igaus+1] = CauchyStress[4]; // yy
      evoit2[6*igaus+2] = CauchyStress[8]; // zz
      //evoit2[6*igaus+3] = CauchyStress[5]; // yz
      //evoit2[6*igaus+4] = CauchyStress[2]; // xz
      evoit2[6*igaus+5] = CauchyStress[1]; // xy
    }

    // add local nodal internal forces at gq
    // -----------------------------------------
    _FORTRAN(gqfintbt)(optcor, delt, expmat->ematpro, gqpoin3[igaus], gqweigt3[igaus], area, 
                       evoit2+6*igaus, bmat1pt, bcmat1pt, bsmat1pt, efintloc.data());
          // input : optcor,delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),area,sigvoitloc,bmat1pt,bcmat1pt,bsmat1pt
          // inoutput : efintloc
  }

  // update hourglass control stresses
  // ---------------------------------
  if(opthgc >= 0) {
    _FORTRAN(updhgcstrsbt)(prmhgc, delt, expmat->ematpro, eveloloc.data(), area, bmat1pt,
                           gamma, zgamma, evoit1);
          // input : prmhgc,delt,ematpro,eveloloc,area,bmat1pt,gamma,zgamma
          // inoutput : hgcvoitloc
  }

  // add the hourglass control forces
  // --------------------------------------
  if(opthgc >= 0) {
    _FORTRAN(gqfhgcbt)(gamma, zgamma, evoit1, efintloc.data());
          // input : gamma,zgamma,hgcvoitloc
          // inoutput : efintloc
  }

  // subtract the local traction forces
  // -------------------------------------
  if(opttrc >= 0) {
    //cerr << "tmftval = " << tmftval << endl;
    _FORTRAN(elefbc3dbrkshl2opt)(area, trac, tmftval, efintloc.data());
          // input : area,trac,tmftval
          // inoutput : efintloc
  }

  // convert local efintloc to global efint
  // -------------------------------------
  efint = locbvec*efintloc;
#else
  cerr << "USE_EIGEN3 is not defined here in BelytschkoTsayShell::Elefintbt1\n";
  exit(-1);
#endif
}
