#include <Element.d/BelytschkoTsayShell.d/BelytschkoTsayShell.h>
#include <Utils.d/dbg_alloca.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Element.d/NonLinearity.d/ExpMat.h>

extern "C" {
  void _FORTRAN(getgqsize)(int&, int&, int*, int*, int*);
  void _FORTRAN(getgq1d)(int&, double*, double*);
  void _FORTRAN(gethistvar3d)(int&, int&, int&, double*, double*, double*, double*, double*, double*, double*);
  void _FORTRAN(elefinthypobt1)(int&, int*, double*, double&, double*, int&, int&, int&,
                                double*, double*, double*, double*, double*, double*,
                                double*, double*, double*, double*);
  void _FORTRAN(elefintevpbt1)(int*, double&, double*, int&, int&, int&, double*, double*,
                               double*, double*, double*, double*, double*, double*, double*, double*, double*);
  void _FORTRAN(elefintj2bt1)(int*, double&, double*, int&, int&, int&, double*, double*,
                              double*, double*, double*, double*, double*, double*, double*, double*);
  void _FORTRAN(elefhgcbt1)(double*, double&, double*, int&, double*, double*, double*,
                            double*, double*);
  void _FORTRAN(elemaslbt)(int&, double*, double*, double*, double*);
  void _FORTRAN(elefbcbt2)(int&, int&, int&, int&, double*, double*, double*, double*);
}

BelytschkoTsayShell::BelytschkoTsayShell(int* nodenums)
{
  optele = 3; // bt shell
  optmhd = 0; // conventional fem
  optdmg = 0; // no damage
  opthgc = 1; // perturbation type hourglass control
  optcri[0] = optcri[1] = 0; // criterion type = none, averaging type = none
  optdmp = 0; // no damping
  nndof  = 6; // number of dofs per node
  ndime  = 3;
  nnode  = 4;
  ngqpt[0] = 1; ngqpt[1] = 1; ngqpt[2] = 3;
  ngqpt4 = 3;
  // damage control parameters
  for(int i = 0; i < 10; ++i) prmdmg[i];
  // hourglass control parameters
  prmhgc[0] = 5e-2;
  prmhgc[1] = 5e-3;
  prmhgc[2] = 5e-2;
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
  _FORTRAN(getgqsize)(optmhd, optele, ngqpt, mgaus, mgqpt);
     // input : optmhd,optele,ngqpt
     // output : mgaus,mgqpt

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
}

BelytschkoTsayShell::~BelytschkoTsayShell()
{
  delete [] evar1;
  delete [] evar2;
  delete [] evoit1;
  delete [] evoit2;
  delete [] evoit3;
}

/*
void
BelytschkoTsayShell::setProp(StructProp *p, bool _myProp)
{
  Element::setProp(p, _myProp);
  ematpro[0] = prop->E;    // Young's modulus
  ematpro[1] = prop->nu;   // Poisson's ratio
  ematpro[2] = prop->rho;  // mass density
  ematpro[3] = 0;          // yield stress (j2) or reference strain rate (evp)
  ematpro[4] = 0;          // hardening modulus (j2) or rate sensitivity parameter (evp)
  ematpro[5] = 0;          // yield stress (evp)
  ematpro[6] = 0;          // yield strain (evp)
  ematpro[7] = 0;          // strain hardening exponent (evp)
  ematpro[18] = 0.833;     // shear correction factor
  ematpro[19] = prop->eh;  // thickness
}
*/

void
BelytschkoTsayShell::setMaterial(NLMaterial *m)
{
  expmat = dynamic_cast<ExpMat *>(m);
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
  // voight rule in xfem code: [xx,yy,zz,yz,xz,xy]
  stress.zero();
  weight = 1.0;
  int k[6] = { 0, 1, 2, 5, 3, 4 };
  for(int i = 0; i < nnode; ++i) {
    for(int j = 0; j < mgqpt[0]; ++j) {
      switch(strInd) {
        case 0 : case 1 : case 2 : case 3 : case 4 : case 5 : // sxx, syy, szz, sxy, syz, sxz
          stress[i] += evoit2[6*j+k[strInd]];
          break;
        case 6 : // von mises stress
          stress[i] += evar2[5*j+0];
          break;
        case 7 : case 8 : case 9 : case 10: case 11: case 12: // exx, eyy, ezz, exy, eyz, exz
          stress[i] += evoit1[6*j+k[strInd-7]];
          break;
        case 13 : // von mises strain
          stress[i] += evar1[5*j+0];
          break;
        case 14 : case 15 : case 16 : // not available
          break;
        case 17 : // damage
          stress[i] += evar1[5*j+1];
          break;
      }
    }
  }
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

  double r1[3], r2[3], r3[3], v1[3], v2[3], v3[3];

  r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
  r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
  r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

  v1[0] = r3[0] - r1[0];
  v1[1] = r3[1] - r1[1];
  v1[2] = r3[2] - r1[2];

  v2[0] = r2[0] - r1[0];
  v2[1] = r2[1] - r1[1];
  v2[2] = r2[2] - r1[2];

  crossprod(v1, v2, v3);

  double area = 0.5*sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
  double mass = area*prop->rho*prop->eh;

  return mass;
}

void
BelytschkoTsayShell::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                                     Vector& gravityForce, int gravflg, GeomState *geomState)
{
  if(prop) {
    int optdom = 0; // 0 : integrate over reference domain (TODO: check)
                    // 1 : integrate over current domain
    int opttrc = 1; // 0 : pressure
                    // 1 : traction
    double* ecord = (double*) dbg_alloca(sizeof(double)*nnode*ndime);
    double* edisp = 0;
    for(int i = 0; i < nnode; ++i) {
      ecord[i*ndime+0] = cs[nn[i]]->x;
      ecord[i*ndime+1] = cs[nn[i]]->y;
      ecord[i*ndime+2] = cs[nn[i]]->z;
    }
    double trac[3];
    for(int i = 0; i < 3; ++i) trac[i] = gravityAcceleration[i]*prop->rho*prop->eh;
    double *efbc = gravityForce.data();

    _FORTRAN(elefbcbt2)(optdom, opttrc, nndof, ngqpt4, ecord, edisp, trac, efbc);
  }
  else gravityForce.zero();
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
      for(int j = 0; j < 6; ++j) edisp[i*nndof+j] = 0; // constant mass matrix
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
BelytschkoTsayShell::getStiffAndForce(GeomState& geomState, CoordSet& cs, FullSquareMatrix& k, double* efint, double delt)
{
  //=======================================================================
  //  compute internal force vector for Belytschko Tsay shell element
  // 
  //  output:
  //  ------
  //  efint(24) : internal force
  //
  // ======================================================================
  k.zero();

  if(prop) { // check for phantom
    // ---------------------------------------------------------------
    // get element nodal displacement and velocity
    // -------------------------------------------
    double* ecord = (double*) dbg_alloca(sizeof(double)*nnode*ndime);
    double* edisp = (double*) dbg_alloca(sizeof(double)*nnode*nndof); // d^{n+1}
    double* evelo = (double*) dbg_alloca(sizeof(double)*nnode*nndof); // v^{n+0.5}
    double* eaccl = (double*) dbg_alloca(sizeof(double)*nnode*nndof); // a^{n}
    for(int i = 0; i < nnode; ++i) {
      ecord[i*ndime+0] = cs[nn[i]]->x;
      ecord[i*ndime+1] = cs[nn[i]]->y;
      ecord[i*ndime+2] = cs[nn[i]]->z;
      edisp[i*nndof+0] = geomState[nn[i]].x - cs[nn[i]]->x;
      edisp[i*nndof+1] = geomState[nn[i]].y - cs[nn[i]]->y;
      edisp[i*nndof+2] = geomState[nn[i]].z - cs[nn[i]]->z;
      mat_to_vec(geomState[nn[i]].R, edisp+i*nndof+3);
      for(int j = 0; j < 6; ++j) evelo[i*nndof+j] = geomState[nn[i]].v[j] + delt*0.5*geomState[nn[i]].a[j];
      for(int j = 0; j < 6; ++j) eaccl[i*nndof+j] = geomState[nn[i]].a[j];
    }

    // ---------------------------------------------------------------
    // get 1d gq rule: through thickness integration
    // --------------
    double* gqpoin3 = (double*) dbg_alloca(sizeof(double)*mgaus[2]);
    double* gqweigt3 = (double*) dbg_alloca(sizeof(double)*mgaus[2]);
    _FORTRAN(getgq1d)(mgaus[2], gqpoin3, gqweigt3);
       // input : mgaus[2]
       // output : gqpoin3,gqweigt3

    // ---------------------------------------------------------------
    // constitutive model
    // ------------------
    switch(expmat->optctv) {
      default:
      case 1 : // hypoelastic
        _FORTRAN(elefinthypobt1)(optdmg, optcri, prmdmg, delt, expmat->ematpro, nndof, mgaus[2], mgqpt[0],
                                 gqpoin3, gqweigt3, ecord, edisp, evelo, evar1, evoit2, evoit3,
                                 evar2, efint);
           // input : optdmg,optcri,prmdmg
           //         delt,ematpro,nndof,mgaus[2],mgqpt[0],gqpoin3,gqweigt3,ecord,edisp,evelo
           // inoutput : evar1(1.eff.strn, 2.damage),evoit2(cauchy, local), evoit3(strain)
           // output : evar2,efint
        break;
      case 2 : // thermo elasto viscoplastic 
        cerr << " *** WARNING: thermo elasto viscoplastic constitutive model not supported for bt shell element"
             << ", using elasto viscoplastic instead\n";
      case 3 : // elasto viscoplastic
        _FORTRAN(elefintevpbt1)(optcri, delt, expmat->ematpro, nndof, mgaus[2], mgqpt[0], gqpoin3, gqweigt3,
                                ecord, edisp, evelo, eaccl, evar1, evoit2, evoit3, evar2, efint);
           // input : optcri,delt,ematpro,nndof,mgaus[2],mgqpt[0],gqpoin3,gqweigt3,ecord,edisp,evelo,eaccl
           // inoutput : evar1,evoit2,evoit3
           // output : evar2,efint
        break;
      case 4: // gurson model
        cerr << " *** WARNING: gurson constitutive model not supported for bt shell element, using j2 instead\n";
      case 5: // j2 explicit
        _FORTRAN(elefintj2bt1)(optcri, delt, expmat->ematpro, nndof, mgaus[2], mgqpt[0], gqpoin3, gqweigt3,
                               ecord, edisp, evelo, evar1, evoit2, evoit3, evar2, efint);
           // input : optcri,delt,ematpro,nndof,mgaus(3),mgqpt[0],gqpoin3,gqweigt3,ecord,edisp,evelo
           // inoutput : evar1,evoit2,evoit3
           // output : evar2,efint
        break;
    }
  
    // ---------------------------------------------------------------
    // hourglass control force
    // -----------------------
    double* efhgc = (double*) dbg_alloca(sizeof(double)*nnode*nndof);
    for(int i = 0; i < nnode*nndof; ++i) efhgc[i] = 0;
    _FORTRAN(elefhgcbt1)(prmhgc, delt, expmat->ematpro, nndof, ecord, edisp, evelo, evoit1, efhgc);
       // input : prmhgc,delt,ematpro,nndof,ecord,edisp,evelo
       // inoutput : evoit1(hgc local)
       // output : efhgc
  
    // add hourglass control force to internal force
    for(int i = 0; i < nnode*nndof; ++i) efint[i] += efhgc[i];

    // ---------------------------------------------------------------
    // damping force
    // -----------------------
    if(optdmp) {
      double dmpfctr = prmdmp[0]; // damping factor: 0.0-1.0
      double freq = prmdmp[1]; // frequency: 1.0d0 / cycle(time)
      double cnst = 2.0 * freq * dmpfctr;
      double* emasl = (double*) dbg_alloca(sizeof(double)*nnode*nndof);
      _FORTRAN(elemaslbt)(nndof, expmat->ematpro, ecord, edisp, emasl); // TODO: unnecessary recomputation
      for(int i = 0; i < nnode*nndof; ++i) efint[i] += cnst*emasl[i]*evelo[i];
    }
  }
}

void
BelytschkoTsayShell::computeDisp(CoordSet&, State &state, const InterpPoint &ip,
                                 double *res, GeomState *gs)
{
  cerr << "BelytschkoTsayShell::computeDisp not implemented\n";
}

void
BelytschkoTsayShell::getFlLoad(CoordSet &, const InterpPoint &ip, double *flF, 
                               double *resF, GeomState *gs)
{
  cerr << "BelytschkoTsayShell::getFlLoad not implemented\n";
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
  int optdom = 1; // 0 : integrate over reference domain
                  // 1 : integrate over current domain
  int opttrc = 0; // 0 : pressure
                  // 1 : traction
  double* ecord = (double*) dbg_alloca(sizeof(double)*nnode*ndime);
  double* edisp = (double*) dbg_alloca(sizeof(double)*nnode*nndof);
  for(int i = 0; i < nnode; ++i) {
    ecord[i*ndime+0] = cs[nn[i]]->x;
    ecord[i*ndime+1] = cs[nn[i]]->y;
    ecord[i*ndime+2] = cs[nn[i]]->z;
    edisp[i*nndof+0] = (*geomState)[nn[i]].x - cs[nn[i]]->x;
    edisp[i*nndof+1] = (*geomState)[nn[i]].y - cs[nn[i]]->y;
    edisp[i*nndof+2] = (*geomState)[nn[i]].z - cs[nn[i]]->z;
    mat_to_vec((*geomState)[nn[i]].R, edisp+i*nndof+3);
  }
  double trac[3] = { pressure, 0, 0 };
  double *efbc = elPressureForce.data();

  _FORTRAN(elefbcbt2)(optdom, opttrc, nndof, ngqpt4, ecord, edisp, trac, efbc);
}

void
BelytschkoTsayShell::getThermalForce(CoordSet& cs, Vector& ndTemps,
                                     Vector &elThermalForce, int glflag, 
                                     GeomState *geomState)
{
  cerr << "BelytschkoTsayShell::getThermalForce not implemented\n";
}
