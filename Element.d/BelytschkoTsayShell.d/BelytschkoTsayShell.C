#include <Element.d/BelytschkoTsayShell.d/BelytschkoTsayShell.h>
#include <Utils.d/dbg_alloca.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Element.d/NonLinearity.d/HypoElasticMat.h>
#include <Element.d/NonLinearity.d/ElastoViscoPlasticMat.h>
#include <Element.d/NonLinearity.d/J2PlasticityMat.h>

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
  optctv = 1; // hypoelastic constitutive law
  optdmg = 0; // no damage
  opthgc = 1; // perturbation type hourglass control
  nndof  = 6; // number of dofs per node
  ndime  = 3;
  nnode  = 4;
  ngqpt[0] = 1; ngqpt[1] = 1; ngqpt[2] = 2;
  ngqpt4 = 2;
  for(int i = 0; i < 10; ++i) prmdmg[i] = prmhgc[i] = optcri[i] = 0;

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
}

BelytschkoTsayShell::~BelytschkoTsayShell()
{
  delete [] evar1;
  delete [] evar2;
  delete [] evoit1;
  delete [] evoit2;
  delete [] evoit3;
}

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
  ematpro[18] = 0.84;      // shear correction factor
  ematpro[19] = prop->eh;  // thickness
}

void
BelytschkoTsayShell::setMaterial(NLMaterial *m)
{
  HypoElasticMat *m1 = dynamic_cast<HypoElasticMat *>(m);
  if(m1) {
    for(int i = 0; i < 3; ++i) ematpro[i] = m1->ematpro[i];
    optctv = 1;
    return;
  }
  ElastoViscoPlasticMat *m3 = dynamic_cast<ElastoViscoPlasticMat *>(m);
  if(m3) {
    for(int i = 0; i < 8; ++i) ematpro[i] = m3->ematpro[i];
    optctv = 3;
    for(int i = 0; i < mgqpt[0]; ++i) evar1[5*i+1] = 293.0; // TODO check
    return;
  }
  J2PlasticityMat *m5 = dynamic_cast<J2PlasticityMat *>(m);
  if(m5) {
    for(int i = 0; i < 5; ++i) ematpro[i] = m5->ematpro[i];
    optctv = 5;
    return;
  }
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
  cerr << "BelytschkoTsayShell::getVonMises not implemented\n";
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
  cerr << "BelytschkoTsayShell::getMass not implemented\n";
  return 0.0;
}

void
BelytschkoTsayShell::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                                     Vector& gravityForce, int gravflg, GeomState *geomState)
{
  cerr << "BelytschkoTsayShell::getGravityForce not implemented\n";
}

FullSquareMatrix
BelytschkoTsayShell::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
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
  _FORTRAN(elemaslbt)(nndof, ematpro, ecord, edisp, emasl);
     // input : nndof,ematpro,ecord,edisp
     // output : emasl
  delete [] ecord;
  delete [] edisp;
  
  FullSquareMatrix ret(nnode*nndof, mel);
  ret.zero();
  for(int i = 0; i < nnode*nndof; ++i) ret[i][i] = emasl[i];

  return ret;
}

FullSquareMatrix
BelytschkoTsayShell::stiffness(CoordSet &cs, double *d, int flg)
{
  cerr << "BelytschkoTsayShell::stiffness not implemented\n";
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
/*
  double ematpro[20]; // material properties
  for(int i = 0; i < 20; ++i) ematpro[i] = 0;
  ematpro[0] = prop->E;    // Young's modulus
  ematpro[1] = prop->nu;   // Poisson's ratio
  ematpro[2] = prop->rho;  // mass density
  ematpro[3] = prop->sigE; // yield stress (j2) or reference strain rate (evp)
  ematpro[4] = prop->Ep;   // hardening modulus (j2) or rate sensitivity parameter (evp)
  ematpro[5] = 0;          // yield stress (evp)
  ematpro[6] = 0;          // yield strain (evp)
  ematpro[7] = 0;          // strain hardening exponent (evp)
  ematpro[18] = 0.84;      // shear correction factor
  ematpro[19] = prop->eh;  // thickness
*/
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
  switch(optctv) {
    default:
    case 1 : // hypoelastic
      _FORTRAN(elefinthypobt1)(optdmg, optcri, prmdmg, delt, ematpro, nndof, mgaus[2], mgqpt[0],
                               gqpoin3, gqweigt3, ecord, edisp, evelo, evar1, evoit2, evoit3,
                               evar2, efint);
         // input : optdmg,optcri,prmdmg
         //         delt,ematpro,nndof,mgaus[2],mgqpt[0],gqpoin3,gqweigt3,ecord,edisp,evelo
         // inoutput : evar1(1.eff.strn, 2.damage),evoit2(cauchy, local), evoit3(strain)
         // output : evar2,efint
      break;
    case 2 : // thermo elasto viscoplastic 
      cerr << " *** WARNING: thermo elasto viscoplastic constitutive model not supported for bt shell element, using elasto viscoplastic instead\n";
    case 3 : // elasto viscoplastic
      _FORTRAN(elefintevpbt1)(optcri, delt, ematpro, nndof, mgaus[2], mgqpt[0], gqpoin3, gqweigt3,
                              ecord, edisp, evelo, eaccl, evar1, evoit2, evoit3, evar2, efint);
         // input : optcri,delt,ematpro,nndof,mgaus[2],mgqpt[0],gqpoin3,gqweigt3,ecord,edisp,evelo,eaccl
         // inoutput : evar1,evoit2,evoit3
         // output : evar2,efint
      break;
    case 4: // gurson model
      cerr << " *** WARNING: gurson constitutive model not supported for bt shell element, using j2 instead\n";
    case 5: // j2 explicit
      _FORTRAN(elefintj2bt1)(optcri, delt, ematpro, nndof, mgaus[2], mgqpt[0], gqpoin3, gqweigt3,
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
  _FORTRAN(elefhgcbt1)(prmhgc, delt, ematpro, nndof, ecord, edisp, evelo, evoit1, efhgc);
     // input : prmhgc,delt,ematpro,nndof,ecord,edisp,evelo
     // inoutput : evoit1(hgc local)
     // output : efhgc

  // add hourglass control force to internal force
  for(int i = 0; i < nnode*nndof; ++i) efint[i] += efhgc[i];
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
