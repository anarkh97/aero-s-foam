// ---------------------------------------------------------------------
// HB - 05-22-05
// ---------------------------------------------------------------------
// 15 nodes wedge element 
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
#include <cstdio>
#include <iostream>
#include <Element.d/Penta15.d/Penta15.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Math.d/FullSquareMatrix.h>
#include <Corotational.d/Penta15Corotator.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NLPentahedral.h>
#include <Corotational.d/MatNLCorotator.h>

extern "C" {
void _FORTRAN(brkcmt)(double&, double&, double*);
}

void rotateConstitutiveMatrix(double *Cin, double *T33, double Cout[6][6]);
double Penta15ShapeFct(double Shape[15], double DShape[15][3], double m[3], double X[15], double Y[15], double Z[15]);
double computePenta15DShapeFct(double dShape[15][3], double X[15], double Y[15], double Z[15], double (*DShape)[3] = 0);
void addBtCBtoK3DSolid(FullSquareMatrix &K, double (*DShape)[3], double C[6][6], double alpha, int nnodes, int* ls);
void addNtDNtoM3DSolid(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls, double (*D)[3] = 0);
void computeStressAndEngStrain3DSolid(double Stress[6], double Strain[6], double C[6][6], double (*DShape)[3], double* U, int nnodes, int* ls=0);
double computeStress3DSolid(double Stress[6],double Strain[6], double C[6][6]);
double computeVonMisesStress(double Stress[6]);
double computeVonMisesStrain(double Strain[6]);
void Penta15ShapeFct(double Shape[15], double m[3]);

extern bool useFull;

double weight3d8[9] = { 0.092592592592593,0.092592592592593,0.092592592592593,
                        0.148148148148148,0.148148148148148,0.148148148148148,
                        0.092592592592593,0.092592592592593,0.092592592592593};
double gauss3d8[9][3] = { {0.166666666666667,0.166666666666667,-0.774596669241483},
                          {0.666666666666667,0.166666666666667,-0.774596669241483},
                          {0.166666666666667,0.666666666666667,-0.774596669241483},
                          {0.166666666666667,0.166666666666667,0.},
                          {0.666666666666667,0.166666666666667,0.},
                          {0.166666666666667,0.666666666666667,0.},
                          {0.166666666666667,0.166666666666667,0.774596669241483},
                          {0.666666666666667,0.166666666666667,0.774596669241483},
                          {0.166666666666667,0.666666666666667,0.774596669241483} };

Penta15::Penta15(int* nodenums)
{
  for(int i=0; i<15; i++)
    nn[i] = nodenums[i];

  cFrame = 0;
  cCoefs = 0;
  mat = 0;
}

Penta15::~Penta15()
{
  if(cCoefs && mat) delete mat;
}

Element *
Penta15::clone()
{
  return new Penta15(*this);
}

void
Penta15::renum(int *table)
{
  for(int i=0; i<15; i++)
    nn[i] = table[nn[i]];
}

void
Penta15::renum(EleRenumMap& table)
{
  for(int i=0; i<15; i++)
    nn[i] = table[nn[i]];
}

void
Penta15::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
                     Vector& elDisp, int strInd, int surface, double *ndTemps,
                     double ylayer, double zlayer, int avgnum)
{
  const int nnodes = 15;
  weight = 1.0;

  double X[15], Y[15], Z[15];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // Flags to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = false; // This can be used to force averaging of the Von Mises stress & strain.

  double elStress[15][7];
  double elStrain[15][7];

  // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[15][3] = {{  0., 0., -1.},{  1.,  0., -1.},{ 0.,  1., -1.},
                                {  0., 0.,  1.},{  1.,  0.,  1.},{ 0.,  1.,  1.},
                                { 0.5, 0., -1.},{ 0.5, 0.5, -1.},{ 0., 0.5, -1.},
                                { 0.5, 0.,  1.},{ 0.5, 0.5,  1.},{ 0., 0.5,  1.},
                                {  0., 0.,  0.},{  1.,  0.,  0.},{ 0.,  1.,  0.}};

  double Shape[15], DShape[15][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta15ShapeFct(Shape, DShape, m, X, Y, Z); 
    computeStressAndEngStrain3DSolid(elStress[inode],elStrain[inode], C, DShape, elDisp.data(), nnodes);

    if(ndTemps) {
      double &Tref  = prop->Ta;
      double &alpha = prop->W;
      double eT     = alpha*(ndTemps[inode]-Tref);
      double thermalStrain[6] = {eT,eT,eT,0.0,0.0,0.0};
      double thermalStress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
      computeStress3DSolid(thermalStress, thermalStrain, C);
      elStress[inode][0] -= thermalStress[0];
      elStress[inode][1] -= thermalStress[1];
      elStress[inode][2] -= thermalStress[2];
      elStress[inode][3] -= thermalStress[3];
      elStress[inode][4] -= thermalStress[4];
      elStress[inode][5] -= thermalStress[5];
    }
    
    if(vmflg) elStress[inode][6] = computeVonMisesStress(elStress[inode]);
    else elStress[inode][6] = 0.0;
 
    if(strainFlg) elStrain[inode][6] = computeVonMisesStrain(elStrain[inode]);
    else elStrain[inode][6] = 0.0;
  }

  // compute average Von Mises stress and/or Von Mises strain: to match old Fortran code 
  if(vmflg && meanVms) {
    double vms = 0.0;
    for(int inode=0; inode<nnodes; inode++) vms += elStress[inode][6];
    vms /= nnodes;
    for(int inode=0; inode<nnodes; inode++) elStress[inode][6] = vms;
  }
  if(strainFlg && meanVms) {
    double vms = 0.0;
    for(int inode=0; inode<nnodes; inode++) vms += elStrain[inode][6];
    vms /= nnodes;
    for(int inode=0; inode<nnodes; inode++) elStrain[inode][6] = vms;
  }
  
  // fill the output array stress with the requested stress or strain component 
  for(int inode=0; inode<nnodes; inode++) {
    if(strInd < 7) 
      stress[inode] = elStress[inode][strInd];
    else if(strInd < 14)
      stress[inode] = elStrain[inode][strInd-7];
    else
      stress[inode] = 0;
  }
}

void
Penta15::getAllStress(FullM& stress, Vector& weight, CoordSet &cs,
                      Vector& elDisp, int strInd,int surface, double *ndTemps)
{
  const int nnodes = 15;
  weight = 1.0;
  
  double X[15], Y[15], Z[15];
  cs.getCoordinates(nn, nnodes, X, Y, Z);
      
  double elStress[15][7];
  double elStrain[15][7];

  // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);

  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[15][3] = {{  0., 0., -1.},{  1.,  0., -1.},{ 0.,  1., -1.},
                                {  0., 0.,  1.},{  1.,  0.,  1.},{ 0.,  1.,  1.},
                                { 0.5, 0., -1.},{ 0.5, 0.5, -1.},{ 0., 0.5, -1.},
                                { 0.5, 0.,  1.},{ 0.5, 0.5,  1.},{ 0., 0.5,  1.},
                                {  0., 0.,  0.},{  1.,  0.,  0.},{ 0.,  1.,  0.}};

  double Shape[15], DShape[15][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta15ShapeFct(Shape, DShape, m, X, Y, Z); 
    computeStressAndEngStrain3DSolid(elStress[inode], elStrain[inode], C, DShape, elDisp.data(), nnodes);

    if(ndTemps) {
      double &Tref  = prop->Ta;
      double &alpha = prop->W;
      double eT     = alpha*(ndTemps[inode]-Tref);
      double thermalStrain[6] = {eT,eT,eT,0.0,0.0,0.0};
      double thermalStress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
      computeStress3DSolid(thermalStress, thermalStrain, C);
      elStress[inode][0] -= thermalStress[0];
      elStress[inode][1] -= thermalStress[1];
      elStress[inode][2] -= thermalStress[2];
      elStress[inode][3] -= thermalStress[3];
      elStress[inode][4] -= thermalStress[4];
      elStress[inode][5] -= thermalStress[5];
    }
  }
  // Store all Stress or all Strain as defined by strInd
  if(strInd == 0) {
    for(int i=0; i<nnodes; ++i)
      for(int j=0; j<6; ++j) 
        stress[i][j] = elStress[i][j];
  } else {
    for(int i=0; i<nnodes; ++i)
      for(int j=0; j<6; ++j) 
        stress[i][j] = elStrain[i][j];
  }       
  
  // Get Element Principals for each node without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};

  for(int i=0; i<nnodes; ++i) {
    for(int j=0; j<6; ++j) {
      svec[j] = stress[i][j];
    }
    // Convert Engineering to Tensor Strains
    if(strInd != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec,pvec);
    for(int j=0; j<3; ++j) {
      stress[i][j+6] = pvec[j];
    }
  }
}

double
Penta15::getMass(CoordSet& cs)
{
  const int nnodes = 15;
  double X[15], Y[15], Z[15];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // integration: loop over gauss pts
  const int ngauss = 9;
  double shape[nnodes];
  double dShape[nnodes][3];
  double volume = 0.0;

  for(int i = 0; i < ngauss; i++) {
    double dOmega = Penta15ShapeFct(shape, dShape, gauss3d8[i], X, Y, Z);
    volume += fabs(dOmega)*weight3d8[i];
  }

  return volume*prop->rho;
}

void
Penta15::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                         Vector& gravityForce, int gravflg, GeomState *geomState)
{
  int nnodes = 15;

  // Lumped
  if (gravflg != 2) {

    double totmas = getMass(cs);

    // divvy up the total body force using same ratio as the corresponding diagonal of the lumped mass matrix to the total mass
    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = totmas*gravityAcceleration[0]*(3.0*factors[3*i+0]);
      gravityForce[3*i+1] = totmas*gravityAcceleration[1]*(3.0*factors[3*i+1]);
      gravityForce[3*i+2] = totmas*gravityAcceleration[2]*(3.0*factors[3*i+2]);
    }
  }
  // Consistent
  else {
    double x[15], y[15], z[15];
    cs.getCoordinates(nn, numNodes(), x, y, z);

    double lforce[15];
    for(int i = 0; i < nnodes; ++i) lforce[i] = 0.0;

    // integration: loop over Gauss pts
    const int ngauss = 9;
    double shape[15], dShape[15][3];

    for(int i = 0; i < ngauss; i++) {
      double dOmega = Penta15ShapeFct(shape, dShape, gauss3d8[i], x, y, z);
      for(int j = 0; j < nnodes; ++j) lforce[j] += fabs(dOmega)*weight3d8[i]*shape[j];
    }

    for(int i = 0; i < nnodes; ++i)
      for(int j = 0; j < 3; ++j)
        gravityForce[3*i+j] = lforce[i]*gravityAcceleration[j]*prop->rho;

  }
}

void
Penta15::getThermalForce(CoordSet &cs, Vector &ndTemps,
                         Vector &elementThermalForce, int glflag,
                         GeomState *geomState)
{
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes = 15;
  const int ndofs = 45;

  // initialize nodal thermal forces
  for(int i=0; i<ndofs; i++) elementThermalForce[i] = 0.0;

  // for nonlinear analyses, the thermal load for this element is now computed in getStiffAndForce
  if(geomState) return;

  double X[15], Y[15], Z[15];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // get material props & constitutive matrix
  double &Tref  = prop->Ta;
  double &alpha = prop->W ;
  double C[6][6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);

  // Integate over the element: F = Int[Bt.ThermaStress]
  // with ThermalStress = C.ThermalStrain, with ThermalStrain = alpha.theta.[1, 1, 1, 0, 0, 0]'
  // where theta = T(M)-Tref = Sum[inode][N[inode]*(ndTemps[inode] - Tref)]
  // N[inode] is the shape fct at node inode
  // M is the position in the real frame, m its associated position in the reference
  // element frame
  // NUMERICAL INTEGRATION BY GAUSS PTS
  int ngauss = 9;
  double Shape[15], DShape[15][3];
  double w, J;

  for(int i = 0; i < ngauss; i++) {
    // compute shape fcts & their derivatives at the Gauss pt
    J = Penta15ShapeFct(Shape, DShape, gauss3d8[i], X, Y, Z);
    w = weight3d8[i]*fabs(J);
    // compute thermal stresses
    double eT = 0.0;
    for(int inode=0; inode<nnodes; inode++) eT += alpha*Shape[inode]*(ndTemps[inode] - Tref);
    double thermalStrain[6] = {eT,eT,eT,0.0,0.0,0.0};
    double thermalStress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    computeStress3DSolid(thermalStress, thermalStrain, C); // thermalStress <- C.thermalStrain
    // sum contribution
    for(int inode=0; inode<nnodes; inode++) {
      elementThermalForce[3*inode  ] += w*(DShape[inode][0]*thermalStress[0] + DShape[inode][1]*thermalStress[3] + DShape[inode][2]*thermalStress[5]);
      elementThermalForce[3*inode+1] += w*(DShape[inode][0]*thermalStress[3] + DShape[inode][1]*thermalStress[1] + DShape[inode][2]*thermalStress[4]);
      elementThermalForce[3*inode+2] += w*(DShape[inode][0]*thermalStress[5] + DShape[inode][1]*thermalStress[4] + DShape[inode][2]*thermalStress[2]);
    }
  }
}

FullSquareMatrix
Penta15::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  const int nnodes = 15;
  const int ndofs = 45;

  FullSquareMatrix M(ndofs, mel);

  if(cmflg) { // consistent mass matrix

    double x[15], y[15], z[15];
    cs.getCoordinates(nn, nnodes, x, y, z);

    M.zero();
    int ls[45] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,
                  1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,
                  2,5,8,11,14,17,20,23,26,29,32,35,38,41,44};

    // integration: loop over Gauss pts
    const int ngauss = 9;
    double shape[nnodes];
    double dShape[nnodes][3];

    for(int i = 0; i < ngauss; i++) {
      double dOmega = Penta15ShapeFct(shape, dShape, gauss3d8[i], x, y, z);
      double w = fabs(dOmega)*weight3d8[i]*prop->rho;
      addNtDNtoM3DSolid(M, shape, w, nnodes, ls);
    }
  }
  else { // lumped mass matrix
    fprintf(stderr," *** In Penta15::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return M;
}

FullSquareMatrix
Penta15::stiffness(CoordSet &cs, double *d, int flg)
{
  const int nnodes = 15;
  const int ndofs = 45;

  double X[15], Y[15], Z[15];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);

  FullSquareMatrix K(ndofs, d);
  K.zero();

  int ls[45] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,
                1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,
                2,5,8,11,14,17,20,23,26,29,32,35,38,41,44};

  // integration: loop over Gauss pts
  const int ngauss = 9;
  double shape[nnodes];
  double dShape[nnodes][3];

  for(int i = 0; i < ngauss; i++) {
    double dOmega = Penta15ShapeFct(shape, dShape, gauss3d8[i], X, Y, Z);
    double w = fabs(dOmega)*weight3d8[i];
    addBtCBtoK3DSolid(K, dShape, C, w, nnodes, ls);
  }

  return K;
}

int
Penta15::numNodes()
{ 
  if(useFull)
    return 15; 
  else
    return 6;
}

int
Penta15::numDofs()
{
  return 45;
}

// treat as a 6-node penta
// this is because xpost does not have a 15-node penta
int
Penta15::getTopNumber()
{
  return 124;
}

int
Penta15::numTopNodes()
{
  return 6;
}

int*
Penta15::nodes(int *p)
{
  if(useFull) {
    if(!p) p = new int[15];
    for(int i=0; i<15; i+=3) {
      p[i] = nn[i]; p[i+1] = nn[i+1]; p[i+2] = nn[i+2];
    }
  }
  else {
    if(!p) p = new int[6];
    for(int i=0; i<6; i++)
      p[i] = nn[i];
  }

  return p;
}

int*
Penta15::dofs(DofSetArray &dsa, int *p)
{
  if(!p) p = new int[45];

  for(int i=0; i<15; i++)
    dsa.number(nn[i], DofSet::XYZdisp, p+3*i);

  return p;
}

void
Penta15::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

void
Penta15::setMaterial(NLMaterial *_mat)
{
  if(cCoefs) { // anisotropic material
    double C[6][6];
    mat = _mat->clone();
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
    if(mat) mat->setTangentMaterial(C);
    Eigen::Map<Eigen::Matrix<double,6,6,Eigen::RowMajor> > _C(&C[0][0]);
  }
  else {
    mat = _mat;
  }
}

int
Penta15::numStates()
{
  int numGaussPoints = 9;
  return (mat) ? numGaussPoints*mat->getNumStates() : 0;
}

Corotator *
Penta15::getCorotator(CoordSet &cs, double *kel, int, int)
{
  if(cCoefs && !mat) {
    double C[6][6];
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
    mat = new StVenantKirchhoffMat(prop->rho, C, prop->Ta, prop->W);
  }
  if(mat) {
#ifdef USE_EIGEN3
    MatNLElement *ele = new NLPentahedral15(nn);
    ele->setMaterial(mat);
    ele->setGlNum(glNum);
    ele->setProp(prop);
    return new MatNLCorotator(ele);
#endif
  }
  else {
    return new Penta15Corotator(nn, prop->E, prop->nu, cs, prop->Ta, prop->W);
  }
  printf("WARNING: Corotator not implemented for element %d\n", glNum+1); return 0;
}
