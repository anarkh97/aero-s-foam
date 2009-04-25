// ---------------------------------------------------------------------
// HB - 05-22-05
// ---------------------------------------------------------------------
// 15 nodes wedge element 
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
// EXPERIMENTAL ...
// ---------------------------------------------------------------------
#include <stdio.h>
#include <iostream>
#include <Element.d/Penta15.d/Penta15.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Corotational.d/PentaCorotator.h>

#define CHECK_JACOBIAN //HB: force check nullity & constant sign of jacobian over el.
//#define PENTA15_DEBUG

extern bool useFull;

extern "C" {
void  _FORTRAN(brkcmt)(double&, double&, double*);

void _FORTRAN(lgauss)(const int &, int &, double *, double *);
}

//HB: for anisotropic elastic constitutive matrix
void rotateConstitutiveMatrix(double *Cin, double *T33, double Cout[6][6]);

// HB: for stiffness matrix with ansitropic constitutive matrix and/or consistent mass matrix
extern void   Penta15ShapeFct(double Shape[15], double dShape[15][3], double m[3]);
extern double Penta15ShapeFct(double Shape[15], double DShape[15][3], double m[3], double X[15], double Y[15], double Z[15]);
extern double computePenta15DShapeFct(double dShape[15][3], double X[15], double Y[15], double Z[15], double (*DShape)[3] = 0);

extern void addBtCBtoK3DSolid(FullSquareMatrix &K, double (*DShape)[3], double C[6][6], double alpha, int nnodes, int* ls);
extern void addNtDNtoM3DSolid(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls, double (*D)[3] = 0);


#ifdef CHECK_JACOBIAN //HB: for checking zero/small & constant sign of jacobian over the el.
extern int checkJacobian(double *J, int *jSign, int elId, const char* mssg= 0, double atol = 0.0, bool stop=true, FILE* file=stderr);
#endif

//HB: for stresses & strains evaluation in case of ansitropic constitutive matrix
extern void computeStressAndEngStrain3DSolid(double Stress[6], double Strain[6], double C[6][6], double (*DShape)[3], double* U, int nnodes, int* ls=0);
extern double computeStress3DSolid(double Stress[6],double Strain[6], double C[6][6]);
extern double computeVonMisesStress(double Stress[6]);
extern double computeVonMisesStrain(double Strain[6]);

Penta15::Penta15(int* nodenums)
{
  for(int i=0; i<15; i++)
    nn[i] = nodenums[i];

  penta15Corotator = 0;

  cFrame = 0; 
  cCoefs = 0;
}

Element *
Penta15::clone()
{
  return(new Penta15(*this));
}

void
Penta15::renum(int *table)
{
 for(int i=0; i<15; i++) { nn[i] = table[nn[i]]; }
}

void
Penta15::getVonMises(Vector& stress,Vector& weight,CoordSet &cs,
                            Vector& elDisp, int strInd, int surface, 
			    double *ndTemps, double ylayer, double zlayer, int avgnum)
{
  fprintf(stderr," *** WARNING: Penta15::getVonMises: NOT implemented. Return null stresses/strains.\n");
  stress.zero();
  weight.zero();
  return;

  //HB: NOT TESTED, VALIDATED YET ...
#ifdef PENTA15_DEBUG 
  fprintf(stderr," *** DEBUG: Get in Penta15::getVonMises.\n");  
#endif

  const int nnodes = 15;
  //const int ndofs  = 45;
  weight = 1.0;
  
  double X[15], Y[15], Z[15];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
     
  // Flags sands17 to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = true; //HB: to force averaging the  Von Mises stress & strain. 
                         //    I don't really know the rational behind that, but its is necessary 
			 //    if we want to recover the same result as the old (fortran based) implementation
  
  double elStress[15][7];
  double elStrain[15][7];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef PENTA15_DEBUG 
    cerr<<" *** DEBUG: Penta15::getVonMises, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[15][3] = {{  0., 0., -1.},{  1.,  0., -1.},{ 0.,  1., -1.},
                                {  0., 0.,  1.},{  1.,  0.,  1.},{ 0.,  1.,  1.},
                                { 0.5, 0., -1.},{ 0.5, 0.5, -1.},{ 0., 0.5, -1.},
                                { 0.5, 0.,  1.},{ 0.5, 0.5,  1.},{ 0., 0.5,  1.},
                                {  0., 0.,  0.},{  1.,  0.,  0.},{ 0.,  1.,  0.}};

  double Shape[15], DShape[15][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta15ShapeFct(Shape, DShape, m, X, Y, Z); 
    computeStressAndEngStrain3DSolid(elStress[inode],elStrain[inode], C, DShape, elDisp.data(), nnodes);

    if(ndTemps) {
      double Tref  = prop->Ta;
      double alpha = prop->W ;
      double eT    = alpha*(ndTemps[inode]-Tref);
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
  
  // compute average Von Mises stress and/or Von Mises strain: to match old Frotran code 
  if(vmflg&meanVms) {
    double vms = 0.0;
    for(int inode=0; inode<nnodes; inode++) vms += elStress[inode][6];
    vms /= nnodes;
    for(int inode=0; inode<nnodes; inode++) elStress[inode][6] = vms;
  }
  if(strainFlg&meanVms) {
    double vms = 0.0;
    for(int inode=0; inode<nnodes; inode++) vms += elStrain[inode][6];
    vms /= nnodes;
    for(int inode=0; inode<nnodes; inode++) elStrain[inode][6] = vms;
  }
  
  // fill the output array stress with the requested stress or strain component 
  for(int inode=0; inode<nnodes; inode++)
    if(strInd < 7) 
      stress[inode] = elStress[inode][strInd];
    else 
      stress[inode] = elStrain[inode][strInd-7];
}

void
Penta15::getAllStress(FullM& stress,Vector& weight,CoordSet &cs,
                      Vector& elDisp, int strInd,int surface, double *ndTemps)
{
  fprintf(stderr," *** WARNING: Penta15::getAllStress: NOT implemented. Return null stress/strain.\n");
  stress.zero();
  return;
  
  //HB: NOT TESTED, VALIDATED YET ...
#ifdef PENTA15_DEBUG 
  fprintf(stderr," *** DEBUG: Get in Brick15::getAllStress.\n");  
#endif

  const int nnodes = 15;
  //const int ndofs  = 45;
  weight = 1.0;
  
  double X[15], Y[15], Z[15];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
      
  double elStress[15][7];
  double elStrain[15][7];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef PENTA15_DEBUG 
    cerr<<" *** DEBUG: Brick15::getVonMises, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[15][3] = {{  0., 0., -1.},{  1.,  0., -1.},{ 0.,  1., -1.},
                                {  0., 0.,  1.},{  1.,  0.,  1.},{ 0.,  1.,  1.},
                                { 0.5, 0., -1.},{ 0.5, 0.5, -1.},{ 0., 0.5, -1.},
                                { 0.5, 0.,  1.},{ 0.5, 0.5,  1.},{ 0., 0.5,  1.},
                                {  0., 0.,  0.},{  1.,  0.,  0.},{ 0.,  1.,  0.}};

  double Shape[15], DShape[15][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta15ShapeFct(Shape, DShape, m, X, Y, Z); 
    computeStressAndEngStrain3DSolid(elStress[inode],elStrain[inode], C, DShape, elDisp.data(), nnodes);

    if(ndTemps) {
      double Tref  = prop->Ta;
      double alpha = prop->W ;
      double eT    = alpha*(ndTemps[inode]-Tref);
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
  
  // Get Element Principals
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};
  for(int j=0; j<6; ++j){ // get average stress/strain  
    for(int i=0; i<nnodes; ++i) 
      svec[j] += stress[i][j];
    svec[j] /= nnodes;
  }
  // Convert Engineering to Tensor Strains
  if(strInd != 0) { svec[3] /= 2; svec[4] /= 2; svec[5] /= 2; }
  pstress(svec,pvec); // compute principal stress (or strain) & direction
  for(int i=0; i<nnodes; ++i) 
    for(int j=0; j<3; ++j) 
      stress[i][j+6] = pvec[j];  
}

double
Penta15::getMass(CoordSet& cs)
{
  fprintf(stderr," *** WARNING: Penta15::getMass: NOT implemented. Return null mass.\n");
  return(0.0);
}

void
Penta15::getGravityForce(CoordSet& cs,double *gravityAcceleration, 
                                Vector& gravityForce, int gravflg, GeomState *geomState)
{
  fprintf(stderr," *** WARNING: Penta15::getGravityForce: NOT implemented. Return null gravity force.\n");
  gravityForce.zero();
}

FullSquareMatrix
Penta15::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
  //int status = 0;
  const int nnodes  = 15;
  const int ndofs   = 45;
  //const int numgauss= 4;

  double X[15], Y[15], Z[15];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  FullSquareMatrix M(ndofs,mel);
  //double *gravityAcceleration = 0, *grvfor = 0;
  //int grvflg = 0;
  //double totmas = 0.0;

  if(cmflg) { //HB: consistent mass matrix
    //fprintf(stderr," *** In Penta15::massMatrix: make consistent mass matrix.\n");
    M.zero();
    int ls[45] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,
                  1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,
                  2,5,8,11,14,17,20,23,15,29,32,35,38,41,44};

    // hard coded order 8 triangle quadrature rule: {r,s,t(=1-r-s),w}
    // This is probably too much for this el.
    double w1 = 0.116786275715379 ;
    double l1 = 0.501415509658179 ;
    double l2 = 0.249286745170910 ;
    double l3 = 1. - l1 - l2;

    double w2 = 0.050844906370207 ;
    double h1 = 0.873821971016996 ;
    double h2 = 0.063089014491502 ;
    double h3 = 1. - h1 - h2;
 
    double w3 = 0.082851075618374 ;
    double k1 = 0.053145049847817 ;
    double k2 = 0.310352451033784 ;
    double k3 = 1. - k1 - k2;
    
    // divide weight by 1/2
    // (unit triangle area = 1/2 and NOT 1)
    w1 *= 0.5; w2 *= 0.5; w3 *= 0.5;

    double TriGPt12[12][4]= {{l1, l2, l3, w1},
                             {l2, l3, l1, w1},
                             {l3, l1, l2, w1},
                             {h1, h2, h3, w2},
                             {h2, h3, h1, w2}, 
                             {h3, h1, h2, w2},
                             {k1, k2, k3, w3},
                             {k1, k3, k2, w3},
                             {k2, k1, k3, w3},
                             {k3, k1, k2, w3},
                             {k2, k3, k1, w3},
                             {k3, k2, k1, w3}};

    // integration: loop over Gauss pts
    double wxy,wz,w;
    int ngpz =  4; // number of (linear) Gauss pts in the (local) z direction
    int ngpxy= 12; // numbder of (triangular) integration pts (in the local x-y plane)
    double m[3], Shape[15], DShape[15][3];
    double dOmega;//det of jacobian
    int jSign = 0;
    for(int iz=1;iz<=ngpz;iz++){ // routine lgauss uses fortran indexing
      _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz); // get z position & weight of the Gauss pt
      for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
        // get x, y  position & weight of the Gauss pt
        m[0] = TriGPt12[ixy][0]; m[1] = TriGPt12[ixy][1]; wxy = TriGPt12[ixy][3];
        dOmega = Penta15ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef PENTA15_DEBUG
        fprintf(stderr," *** In Penta15::massMatrix: iz = %d, ixy = %d, J = %e\n",iz,ixy+1,dOmega);
#endif
#ifdef CHECK_JACOBIAN
        checkJacobian(&dOmega, &jSign, getGlNum()+1, "Penta15::massMatrix");
#endif
        w = fabs(dOmega)*wxy*wz*prop->rho;
        addNtDNtoM3DSolid(M, Shape, w, nnodes, ls);
      }
    }
  } else { // Lumped mass matrix
    fprintf(stderr," *** In Penta15::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return(M);
}

FullSquareMatrix
Penta15::stiffness(CoordSet &cs, double *d, int flg)
{
  double X[15], Y[15], Z[15];
                                                                                                                                         
  //int status = 0;
  const int nnodes  = 15;
  const int ndofs   = 45;
  //const int numgauss= 4;

  cs.getCoordinates(nn, nnodes, X, Y, Z);

  FullSquareMatrix K(ndofs,d);
  K.zero();
  //double *gravityAcceleration = 0, *grvfor = 0;
  //int grvflg = 0;
  //double totmas = 0.0;
  double C[6][6];                                                                                                                                       
  if(cCoefs) {  // PJSA 3-30-05: orthotropic material
    // transform local constitutive matrix to global
    rotateConstitutiveMatrix(cCoefs, cFrame, C); 
  }
  else { // isotropic material
   _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C); 
  }

  int ls[45] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,
                1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,
                2,5,8,11,14,17,20,23,15,29,32,35,38,41,44};

  // hard coded order 4 triangle quadrature rule: {r,s,t(=1-r-s),w}
  // This is probably too much for this el.
  double a1 = 0.445978490915965;
  double b1 = 0.091576213509771;
  double a2 = 1-2.*a1;
  double b2 = 1-2.*b1;
  double w1 = 0.111690797839005;
  double w2 = 0.054975871827661;
  double TriGPt6[6][4] = {{a1, a1, 1.-a1-a1, w1},
                          {a2, a1, 1.-a2-a1, w1},
                          {a1, a2, 1.-a1-a2, w1},
                          {b1, b1, 1.-b1-b1, w2},
                          {b2, b1, 1.-b2-b1, w2},
                          {b1, b2, 1.-b1-b2, w2}};

  // integration: loop over Gauss pts
  double wxy,wz,w;
  int ngpz = 3; // number of (linear) Gauss pts in the (local) z direction
  int ngpxy= 6; // numbder of (triangular) integration pts (in the local x-y plane)
  double m[3], Shape[15], DShape[15][3];
  double dOmega;//det of jacobian
  int jSign = 0;
#ifdef PENTA15_DEBUG
  double Vol = 0.0;
#endif
  for(int iz=1;iz<=ngpz;iz++){ // routine lgauss uses fortran indexing
    _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz); // get z position & weight of the Gauss pt
    for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
      // get x, y  position & weight of the Gauss pt
      m[0] = TriGPt6[ixy][0]; m[1] = TriGPt6[ixy][1]; wxy = TriGPt6[ixy][3];
      dOmega = Penta15ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef PENTA15_DEBUG
      fprintf(stderr," *** In Penta15::stiffness: iz = %d, ixy = %d, J = %e\n",iz,ixy+1,dOmega);
      Vol += dOmega*wxy*wz;
#endif
#ifdef CHECK_JACOBIAN
      checkJacobian(&dOmega, &jSign, getGlNum()+1, "Penta::stiffness");
#endif
      w = fabs(dOmega)*wxy*wz;
      addBtCBtoK3DSolid(K, DShape, C, w, nnodes, ls);
    }
  }
#ifdef PENTA15_DEBUG
  fprintf(stderr," *** In Penta15::stiffness: volume = %e\n",Vol);
#endif
  return(K);
}

int
Penta15::numNodes() { 
  if(useFull)
    return(15); 
  else
    return(6);
}

int
Penta15::numDofs() { return(45); }

int*
Penta15::nodes(int *p)
{
  if(useFull)
    {
      if(!p) p = new int[15];
      for(int i=0; i<15; i+=3) {
	p[i] = nn[i]; p[i+1] = nn[i+1]; p[i+2] = nn[i+2];
      }
      return(p);
    }
  else
    {
      if(!p) p = new int[6];
      for(int i=0; i<6; i+=1) {
	p[i] = nn[i];
      }
      return(p);
    }
}

int*
Penta15::dofs(DofSetArray &dsa, int *p)
{
  if(!p) p = new int[45];
  for(int i=0; i<15; i++)
    dsa.number(nn[i],DofSet::XYZdisp, p+3*i);

  return(p);
}

void
Penta15::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

int
Penta15::getTopNumber() { return(197); } 

int
Penta15::numTopNodes() { return(15); } 

//-------------------------------------------------------------------
void
Penta15::getThermalForce(CoordSet &cs, Vector &ndTemps,
                                Vector &elementThermalForce, int glflag, 
				GeomState *geomState)
{
  fprintf(stderr," *** WARNING: Penta15::getThermalForce: NOT implemented. Return NULL force.\n");
  elementThermalForce.zero();
  return;
 
  //HB: NOT TESTED, VALIDATED YET ...
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes= 15;	
  const int ndofs = 45;

  // extract nodes coordinates
  double X[15], Y[15], Z[15];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // initialize nodal thermal forces
  for(int i=0; i<ndofs; i++) elementThermalForce[i] = 0.0;

  // get material props & constitutive matrix
  double Tref  = prop->Ta;
  double alpha = prop->W ;
  //double coef  = prop->E*(1.-2.*prop->nu);
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
    //cerr<<" *** DEBUG: in Pentahedral::stiffness, anisotropic/orthotropic material\n";
     rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  // Integate over the element: F = Int[Bt.ThermaStress]
  // with ThermalStress = C.ThermalStrain, with ThermalStrain = alpha.theta.[1, 1, 1, 0, 0, 0]'
  // where theta = T(M)-Tref = Sum[inode][N[inode]*(ndTemps[inode] - Tref)]
  // N[inode] is the shape fct at node inode
  // M is the position in the real frame, m its associated position in the reference
  // element frame
  // !!! USE BRUTE FORCE: NUMERICAL INETGRATION BY GAUSS PTS !!!
  double Shape[15], DShape[15][3], m[3];
  double wxy,wz,w,J;
  int ngpz = 3; // number of (linear) Gauss pts in the (local) z direction
  int ngpxy= 6; // numbder of (triangular) integration pts (in the local x-y plane)

  // hard coded order 4 triangle quadrature rule: {r,s,t(=1-r-s),w}
  double a1 = 0.445978490915965;
  double b1 = 0.091576213509771;
  double a2 = 1-2.*a1;
  double b2 = 1-2.*b1;
  double w1 = 0.111690797839005;
  double w2 = 0.054975871827661;
  double TriGPt6[6][4] = {{a1, a1, 1.-a1-a1, w1},
                          {a2, a1, 1.-a2-a1, w1},
                          {a1, a2, 1.-a1-a2, w1},
                          {b1, b1, 1.-b1-b1, w2},
                          {b2, b1, 1.-b2-b1, w2},
                          {b1, b2, 1.-b1-b2, w2}};
 
  if (geomState)  { // GEOMETRICAL NONLINEAR ANALYSIS
     fprintf(stderr," *** ERROR: Penta15::getThermalForce NOT supported for geometric nonlinear analysis. Abort.\n");
     exit(-1);
  } else { // GEOMETRICAL LINEAR ANALYSIS
    // integration: loop over Gauss pts
    for(int iz=1;iz<=ngpz;iz++){        // z Gauss pts
      // get z position & weight of the Gauss pt
      _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz);
      for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
        // get x, y  position & weight of the Gauss pt
        m[0] = TriGPt6[ixy][0]; m[1] = TriGPt6[ixy][1]; wxy = TriGPt6[ixy][3];
        // compute shape fcts & their derivatives at the Gauss pt
        J = Penta15ShapeFct(Shape, DShape, m, X, Y, Z);
        w = wxy*wz*fabs(J);
        // compute thermal stresses
        double eT = 0.0;
        for(int inode=0; inode<nnodes; inode++) eT += alpha*Shape[inode]*(ndTemps[inode] - Tref);
        double thermalStrain[6] = {eT,eT,eT,0.0,0.0,0.0};
	double thermalStress[6] = {0.0,0.0,0.0,0.0,0.0,0.0}; 
        computeStress3DSolid(thermalStress, thermalStrain, C); // thermalStress <- C.thermalStrain
       // sum contribution
        for(int inode=0; inode<nnodes; inode++){
          elementThermalForce[3*inode  ] += w*(DShape[inode][0]*thermalStress[0] + DShape[inode][1]*thermalStress[3] + DShape[inode][2]*thermalStress[5]);
          elementThermalForce[3*inode+1] += w*(DShape[inode][0]*thermalStress[3] + DShape[inode][1]*thermalStress[1] + DShape[inode][2]*thermalStress[4]);
          elementThermalForce[3*inode+2] += w*(DShape[inode][0]*thermalStress[5] + DShape[inode][1]*thermalStress[4] + DShape[inode][2]*thermalStress[2]);
        }  
      }
    }
  }
  
  //cerr<<" ---------------------------------------"<<endl;
  //for(int i=0; i<78; i++){
  //  cerr<<" elementThermalForce["<<i<<"] = "<<elementThermalForce[i]<<endl;
  //}
}

//---------------------------------------------------------------------------------
Corotator *
Penta15::getCorotator(CoordSet &cs, double *kel, int , int )
{
 fprintf(stderr," *** WARNING: Penta15::getCorotator: NOT implemented. Abort.\n");
 exit(-1);
 return((Corotator*)0); 
}

void Penta15::buildCorotator(CoordSet &cs)  
{
  fprintf(stderr," *** WARNING: Penta15::buildCorotator: NOT implemented. Abort.\n");
  exit(-1);
  //if (!brickCorotator)
  //  brickCorotator = new PentaCorotator(nn, prop->E, prop->nu, cs);
}

