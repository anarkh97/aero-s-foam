// ---------------------------------------------------------------------
// HB - 05-22-05
// ---------------------------------------------------------------------
// 32 nodes brick element 
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
// EXPERIMENTAL ...
// ---------------------------------------------------------------------
#include <cstdio>
#include <iostream>
#include <Element.d/Brick32.d/Brick32.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Corotational.d/BrickCorotator.h>

#define CHECK_JACOBIAN //HB: force check nullity & constant sign of jacobian over el.
//#define BRICK32_DEBUG 

extern bool useFull;

extern "C" {
void  _FORTRAN(brkcmt)(double&, double&, double*);

void _FORTRAN(lgauss)(const int &, int &, double *, double *);
}

//HB: for anisotropic elastic constitutive matrix
void rotateConstitutiveMatrix(double *Cin, double *T33, double Cout[6][6]);

// HB: for stiffness matrix with ansitropic constitutive matrix and/or consistent mass matrix
extern void   Hexa32ShapeFct(double Shape[32], double dShape[32][3], double m[3]);
extern double Hexa32ShapeFct(double Shape[32], double DShape[32][3], double m[3], double X[32], double Y[32], double Z[32]);
extern double computeHexa32DShapeFct(double dShape[32][3], double X[32], double Y[32], double Z[32], double (*DShape)[3] = 0);

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
 
Brick32::Brick32(int* nodenums)
{
  for(int i=0; i<32; i++)
    nn[i] = nodenums[i];

  brick32Corotator = 0;

  cFrame = 0; 
  cCoefs = 0;
  mat = 0;
}

Element *
Brick32::clone()
{
  return(new Brick32(*this));
}

void
Brick32::renum(int *table)
{
 for(int i=0; i<32; i++) { nn[i] = table[nn[i]]; }
}

void
Brick32::renum(EleRenumMap& table)
{
 for(int i=0; i<32; i++) { nn[i] = table[nn[i]]; }
}

// Stress evaluation in case of isotropic & anisotropic elastic constitutive matrix 
void
Brick32::getVonMises(Vector& stress,Vector& weight,CoordSet &cs,
                     Vector& elDisp, int strInd, int surface, 
		     double *ndTemps, double ylayer, double zlayer, int avgnum)
{
  fprintf(stderr," *** WARNING: Brick32::getVonMises: NOT implemented. Return null stresses/strains.\n");
  stress.zero();
  weight.zero();
  return;

  //HB: NOT TESTED, VALIDATED YET ...
#ifdef BRICK32_DEBUG 
  fprintf(stderr," *** DEBUG: Get in Brick32::getVonMises.\n");  
#endif

  const int nnodes = 32;
  //const int ndofs  = 96;
  weight = 1.0;
  
  double X[32], Y[32], Z[32];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
     
  // Flags sands17 to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = false; //HB: to force averaging the  Von Mises stress & strain. 
                         //    I don't really know the rational behind that, but its is necessary 
			 //    if we want to recover the same result as the old (fortran based) implementation
  
  double elStress[32][7];
  double elStrain[32][7];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef BRICK32_DEBUG 
    cerr<<" *** DEBUG: Brick32::getVonMises, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[32][3] = {{-1.   ,-1.,-1.},{ 1.   ,-1.,-1.},{ 1., 1.   ,-1.},{-1., 1.   ,-1.},
                                {-1.   ,-1., 1.},{ 1.   ,-1., 1.},{ 1., 1.   , 1.},{-1., 1.   , 1.},
                                {-1./3.,-1.,-1.},{ 1./3.,-1.,-1.},{ 1.,-1./3.,-1.},{ 1., 1./3.,-1.},
		                { 1./3., 1.,-1.},{-1./3., 1.,-1.},{-1., 1./3.,-1.},{-1.,-1./3.,-1.},
                                {-1./3.,-1., 1.},{ 1./3.,-1., 1.},{ 1.,-1./3., 1.},{ 1., 1./3., 1.},
                                { 1./3., 1., 1.},{-1./3., 1., 1.},{-1., 1./3., 1.},{-1.,-1./3., 1.},
                                {-1.,-1.,-1./3.},{-1.,-1., 1./3.},{ 1.,-1.,-1./3.},{ 1.,-1., 1./3.},
                                { 1., 1.,-1./3.},{ 1., 1., 1./3.},{-1., 1.,-1./3.},{-1., 1., 1./3.}};

  double Shape[32], DShape[32][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Hexa32ShapeFct(Shape, DShape, m, X, Y, Z); 
    computeStressAndEngStrain3DSolid(elStress[inode], elStrain[inode], C, DShape, elDisp.data(), nnodes);

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

// Stress evaluation in case of isotropic & anisotropic elastic constitutive matrix 
void
Brick32::getAllStress(FullM& stress,Vector& weight,CoordSet &cs,
                      Vector& elDisp, int strInd,int surface, double *ndTemps)
{
  fprintf(stderr," *** WARNING: Brick32::getAllStress: NOT implemented. Return null stress/strain.\n");
  stress.zero();
  return;

  //HB: NOT TESTED, VALIDATED YET ...
#ifdef BRICK32_DEBUG 
  fprintf(stderr," *** DEBUG: Get in Brick32::getAllStress.\n");  
#endif

  const int nnodes = 32;
  //const int ndofs  = 96;
  weight = 1.0;
  
  double X[32], Y[32], Z[32];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
      
  double elStress[32][7];
  double elStrain[32][7];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef BRICK32_DEBUG 
    cerr<<" *** DEBUG: Brick32::getVonMises, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[32][3] = {{-1.   ,-1.,-1.},{ 1.   ,-1.,-1.},{ 1., 1.   ,-1.},{-1., 1.   ,-1.},
                                {-1.   ,-1., 1.},{ 1.   ,-1., 1.},{ 1., 1.   , 1.},{-1., 1.   , 1.},
                                {-1./3.,-1.,-1.},{ 1./3.,-1.,-1.},{ 1.,-1./3.,-1.},{ 1., 1./3.,-1.},
		                { 1./3., 1.,-1.},{-1./3., 1.,-1.},{-1., 1./3.,-1.},{-1.,-1./3.,-1.},
                                {-1./3.,-1., 1.},{ 1./3.,-1., 1.},{ 1.,-1./3., 1.},{ 1., 1./3., 1.},
                                { 1./3., 1., 1.},{-1./3., 1., 1.},{-1., 1./3., 1.},{-1.,-1./3., 1.},
                                {-1.,-1.,-1./3.},{-1.,-1., 1./3.},{ 1.,-1.,-1./3.},{ 1.,-1., 1./3.},
                                { 1., 1.,-1./3.},{ 1., 1., 1./3.},{-1., 1.,-1./3.},{-1., 1., 1./3.}};

  double Shape[32], DShape[32][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Hexa32ShapeFct(Shape, DShape, m, X, Y, Z); 
    computeStressAndEngStrain3DSolid(elStress[inode], elStrain[inode], C, DShape, elDisp.data(), nnodes);

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
  
  // Get Element Principals for each node without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};

  for (int i=0; i<nnodes; ++i) {
    for (int j=0; j<6; ++j) {
      svec[j] = stress[i][j];
    }
// Convert Engineering to Tensor Strains
    if(strInd != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec,pvec);
    for (int j=0; j<3; ++j) {
      stress[i][j+6] = pvec[j];
    }
  }  
}

double
Brick32::getMass(CoordSet& cs)
{
  const int nnodes = 32;
  double *X = new double[nnodes], *Y = new double[nnodes], *Z = new double[nnodes];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  const int numgauss= 4;
  double wx, wy, wz;
  double m[3], Shape[32], DShape[32][3];
  double dOmega; // det of jacobian
  double v = 0.0; // volume
  for(int i = 1; i <= numgauss; i++) {
    _FORTRAN(lgauss)(numgauss, i, &m[0], &wx);
    for(int j = 1; j <= numgauss; j++) {
      _FORTRAN(lgauss)(numgauss, j, &m[1], &wy);
      for(int k = 1; k <= numgauss; k++) {
        _FORTRAN(lgauss)(numgauss, k, &m[2], &wz);
        dOmega = Hexa32ShapeFct(Shape, DShape, m, X, Y, Z);
        v += fabs(dOmega)*wx*wy*wz;
      }
    }
  }

  delete [] X; delete [] Y; delete [] Z;

  return v*prop->rho;
}

void
Brick32::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                         Vector& gravityForce, int gravflg, GeomState *geomState)
{
  int nnodes = 32;

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
    double x[32], y[32], z[32];
    cs.getCoordinates(nn, numNodes(), x, y, z);

    double lforce[32];
    for(int i=0; i<nnodes; ++i) lforce[i] = 0.0;

    int numgauss = 4;
    // integration: loop over Gauss pts
    double wx, wy, wz, w;
    double m[3], Shape[32], DShape[32][3];
    double dOmega; // det of jacobian
    for(int i = 1; i <= numgauss; i++) {
      _FORTRAN(lgauss)(numgauss, i, &m[0], &wx);
      for(int j = 1; j <= numgauss; j++) {
        _FORTRAN(lgauss)(numgauss, j, &m[1], &wy);
        for(int k = 1; k <= numgauss; k++) {
          _FORTRAN(lgauss)(numgauss, k, &m[2], &wz);

          dOmega = Hexa32ShapeFct(Shape, DShape, m, x, y, z);
          w = fabs(dOmega)*wx*wy*wz*prop->rho;

          for (int n = 0; n < nnodes; ++n)
            lforce[n] += w*Shape[n];
        }
      }
    }

    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = lforce[i]*gravityAcceleration[0];
      gravityForce[3*i+1] = lforce[i]*gravityAcceleration[1];
      gravityForce[3*i+2] = lforce[i]*gravityAcceleration[2];
    }
  }
}

FullSquareMatrix
Brick32::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
  //int status = 0;
  const int nnodes  = 32;
  const int ndofs   = 96;
  const int numgauss= 4;

  double X[32], Y[32], Z[32];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  FullSquareMatrix M(ndofs,mel);

  if(cmflg) { //HB: consistent mass matrix
    //fprintf(stderr," *** In Brick32::massMatrix: make consistent mass matrix.\n");
    M.zero();
    int ls[96] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,
                  1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,79,82,85,88,91,94,
                  2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95};
                                                                                                                                       
    // integration: loop over Gauss pts
    double wx,wy,wz,w;
    double m[3], Shape[32], DShape[32][3];
    double dOmega;//det of jacobian
    int jSign = 0;                                                                                                                                        
    for(int i=1;i<=numgauss;i++){
      _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
      for(int j=1;j<=numgauss;j++){
        _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
        for(int k=1;k<=numgauss;k++){
          _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
          dOmega = Hexa32ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef BRICK32_DEBUG
          fprintf(stderr," *** In Brick32::massMatrix: i = %d, j = %d, k = %d, J = %e\n",i,j,k,dOmega);
#endif
#ifdef CHECK_JACOBIAN
          checkJacobian(&dOmega, &jSign, getGlNum()+1, "Brick32::massMatrix");
#endif
          w = fabs(dOmega)*wx*wy*wz*prop->rho;
          addNtDNtoM3DSolid(M, Shape, w, nnodes, ls);
        }
      }
    }
  } else { // Lumped mass matrix
    fprintf(stderr," *** In Brick32::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return(M);
}

FullSquareMatrix
Brick32::stiffness(CoordSet &cs, double *d, int flg)
{
  //int status = 0;
  const int nnodes  = 32;
  const int ndofs   = 96;
  const int numgauss= 4;

  double X[32], Y[32], Z[32];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

#ifdef BRICK32_DEBUG
  fprintf(stderr," Brick32 nodes = \n");
  for(int i=0; i<32; i++)
    fprintf(stderr,"  node[%2d] = %6d = %e  %e  %e \n",i+1, nn[i], X[i], Y[i], Z[i]);
  fprintf(stderr,"\n");
#endif

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

  int ls[96] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,
                1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,79,82,85,88,91,94,
                2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95};
                                                                                                                                         
  // integration: loop over Gauss pts
  double wx,wy,wz,w;
  double m[3], Shape[32], DShape[32][3];
  double dOmega;//det of jacobian
  int jSign = 0;                                                                                                                                          
  for(int i=1;i<=numgauss;i++){
    _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
    for(int j=1;j<=numgauss;j++){
      _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
      for(int k=1;k<=numgauss;k++){
        _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
        dOmega = Hexa32ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef BRICK32_DEBUG
        fprintf(stderr," *** In Brick32::stiffness: i = %d, j = %d, k = %d, J = %e\n",i,j,k,dOmega);
#endif
#ifdef CHECK_JACOBIAN
        checkJacobian(&dOmega, &jSign, getGlNum()+1, "Brick32::stiffness");
#endif
        w = fabs(dOmega)*wx*wy*wz;
        addBtCBtoK3DSolid(K, DShape, C, w, nnodes, ls);
      }
    }
  }

  return(K);
}

int
Brick32::numNodes() { 
  if(useFull)
    return(32); 
  else
    return(8);   // to ignore effect of mid-size nodes in dec
}

int
Brick32::numDofs() { return(96); }

int*
Brick32::nodes(int *p)
{
  if(useFull)
    {
      if(!p) p = new int[32];
      for(int i=0; i<32; i+=2) {
	p[i] = nn[i]; p[i+1] = nn[i+1];
      }
      return(p);
    }
  else
    {
      if(!p) p = new int[8];
      for(int i=0; i<8; i+=1) {
	p[i] = nn[i]; 
      }
      return(p);
    }
}

int*
Brick32::dofs(DofSetArray &dsa, int *p)
{
  if(!p) p = new int[96];
  for(int i=0; i<32; i++)
    dsa.number(nn[i],DofSet::XYZdisp, p+3*i);

  return(p);
}

void
Brick32::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

int
Brick32::getTopNumber() { return(191); } 

int Brick32::numTopNodes() { return(32); } 


void
Brick32::getThermalForce(CoordSet &cs, Vector &ndTemps,
                         Vector &elementThermalForce, int glflag, 
			 GeomState *geomState)
{
  fprintf(stderr," *** WARNING: Brick32::getThermalForce NOT implemented. Return NULL force.\n");
  elementThermalForce.zero();
  return;
 
  //HB: NOT TESTED, VALIDATED YET ...
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes= 32;
  const int ndofs = 60;
  const int numgauss= 4;

  // extract nodes coordinates
  double X[32], Y[32], Z[32];                                                                                                                             
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
#ifdef BRICK32_DEBUG 
    cerr<<" *** DEBUG: Brick32::getThermalForc, anisotropic/orthotropic material\n";
#endif
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
  double Shape[32], DShape[32][3], m[3];
  double wx,wy,wz,w,J; 
  if (geomState)  {
    fprintf(stderr," *** ERROR: Brick32::getThermalForce NOT supported for geometric nonlinear analysis. Abort.\n");
    exit(-1);
  }
  else {
    // integration: loop over Gauss pts
    for(int i=1;i<=numgauss;i++){
      _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
      for(int j=1;j<=numgauss;j++){
        _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
        for(int k=1;k<=numgauss;k++){
          _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
          J = Hexa32ShapeFct(Shape, DShape, m, X, Y, Z);
          w = fabs(J)*wx*wy*wz;
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
  }
  //cerr<<" -------------------------"<<endl;
  //for(int i=0; i<96; i++) 
   //cerr<<" elementThermalForce["<<i+1<<"] = "<<elementThermalForce[i]<<endl;  
}

#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NLHexahedral.h>
#include <Corotational.d/MatNLCorotator.h>

void
Brick32::setMaterial(NLMaterial *_mat)
{
  mat = _mat;
}

int
Brick32::numStates()
{
  int numGaussPoints = 64;
  return (mat) ? numGaussPoints*mat->getNumStates(): 0;
}

Corotator*
Brick32::getCorotator(CoordSet &cs, double *kel, int , int )
{
#ifdef USE_EIGEN3
  if(!mat)
    mat = new StVenantKirchhoffMat(prop->rho, prop->E, prop->nu);
  MatNLElement *ele = new NLHexahedral32(nn);
  ele->setMaterial(mat);
  ele->setGlNum(glNum);
  return new MatNLCorotator(ele);
#else
  printf("WARNING: Corotator not implemented for element %d\n", glNum+1); return 0;
#endif
} 
