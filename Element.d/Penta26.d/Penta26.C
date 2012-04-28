// ---------------------------------------------------------------------
// HB - 05-22-05
// ---------------------------------------------------------------------
// 26 nodes wedge element 
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
// EXPERIMENTAL ...
// ---------------------------------------------------------------------
#include <cstdio>
#include <iostream>
#include <Element.d/Penta26.d/Penta26.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>

#define CHECK_JACOBIAN //HB: force check nullity & constant sign of jacobian over el.
//#define PENTA26_DEBUG

extern bool useFull;

extern "C" {
void  _FORTRAN(brkcmt)(double&, double&, double*);
}

//HB: for anisotropic elastic constitutive matrix
void rotateConstitutiveMatrix(double *Cin, double *T33, double Cout[6][6]);

// HB: for stiffness matrix with ansitropic constitutive matrix and/or consistent mass matrix
extern void   Penta26ShapeFct(double Shape[26], double dShape[26][3], double m[3]);
extern double Penta26ShapeFct(double Shape[26], double DShape[26][3], double m[3], double X[26], double Y[26], double Z[26]);
extern double computePenta26DShapeFct(double dShape[26][3], double X[26], double Y[26], double Z[26], double (*DShape)[3] = 0);

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

double gauss3d9[18][3] = {
     { 0.166666666666667, 0.166666666666667, -0.774596669241483 },
     { 0.166666666666667, 0.666666666666667, -0.774596669241483 },
     { 0.666666666666667, 0.166666666666667, -0.774596669241483 },
     { 0.000000000000000, 0.500000000000000, -0.774596669241483 },
     { 0.500000000000000, 0.000000000000000, -0.774596669241483 },
     { 0.500000000000000, 0.500000000000000, -0.774596669241483 },
     { 0.166666666666667, 0.166666666666667, 0.0 },
     { 0.166666666666667, 0.666666666666667, 0.0 },
     { 0.666666666666667, 0.166666666666667, 0.0 },
     { 0.000000000000000, 0.500000000000000, 0.0 },
     { 0.500000000000000, 0.000000000000000, 0.0 },
     { 0.500000000000000, 0.500000000000000, 0.0 },
     { 0.166666666666667, 0.166666666666667, 0.774596669241483 },
     { 0.166666666666667, 0.666666666666667, 0.774596669241483 },
     { 0.666666666666667, 0.166666666666667, 0.774596669241483 },
     { 0.000000000000000, 0.500000000000000, 0.774596669241483 },
     { 0.500000000000000, 0.000000000000000, 0.774596669241483 },
     { 0.500000000000000, 0.500000000000000, 0.774596669241483 }
  };

double weight3d9[18] = { 
     0.083333333333333, 0.083333333333333, 0.083333333333333,
     0.009259259259259, 0.009259259259259, 0.009259259259259,
     0.133333333333333, 0.133333333333333, 0.133333333333333,
     0.014814814814815, 0.014814814814815, 0.014814814814815,
     0.083333333333333, 0.083333333333333, 0.083333333333333,
     0.009259259259259, 0.009259259259259, 0.009259259259259
  };

                                                                                                                                         
Penta26::Penta26(int* nodenums)
{
  for(int i=0; i<26; i++)
    nn[i] = nodenums[i];

  penta26Corotator = 0;

  cFrame = 0; 
  cCoefs = 0;
  mat = 0;
}

Element *
Penta26::clone()
{
  return(new Penta26(*this));
}

void
Penta26::renum(int *table)
{
 for(int i=0; i<26; i++) { nn[i] = table[nn[i]]; }
}

void
Penta26::getVonMises(Vector& stress,Vector& weight,CoordSet &cs,
                            Vector& elDisp, int strInd, int surface, 
			    double *ndTemps, double ylayer, double zlayer, int avgnum)
{
  fprintf(stderr," *** WARNING: Penta26::getVonMises: NOT implemented. Return null stresses/strains.\n");
  stress.zero();
  weight.zero();
  return;

  //HB: NOT TESTED, VALIDATED YET ...
#ifdef PENTA26_DEBUG 
  fprintf(stderr," *** DEBUG: Get in Penta26::getVonMises.\n");  
#endif

  const int nnodes = 26;
  //const int ndofs  = 78;
  weight = 1.0;
  
  double X[26], Y[26], Z[26];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
     
  // Flags sands17 to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = true; //HB: to force averaging the  Von Mises stress & strain. 
                         //    I don't really know the rational behind that, but its is necessary 
			 //    if we want to recover the same result as the old (fortran based) implementation
  
  double elStress[26][7];
  double elStrain[26][7];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef PENTA26_DEBUG 
    cerr<<" *** DEBUG: Penta26::getVonMises, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[26][3] = {{  0. ,  0. ,  -1.},{  1. ,  0. ,  -1.},{  0. ,  1. ,  -1.},
                                {  0. ,  0. ,   1.},{  1. ,  0. ,   1.},{  0. ,  1. ,   1.},
                                {1./3.,  0. ,  -1.},{2./3.,  0. ,  -1.},{2./3.,1./3.,  -1.},
                                {1./3.,2./3.,  -1.},{  0. ,2./3.,  -1.},{  0. ,1./3.,  -1.},
                                {1./3.,  0. ,   1.},{2./3.,  0. ,   1.},{2./3.,1./3.,   1.},
                                {1./3.,2./3.,   1.},{  0. ,2./3.,   1.},{  0. ,1./3.,   1.},
                                {  0. ,  0. ,-1./3},{  0. ,  0. , 1./3},{  1. ,  0. ,-1./3},
                                {  1. ,  0. , 1./3},{  0. ,  1. ,-1./3},{  0. ,  1. , 1./3},
                                {1./3., 1./3. ,-1.},{1./3., 1./3. , 1.}};

  double Shape[26], DShape[26][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta26ShapeFct(Shape, DShape, m, X, Y, Z); 
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

void
Penta26::getAllStress(FullM& stress,Vector& weight,CoordSet &cs,
                      Vector& elDisp, int strInd,int surface, double *ndTemps)
{
  fprintf(stderr," *** WARNING: Penta26::getAllStress: NOT implemented. Return null stress/strain.\n");
  stress.zero();
  return;
  
  //HB: NOT TESTED, VALIDATED YET ...
#ifdef PENTA26_DEBUG 
  fprintf(stderr," *** DEBUG: Get in Brick26::getAllStress.\n");  
#endif

  const int nnodes = 26;
  //const int ndofs  = 78;
  weight = 1.0;
  
  double X[26], Y[26], Z[26];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
      
  double elStress[26][7];
  double elStrain[26][7];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef PENTA26_DEBUG 
    cerr<<" *** DEBUG: Brick26::getVonMises, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[26][3] = {{  0. ,  0. ,  -1.},{  1. ,  0. ,  -1.},{  0. ,  1. ,  -1.},
                                {  0. ,  0. ,   1.},{  1. ,  0. ,   1.},{  0. ,  1. ,   1.},
                                {1./3.,  0. ,  -1.},{2./3.,  0. ,  -1.},{2./3.,1./3.,  -1.},
                                {1./3.,2./3.,  -1.},{  0. ,2./3.,  -1.},{  0. ,1./3.,  -1.},
                                {1./3.,  0. ,   1.},{2./3.,  0. ,   1.},{2./3.,1./3.,   1.},
                                {1./3.,2./3.,   1.},{  0. ,2./3.,   1.},{  0. ,1./3.,   1.},
                                {  0. ,  0. ,-1./3},{  0. ,  0. , 1./3},{  1. ,  0. ,-1./3},
                                {  1. ,  0. , 1./3},{  0. ,  1. ,-1./3},{  0. ,  1. , 1./3},
                                {1./3., 1./3. ,-1.},{1./3., 1./3. , 1.}};

  double Shape[26], DShape[26][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta26ShapeFct(Shape, DShape, m, X, Y, Z); 
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
Penta26::getMass(CoordSet& cs)
{
  const int nnodes = 26;
  const int ngauss = 18;
  double *x = new double[nnodes];
  double *y = new double[nnodes];
  double *z = new double[nnodes];
  double *shape = new double[nnodes];
  double (*dShape)[3] = new double[nnodes][3];

  // get global coordinates of the nodes
  cs.getCoordinates(nn, nnodes, x, y, z);

  // integration: loop over gauss pts
  double volume = 0.0;
  for(int i = 0; i < ngauss; i++){
    double dOmega = Penta26ShapeFct(shape, dShape, gauss3d9[i], x, y, z);
    volume += fabs(dOmega)*weight3d9[i];
  }

  delete [] x; delete [] y; delete [] z;
  delete [] shape; delete [] dShape;

  return volume*prop->rho;

}

void
Penta26::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                         Vector& gravityForce, int gravflg, GeomState *geomState)
{
  int nnodes = 26;

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

    const int ngauss = 18;
    double *x = new double[nnodes];
    double *y = new double[nnodes];
    double *z = new double[nnodes];
    double *shape = new double[nnodes];
    double (*dShape)[3] = new double[nnodes][3];
    double *lforce = new double[nnodes];

    // get global coordinates of the nodes
    cs.getCoordinates(nn, nnodes, x, y, z);

    for(int i = 0; i < nnodes; ++i) lforce[i] = 0.0;

    // integration: loop over gauss pts
    for(int i = 0; i < ngauss; i++){
      double dOmega = Penta26ShapeFct(shape, dShape, gauss3d9[i], x, y, z);
      for(int j = 0; j < nnodes; ++j) lforce[j] += fabs(dOmega)*weight3d9[i]*shape[j];
    }

    for(int i = 0; i < nnodes; ++i)
      for(int j = 0; j < 3; ++j)
        gravityForce[3*i+j] = lforce[i]*gravityAcceleration[j]*prop->rho;

    delete [] x; delete [] y; delete [] z;
    delete [] shape; delete [] dShape;
    delete [] lforce;
  }
}

FullSquareMatrix
Penta26::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
  const int ndofs = 78;
  FullSquareMatrix M(ndofs, mel);

  if(cmflg) { // consistent mass matrix
    const int nnodes = 26;
    const int ngauss = 18;
    double *x = new double[nnodes];
    double *y = new double[nnodes];
    double *z = new double[nnodes];
    double *shape = new double[nnodes];
    double (*dShape)[3] = new double[nnodes][3];
    int ls[78] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,
                  1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,
                  2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77};

    // get global coordinates of the nodes
    cs.getCoordinates(nn, nnodes, x, y, z);

    M.zero();

    // integration: loop over gauss pts
    for(int i = 0; i < ngauss; i++){
      double dOmega = Penta26ShapeFct(shape, dShape, gauss3d9[i], x, y, z);
      double w = fabs(dOmega)*weight3d9[i]*prop->rho;
      addNtDNtoM3DSolid(M, shape, w, nnodes, ls);
    }

    delete [] x; delete [] y; delete [] z;
    delete [] shape; delete [] dShape;
  }
  else { // lumped mass matrix
    fprintf(stderr," *** In Penta26::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return M;
}

FullSquareMatrix
Penta26::stiffness(CoordSet &cs, double *d, int flg)
{
  const int nnodes = 26;
  const int ngauss = 18;
  const int ndofs = 78;

  double *x = new double[nnodes];
  double *y = new double[nnodes];
  double *z = new double[nnodes];
  double *shape = new double[nnodes];
  double (*dShape)[3] = new double[nnodes][3];
  FullSquareMatrix K(ndofs, d);
  int ls[78] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,
                1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,
                2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77};

  double C[6][6];
  if(cCoefs)
    rotateConstitutiveMatrix(cCoefs, cFrame, C);  // orthotropic material
  else
   _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C); // isotropic material

  // get global coordinates of the nodes
  cs.getCoordinates(nn, nnodes, x, y, z);

  K.zero();

  // integration: loop over gauss pts
  for(int i = 0; i < ngauss; i++){
    double J = Penta26ShapeFct(shape, dShape, gauss3d9[i], x, y, z);
    double w = fabs(J)*weight3d9[i];
    addBtCBtoK3DSolid(K, dShape, C, w, nnodes, ls);
  }

  delete [] x; delete [] y; delete [] z;
  delete [] shape; delete [] dShape;

  return K;
}

int
Penta26::numNodes() { 
  if(useFull)
    return(26); 
  else
    return(8);
}

int
Penta26::numDofs() { return(78); }

int*
Penta26::nodes(int *p)
{
  if(useFull)
    {
      if(!p) p = new int[26];
      for(int i=0; i<26; i+=2) {
	p[i] = nn[i]; p[i+1] = nn[i+1];
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
Penta26::dofs(DofSetArray &dsa, int *p)
{
  if(!p) p = new int[78];
  for(int i=0; i<26; i++)
    dsa.number(nn[i],DofSet::XYZdisp, p+3*i);

  return(p);
}

void
Penta26::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

// Treat as Penta6
//int Penta26::getTopNumber() { return(124); }  
//int Penta26::numTopNodes() { return(6); }  
// Treat as Penta26
int Penta26::getTopNumber() { return(192); }  
int Penta26::numTopNodes() { return(26); }  

//-------------------------------------------------------------------
void
Penta26::getThermalForce(CoordSet &cs, Vector &ndTemps,
                                Vector &elementThermalForce, int glflag, 
				GeomState *geomState)
{
  //HB: NOT TESTED, VALIDATED YET ...
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes= 26;	
  const int ndofs = 78;
                                                                                                                   
  // extract nodes coordinates
  double X[26], Y[26], Z[26];                                                                                                                              
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
  double Shape[26], DShape[26][3];
  double w,J;
  int ngauss = 18;

  if (geomState)  { // GEOMETRICAL NONLINEAR ANALYSIS
     fprintf(stderr," *** ERROR: Penta26::getThermalForce NOT supported for geometric nonlinear analysis. Abort.\n");
     exit(-1);
  }
  else { // GEOMETRICAL LINEAR ANALYSIS
    // integration: loop over Gauss pts
    for(int i = 0; i < ngauss; i++) {
      // compute shape fcts & their derivatives at the Gauss pt
      J = Penta26ShapeFct(Shape, DShape, gauss3d9[i], X, Y, Z);
      w = weight3d9[i]*fabs(J);
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

#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NLPentahedral.h>
#include <Corotational.d/MatNLCorotator.h>

void
Penta26::setMaterial(NLMaterial *_mat)
{
  mat = _mat;
}

int
Penta26::numStates()
{
  int numGaussPoints = 18;
  return (mat) ? numGaussPoints*mat->getNumStates(): 0;
}

Corotator*
Penta26::getCorotator(CoordSet &cs, double *kel, int , int )
{
  if(!mat)
    mat = new StVenantKirchhoffMat(prop->rho, prop->E, prop->nu);
  MatNLElement *ele = new NLPentahedral26(nn);
  ele->setMaterial(mat);
  ele->setGlNum(glNum);
  return new MatNLCorotator(ele);
}

