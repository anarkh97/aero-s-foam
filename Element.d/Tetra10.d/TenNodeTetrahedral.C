#include	<cstdio>

#include	<Element.d/Tetra10.d/TenNodeTetrahedral.h>
#include        <Math.d/matrix.h>
#include	<Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>
#include        <Utils.d/pstress.h>

#include <Driver.d/PolygonSet.h>

#define USE_NEW_TET10_STIFF //HB: force using new implementation of the stiffness matrix
                            //    that deals with ansitropic constitutive matrix
#define CHECK_JACOBIAN //HB: force check nullity & constant sign of jacobian over el.
#define TETRA10_DEBUG

extern "C"      {
void	_FORTRAN(mstf25)(double*, double*, double*, double&, double&, double*, 
		const int&, double*, double*, int &);

void    _FORTRAN(sands25)(const int&,    double*,    double*,    double*,    
                             double&,    double&,    double*,    double*,    
                             double*, const int&, const int&, const int&, 
                          const int&, const int&, const int&);

#ifdef USE_NEW_TET10_STIFF //HB: used in the new version of the stiffness matrix
void  _FORTRAN(brkcmt)(double&, double&, double*);
#endif
}

// HB: for stiffness matrix with ansitropic constitutive matrix and/or consistent mass matrix
extern void rotateConstitutiveMatrix(double *Cin, double *T33, double Cout[6][6]);
extern void Tetra10ShapeFct(double Shape[10], double dShape[10][3], double m[3]);
extern double Tetra10ShapeFct(double Shape[10], double DShape[10][3], double m[3], double X[10], double Y[10], double Z[10]);
extern double computeTet10DShapeFct(double dShape[10][3], double X[10], double Y[10], double Z[10], double (*DShape)[3] = 0);
extern void addBtCBtoK3DSolid(FullSquareMatrix &K, double (*DShape)[3], double C[6][6], double alpha, int nnodes, int* ls);
void addNtDNtoM3DSolid(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls, double (*D)[3] = 0);
                                                                                                                             
#ifdef CHECK_JACOBIAN //HB: for checking zero/small & constant sign of jacobian over the el.
extern int checkJacobian(double *J, int *jSign, int elId, const char* mssg= 0, double atol = 0.0, bool stop=true, FILE* file=stderr);
#endif

//HB: for stresses & strains evaluation in case of ansitropic constitutive matrix
extern void computeStressAndEngStrain3DSolid(double Stress[6], double Strain[6], double C[6][6], double (*DShape)[3], double* U, int nnodes, int* ls=0);
extern double computeStress3DSolid(double Stress[6],double Strain[6], double C[6][6]);
extern double computeVonMisesStress(double Stress[6]);
extern double computeVonMisesStrain(double Strain[6]);
                                                                                                                             

typedef double Coord[3];
extern void compute_coupling(Coord *p, double A[30][10]);

extern bool useFull;

TenNodeTetrahedral::TenNodeTetrahedral(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nn[2] = nodenums[2];
  nn[3] = nodenums[3];
  nn[4] = nodenums[4];
  nn[5] = nodenums[5];
  nn[6] = nodenums[6];
  nn[7] = nodenums[7];
  nn[8] = nodenums[8];
  nn[9] = nodenums[9];

  cFrame = 0;
  cCoefs = 0;
  mat = 0;
}

Element *
TenNodeTetrahedral::clone()
{
 return new TenNodeTetrahedral(*this);
}

void
TenNodeTetrahedral::renum(int *table)
{
     nn[0] = table[nn[0]];
     nn[1] = table[nn[1]];
     nn[2] = table[nn[2]];
     nn[3] = table[nn[3]];
     nn[4] = table[nn[4]];
     nn[5] = table[nn[5]];
     nn[6] = table[nn[6]];
     nn[7] = table[nn[7]];
     nn[8] = table[nn[8]];
     nn[9] = table[nn[9]];
}

void
TenNodeTetrahedral::getVonMises(Vector& stress,Vector& weight,CoordSet &cs,
                                Vector& elDisp, int strInd,int surface,double* ndTemps, double ylayer, double zlayer, int avgnum)
{
 if(cCoefs){
   getVonMisesAniso(stress, weight, cs,
 		    elDisp, strInd, surface, ndTemps,
 		    ylayer, zlayer);
   return;		  
 }

  if(stress.size() != 10) fprintf(stderr,"ERROR in stress vector.\n");
  if(weight.size() != 10) fprintf(stderr,"ERROR in weight vector.\n");

  weight = 1.0;

  double x[10], y[10], z[10]; 
  cs.getCoordinates(nn, numNodes(), x, y, z);

  int vmflg = 0;
  if(strInd == 6) vmflg  = 1; // Flags sands25 to calculate Von Mises stress.

  const int maxgus = 15; // maximum gauss points 
  const int maxstr = 7; 

  double elStress[maxgus][maxstr];
  double elStrain[maxgus][maxstr];

  int strainFlg = 0;

  const int elm    = 1;
  const int msize  = 1;
  const int outerr = 6;

  double nu = prop->nu;
  double E  = prop->E;

  _FORTRAN(sands25)(elm,   x,                 y,                 z,         E,
                     nu, elDisp.data(),(double*)elStress, (double*)elStrain,
                 maxgus,
                 maxstr,msize,            outerr,             vmflg, strainFlg);
  
   double Tref  = prop->Ta;
   double alpha = prop->W;

   int i;
   if((1.0 - 2.0*nu) == 0.0) 
     fprintf(stderr,"Division by Zero! in Tetra10\n");
   double thermalStress[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
   if(strInd < 7) {
     if(strInd == 0 || strInd == 1 || strInd == 2) {
       for(i=0; i<10; ++i)
         thermalStress[i] = E*alpha*(ndTemps[i] - Tref)/(1.0 - 2.0*nu);
     }
     for(i=0; i<10; ++i)
       stress[i] = elStress[i][strInd] - thermalStress[i];
   } else {
     for(i=0; i<10; ++i)
       stress[i] = elStrain[i][strInd-7];
   }

   //fprintf(stderr,"stress*stress = %e\n",stress*stress);
   fflush(stderr);
}

void
TenNodeTetrahedral::getAllStress(FullM& stress,Vector& weight,CoordSet &cs,
                                 Vector& elDisp, int strInd,int surface,double* ndTemps)
{
  if(cCoefs) {
    getAllStressAniso(stress, weight, cs,
 	 	      elDisp, strInd, surface, ndTemps);
    return;		     
  }		  
  weight = 1.0;

  double x[10], y[10], z[10]; 
  cs.getCoordinates(nn, numNodes(), x, y, z);

  int vmflg = 0;

  const int maxgus = 15; // maximum gauss points 
  const int maxstr = 7; 

  double elStress[maxgus][maxstr];
  double elStrain[maxgus][maxstr];

  int strainFlg = 0;

  const int elm    = 1;
  const int msize  = 1;
  const int outerr = 6;

  double nu = prop->nu;
  double E  = prop->E;

  _FORTRAN(sands25)(elm,   x,                 y,                 z,         E,
                     nu, elDisp.data(),(double*)elStress, (double*)elStrain,
                 maxgus,
                 maxstr,msize,            outerr,             vmflg, strainFlg);
  
   double Tref  = prop->Ta;
   double alpha = prop->W;

// Store all Stress or all Strain as defined by strInd
        int i,j;
        if(strInd == 0) {
          if((1.0 - 2.0*nu) == 0.0) {
            fprintf(stderr,"Division by Zero! in Tetra10\n");
          }

          double thermalStress[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	  if(ndTemps) 
            for(i=0; i<10; ++i) 
              thermalStress[i] = E*alpha*(ndTemps[i] - Tref)/(1.0 - 2.0*nu);
          
          for (i=0; i<10; ++i) {
            for (j=0; j<3; ++j) {
              stress[i][j] = elStress[i][j] - thermalStress[i];
            }
            for (j=3; j<6; ++j) {
              stress[i][j] = elStress[i][j];
            }
          }
        }
        else {
          for (i=0; i<10; ++i) {
            for (j=0; j<6; ++j) {
              stress[i][j] = elStrain[i][j];
            }
          }
        }

// Get Element Principals for each node without averaging
        double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
        double pvec[3] = {0.0,0.0,0.0};

        for (i=0; i<10; ++i) {
          for (j=0; j<6; ++j) {
            svec[j] = stress[i][j];
          }
// Convert Engineering to Tensor Strains
          if(strInd != 0) {
            svec[3] /= 2;
            svec[4] /= 2;
            svec[5] /= 2;
          }
          pstress(svec,pvec);
          for (j=0; j<3; ++j) {
            stress[i][j+6] = pvec[j];
          }
        }
}

double
TenNodeTetrahedral::getMass(CoordSet& cs)
{
  double x[10], y[10], z[10]; 
  cs.getCoordinates(nn, numNodes(), x, y, z);

  const int numgauss = 15;
  extern double dp[15][10][3]; // contains the values of the Tet10 shape fct at the 15 integration pts

  // integration: loop over Gauss pts
  // reuse the 15 pts integration rule (order 5) -> reuse arrays vp1 & dp 
  const double weight[15] = {1.975308731198311E-02, 1.198951396316977E-02,
                             1.198951396316977E-02, 1.198951396316977E-02,
                             1.198951396316977E-02, 1.151136787104540E-02,
                             1.151136787104540E-02, 1.151136787104540E-02,
                             1.151136787104540E-02, 8.818342350423336E-03,
                             8.818342350423336E-03, 8.818342350423336E-03,
                             8.818342350423336E-03, 8.818342350423336E-03,
                             8.818342350423336E-03};
  double dOmega; // det of jacobian
  double volume = 0.0;
  for(int i = 0; i < numgauss; i++) {
    dOmega = computeTet10DShapeFct(dp[i], x, y, z);
    volume += fabs(dOmega)*weight[i];
  }

  return volume*prop->rho;
}

void
TenNodeTetrahedral::getGravityForce(CoordSet& cs, double *gravityAcceleration,  
                                    Vector& gravityForce, int gravflg, GeomState *geomState)
{
  int nnodes = 10;

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

    double x[10], y[10], z[10];
    cs.getCoordinates(nn, numNodes(), x, y, z);

    double lforce[10];
    for(int i = 0; i < nnodes; ++i) lforce[i] = 0.0;

    const int numgauss = 15;
    extern double dp[15][10][3]; // arrays vp1 & dp contain the values of the Tet10 shape fct  & their 
    extern double vp1[15][10];   // derivatives w.r.t reference coordinate system at the 15 integration pts

    // integration: loop over Gauss pts
    // reuse the 15 pts integration rule (order 5) -> reuse arrays vp1 & dp 
    const double weight[15] = {1.975308731198311E-02, 1.198951396316977E-02,
                               1.198951396316977E-02, 1.198951396316977E-02,
                               1.198951396316977E-02, 1.151136787104540E-02,
                               1.151136787104540E-02, 1.151136787104540E-02,
                               1.151136787104540E-02, 8.818342350423336E-03,
                               8.818342350423336E-03, 8.818342350423336E-03,
                               8.818342350423336E-03, 8.818342350423336E-03,
                               8.818342350423336E-03};
    double dOmega; //det of jacobian
    for(int i = 0; i < numgauss; i++) {
      dOmega = computeTet10DShapeFct(dp[i], x, y, z);
      double w = fabs(dOmega)*weight[i]*prop->rho;

      for(int n = 0; n < nnodes; ++n)
        lforce[n] += w*vp1[i][n];
    }

    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = lforce[i]*gravityAcceleration[0];
      gravityForce[3*i+1] = lforce[i]*gravityAcceleration[1];
      gravityForce[3*i+2] = lforce[i]*gravityAcceleration[2];
    }
  }
}

// ---------------------------------------------------------------------------------
//HB (06-11-05): add anisotropic/orthotropic material case
void
TenNodeTetrahedral::getThermalForce(CoordSet &cs,Vector &ndTemps,
                                    Vector &elementThermalForce, 
				    int glflag, GeomState *geomState)
{
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes= 10;                                                                                           
  const int ndofs = 30;

  double X[10], Y[10], Z[10];
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
 #ifdef TETRA10_DEBUG 
   fprintf(stderr," *** DEBUG: in TenNodeTetrahedral::getThermalForce, anisotropic/orthotropic material\n");
 #endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);

  const int numgauss = 15;     // use 15 pts integration rule (order 5)
  extern double dp[15][10][3]; // arrays vp1 & dp contain the values of the Tet10 shape fct  & their 
  extern double vp1[15][10];   // derivatives w.r.t reference coordinate system at the 15 integration pts
 
  if(geomState) { // GEOMETRICAL NONLINEAR ANALYSIS
    fprintf(stderr," *** ERROR: TenNodeTetrahedral::getThermalForce: NOT supported for Non-Linear Geometric analysis. Abort.\n");
    exit(-1);
  } else { // GEOMETRICAL LINEAR ANALYSIS
    // BRUTE FORCE: NUMERICAL INETGRATION BY GAUSS PTS 
    // integration: loop over Gauss pts
    // reuse the 15 pts integration rule (order 5) -> reuse arrays vp1 & dp 
    const double weight[15] = {1.975308731198311E-02, 1.198951396316977E-02,
                               1.198951396316977E-02, 1.198951396316977E-02,
                               1.198951396316977E-02, 1.151136787104540E-02,
                               1.151136787104540E-02, 1.151136787104540E-02,
                               1.151136787104540E-02, 8.818342350423336E-03,
                               8.818342350423336E-03, 8.818342350423336E-03,
                               8.818342350423336E-03, 8.818342350423336E-03,
                               8.818342350423336E-03};
    double DShape[10][3];
    double w, J;
    int jSign = 0;
    for(int i=0;i<numgauss;i++){
      J = computeTet10DShapeFct(dp[i],X,Y,Z,DShape);
      double* Shape = &vp1[i][0];
#ifdef CHECK_JACOBIAN
      checkJacobian(&J, &jSign, getGlNum()+1, "TenNodeTetrahedral::getThermalForce");
#endif
      w = 0.5*fabs(J)*weight[i]; //HB: need to check why the 0.5 is needed to recover the old implementation ...
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
  
   /* OLD CODE (isotropic material)
   double elC[30][10];

   compute_coupling(p,elC);

   double alpha = prop->W;
   double E     = prop->E;
   double nu    = prop->nu;
   double Tref  = prop->Ta;

   double coef = (E*alpha/(1.0-2.0*nu));
   int i,j;
   for(i=0; i<30; ++i) {
     elementThermalForce[i] = 0.0;
     for(j=0; j<10; ++j)
       elementThermalForce[i] += coef*elC[i][j]*(ndTemps[j] - Tref);
   }
  */
  }
//  for (int kk = 0; kk < 30; kk++)
//    fprintf(stderr, "eForce[%d] = %e\n", kk, elementThermalForce[kk]);
}

FullSquareMatrix
TenNodeTetrahedral::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
  double X[10], Y[10], Z[10];
  cs.getCoordinates(nn, numNodes(), X, Y, Z);

  const int numgauss = 15;     // use 15 pts integration rule (order 5)
  extern double dp[15][10][3]; // arrays vp1 & dp contain the values of the Tet10 shape fct  & their 
  extern double vp1[15][10];   // derivatives w.r.t reference coordinate system at the 15 integration pts

  const int nnodes = 10;
  const int ndofs = 30;
  FullSquareMatrix M(ndofs,mel);
                                                                                                                             
  if(cmflg) { //HB: consistent mass matrix
    //fprintf(stderr," *** In TenNodeTetrahedral::massMatrix: make consistent mass matrix.\n");
    M.zero();
    int ls[30] = {0, 3, 6, 9,12,15,18,21,24,27,
                  1, 4, 7,10,13,16,19,22,25,28,
                  2, 5, 8,11,14,17,20,23,26,29};
                                                                                                                             
    // integration: loop over Gauss pts
    // reuse the 15 pts integration rule (order 5) -> reuse arrays vp1 & dp 
    const double weight[15] = {1.975308731198311E-02, 1.198951396316977E-02,
                               1.198951396316977E-02, 1.198951396316977E-02,
                               1.198951396316977E-02, 1.151136787104540E-02,
                               1.151136787104540E-02, 1.151136787104540E-02,
                               1.151136787104540E-02, 8.818342350423336E-03,
                               8.818342350423336E-03, 8.818342350423336E-03,
                               8.818342350423336E-03, 8.818342350423336E-03,
                               8.818342350423336E-03};
    double dOmega;//det of jacobian
    int jSign = 0;
    for(int i = 0; i < numgauss; i++){
      dOmega = computeTet10DShapeFct(dp[i], X, Y, Z);
#ifdef CHECK_JACOBIAN
      checkJacobian(&dOmega, &jSign, getGlNum()+1, "TenNodeTetrahedral::massMatrix");
#endif
      double w = fabs(dOmega)*weight[i]*prop->rho;
      addNtDNtoM3DSolid(M, vp1[i], w, nnodes, ls);
    }
  } 
  else { // Lumped mass matrix
    fprintf(stderr," *** In TenNodeTetrahedral::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }
                                                                                                                             
  return(M);
}

#ifndef USE_NEW_TET10_STIFF
// HB (06/19/05): replaced by new implementation (see below) to deal with anisotropic constitutive matrix
FullSquareMatrix
TenNodeTetrahedral::stiffness(CoordSet &cs, double *d, int flg)
{

  double x[10], y[10], z[10];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  extern double dp[15][10][3];
  extern double vp1[15][10];

  int status;
  _FORTRAN(mstf25)(x, y, z, prop->E, prop->nu, (double*) d, 
                   numDofs(),(double *) dp, (double *)vp1,
                   status);

   if(status != 0) {
          fprintf(stderr, "Error in 10 node tetrahedron:");
	  int i;
          for(i = 0; i < 10; ++i)
             fprintf(stderr, " %d",nn[i]+1);
          fprintf(stderr, "\n");
   }

  FullSquareMatrix ret(30,d);
	    
  return ret;
}
#else
//---------------------------------------------------------------------------------
//HB (06/19/05)  new implementation of the Tetra4 stiffness matrix to deal
//               with anisotropic constitutive matrix
FullSquareMatrix
TenNodeTetrahedral::stiffness(CoordSet &cs,double *d, int flg)
{
  //fprintf(stderr, " *** WARNING: use new Tetrahedral::stiffness method (HB)\n");
  double X[10], Y[10], Z[10];
  cs.getCoordinates(nn, numNodes(), X, Y, Z);

  const int numgauss = 15;    // use 15 pts integration rule (order 5) (order 2 is exact for stiffness if linear mapping)
  extern double dp[15][10][3]; // arrays vp1 & dp contain the values of the Tet10 shape fct  & their
  //extern double vp1[15][10];   // derivatives w.r.t reference coordinate system at the 15 integration pts
                                                                                                                                        
  const int nnodes= 10;
  const int ndofs = 30;
  int status = 0;
  FullSquareMatrix K(ndofs,d);
  K.zero();
                                                                                                                                        
  // get constitutive matrix
  double C[6][6];
  if(cCoefs) {  //anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef TET10_DEBUG    
    cerr<<" *** DEBUG: in TenNodeTetrahedral::stiffness, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);

  int ls[30] = {0, 3, 6, 9,12,15,18,21,24,27,
                1, 4, 7,10,13,16,19,22,25,28,
                2, 5, 8,11,14,17,20,23,26,29};
                                                                                                                                        
  // integration: loop over Gauss pts
  // reuse the 15 pts integration rule (order 5) -> reuse arrays vp1 & dp
  const double weight[15] = {1.975308731198311E-02, 1.198951396316977E-02,
                             1.198951396316977E-02, 1.198951396316977E-02,
                             1.198951396316977E-02, 1.151136787104540E-02,
                             1.151136787104540E-02, 1.151136787104540E-02,
                             1.151136787104540E-02, 8.818342350423336E-03,
                             8.818342350423336E-03, 8.818342350423336E-03,
                             8.818342350423336E-03, 8.818342350423336E-03,
                             8.818342350423336E-03};
  double DShape[10][3];
  double dOmega;//det of jacobian
  int jSign = 0;
  for(int i=0;i<numgauss;i++){
    dOmega = computeTet10DShapeFct(dp[i],X,Y,Z,DShape);
#ifdef CHECK_JACOBIAN
    checkJacobian(&dOmega, &jSign, getGlNum()+1, "TenNodeTetrahedral::stiffness");
#endif
    double w = fabs(dOmega)*weight[i];
    addBtCBtoK3DSolid(K, DShape, C, w, nnodes, ls);
  } 

  if(status != 0) {
    fprintf(stderr, " *** FATAL ERROR in TenNodeTetrahedral::stiffness");
    fprintf(stderr, " -> Tetra10's nodes: \n");
    for(int i = 0; i < nnodes; ++i)
       fprintf(stderr, " %d",nn[i]+1);
    fprintf(stderr, "\n");
    exit(-1);
  }
                                                                                                                                        
  return(K);
}
#endif

// WARNING : BEFORE DEC integration into FEM this guy reported always 10
int
TenNodeTetrahedral::numNodes()
{
  if(useFull) return(10);
  return 4;
}

int*
TenNodeTetrahedral::nodes(int *p)
{
  if(useFull)
    {

      if(p == 0) p = new int[10];
      p[0] = nn[0];
      p[1] = nn[1];
      p[2] = nn[2];
      p[3] = nn[3];
      p[4] = nn[4];
      p[5] = nn[5];
      p[6] = nn[6];
      p[7] = nn[7];
      p[8] = nn[8];
      p[9] = nn[9];
      
    }
  else
    {
      if(p == 0) p = new int[4];
      p[0] = nn[0];
      p[1] = nn[1];
      p[2] = nn[2];
      p[3] = nn[3];
    }
  
  return p;
}

int
TenNodeTetrahedral::numDofs()
{
  return 30;
}

int*
TenNodeTetrahedral::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[30];

  dsa.number(nn[0],DofSet::XYZdisp,p);
  dsa.number(nn[1],DofSet::XYZdisp,p+3);
  dsa.number(nn[2],DofSet::XYZdisp,p+6);
  dsa.number(nn[3],DofSet::XYZdisp,p+9);
  dsa.number(nn[4],DofSet::XYZdisp,p+12);
  dsa.number(nn[5],DofSet::XYZdisp,p+15);
  dsa.number(nn[6],DofSet::XYZdisp,p+18);
  dsa.number(nn[7],DofSet::XYZdisp,p+21);
  dsa.number(nn[8],DofSet::XYZdisp,p+24);
  dsa.number(nn[9],DofSet::XYZdisp,p+27);

  return p;
}

void
TenNodeTetrahedral::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

int
TenNodeTetrahedral::getTopNumber()
{
  return 125;//return 11;
}

/*int
TenNodeTetrahedral::facelist(PolygonSet &ps,int *list)
{
 if(list != 0) {
   list[0] = ps.addTri6(nn[0],nn[2],nn[1],nn[6],nn[5],nn[4]);
   list[1] = ps.addTri6(nn[0],nn[3],nn[2],nn[7],nn[9],nn[6]);
   list[2] = ps.addTri6(nn[0],nn[1],nn[3],nn[4],nn[8],nn[7]);
   list[3] = ps.addTri6(nn[1],nn[2],nn[3],nn[5],nn[9],nn[8]);
 }
 return 4;
}
*/

//---------------------------------------------------------------------------------
//HB (06/11/05)
// Stress evaluation in case of anisotropic elastic constitutive matrix 
void
TenNodeTetrahedral::getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs,
			             Vector &elDisp, int strInd, int surface, double *ndTemps,
			             double ylayer, double zlayer, int avgnum)
{
#ifdef TETRA10_DEBUG 
  fprintf(stderr," *** DEBUG: Get in TenNodeTetrahedral::getVonMisesAniso.\n");  
#endif

  const int nnodes = 10;
  //const int ndofs  = 30;
  weight = 1.0;
  
  double X[10], Y[10], Z[10];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
     
  // Flags sands17 to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = false; //HB: to force averaging the  Von Mises stress & strain. 
                         //    I don't really know the rational behind that, but its is necessary 
			 //    if we want to recover the same result as the old (fortran based) implementation
  
  double elStress[10][7];
  double elStrain[10][7];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef TETRA10_DEBUG 
    cerr<<" *** DEBUG: TenNodeTetrahedral::getVonMisesAniso, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[10][3] = {{0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0},
                                {0.5,0.0,0.0},{0.5,0.5,0.0},{0.0,0.5,0.0},
				{0.0,0.0,0.5},{0.5,0.0,0.5},{0.0,0.5,0.5}};

  double Shape[10], DShape[10][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Tetra10ShapeFct(Shape, DShape, m, X, Y, Z); 
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


//---------------------------------------------------------------------------------
//HB (06/11/05)
// Stress evaluation in case of anisotropic elastic constitutive matrix 
void
TenNodeTetrahedral::getAllStressAniso(FullM& stress,Vector& weight,CoordSet &cs,
                                      Vector& elDisp, int strInd,int surface ,double* ndTemps)
{
#ifdef TETRA10_DEBUG 
  fprintf(stderr," *** DEBUG: Get in TenNodeTetrahedral::getAllStressAniso.\n");  
#endif

  const int nnodes = 10;
  //const int ndofs  = 30;
  weight = 1.0;
  
  double X[10], Y[10], Z[10];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
      
  double elStress[10][7];
  double elStrain[10][7];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef TETRA10_DEBUG 
    cerr<<" *** DEBUG: TenNodeTetrahedral::getVonMisesAniso, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[10][3] = {{0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0},
                                {0.5,0.0,0.0},{0.5,0.5,0.0},{0.0,0.5,0.0},
				{0.0,0.0,0.5},{0.5,0.0,0.5},{0.0,0.5,0.5}};

  double Shape[10], DShape[10][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Tetra10ShapeFct(Shape, DShape, m, X, Y, Z); 
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
  
  // Get Element Principals without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};
  for(int i=0; i<nnodes; ++i){
    for(int j=0; j<6; ++j) 
      svec[j] = stress[i][j];
 
    // Convert Engineering to Tensor Strains
    if(strInd != 0) { svec[3] /= 2; svec[4] /= 2; svec[5] /= 2; }
    pstress(svec,pvec); // compute principal stress (or strain) & direction
    for(int j=0; j<3; ++j) 
      stress[i][j+6] = pvec[j];
  }
}

#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NLTetrahedral.h>
#include <Corotational.d/MatNLCorotator.h>

void
TenNodeTetrahedral::setMaterial(NLMaterial *_mat)
{
  mat = _mat;
}

int
TenNodeTetrahedral::numStates()
{
  int numGaussPoints = 15;
  return (mat) ? numGaussPoints*mat->getNumStates(): 0;
}

Corotator*
TenNodeTetrahedral::getCorotator(CoordSet &cs, double *kel, int , int )
{
  if(!mat)
    mat = new StVenantKirchhoffMat(prop->rho, prop->E, prop->nu);
  MatNLElement *ele = new NLTetrahedral10(nn);
  ele->setMaterial(mat);
  ele->setGlNum(glNum);
  return new MatNLCorotator(ele);
}
