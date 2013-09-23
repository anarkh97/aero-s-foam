#include <cstdio>
#include <iostream>
#include <cmath>

#include <Element.d/Tetra.d/Tetrahedral.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Math.d/FullSquareMatrix.h>
#include <Corotational.d/TetCorotator.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/NLTetrahedral.h>
#include <Corotational.d/MatNLCorotator.h>


#define USE_NEW_TET4_STIFF  //HB: force using new implementation of the stiffness matrix 
                            //    that deals with ansitropic constitutive matrix
#define CHECK_JACOBIAN      //HB: force check nullity & constant sign of jacobian over el.
//#define TETRA4_DEBUG
extern "C"      {

void	_FORTRAN(mstf23)(double*, double*, double*, double&, double&, double*, 
			 const int&, int &);

void	_FORTRAN(mass23)(double*, double*, double*, double&, double*, 
			 const int&, double*, double*, const int&, 
			 double&, const int&);

void    _FORTRAN(sands23)(const int&, double*, double*, double*, double&, 
			  double&, double*, double*, double*, const int&, 
                          const int&,
			  const int&, const int&, const int&, const int&);

void  _FORTRAN(brkcmt)(double&, double&, double*);
}

typedef double Coord[3];
extern void compute_coupling(Coord *p, double A[12][4]);

// HB: for stiffness matrix with ansitropic constitutive matrix and/or consistent mass matrix
extern void rotateConstitutiveMatrix(double *Cin, double *T33, double Cout[6][6]);
extern double Tetra4ShapeFct(double Shape[4], double DShape[4][3], double m[3], double X[4], double Y[4], double Z[4]);
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
                                                                                               

Tetrahedral::Tetrahedral(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nn[2] = nodenums[2];
  nn[3] = nodenums[3];

  tetCorotator = 0;

  cFrame = 0;
  cCoefs = 0;
  mat = 0;
}

Element *
Tetrahedral::clone()
{
 return new Tetrahedral(*this);
}

void
Tetrahedral::renum(int *table)
{
     nn[0] = table[nn[0]];
     nn[1] = table[nn[1]];
     nn[2] = table[nn[2]];
     nn[3] = table[nn[3]];
}

void
Tetrahedral::renum(EleRenumMap& table)
{
     nn[0] = table[nn[0]];
     nn[1] = table[nn[1]];
     nn[2] = table[nn[2]];
     nn[3] = table[nn[3]];
}

void
Tetrahedral::getVonMises(Vector& stress,Vector& weight,CoordSet &cs,
                         Vector& elDisp, int strInd, int surface , double* ndTemps,
			 double ylayer, double zlayer, int avgnum)
{
        if(cCoefs){
          getVonMisesAniso(stress, weight, cs,
			   elDisp, strInd, surface, ndTemps,
			   ylayer, zlayer, avgnum);
	  return;		 
        }
 	weight = 1.0;

	Node &nd1 = cs.getNode(nn[0]);
	Node &nd2 = cs.getNode(nn[1]);
	Node &nd3 = cs.getNode(nn[2]);
	Node &nd4 = cs.getNode(nn[3]);

  	double x[4], y[4], z[4]; 

  	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

        // Flags sands23 to calculate Von Mises stress.
        int vmflg = 0;
        if(strInd == 6) vmflg  = 1;
	int strainFlg = 0;
	if(strInd == 13) strainFlg = 1;

        int maxgus = 4; // maximum gauss points 
        int maxstr = 7; 
        int elm    = 1;
	int outerr = 6;

        double elStress[4][7];
        double elStrain[4][7];

        const int one = 1;

        double Tref  = prop->Ta; // ambient temperature
        double alpha = prop->W;  // thermal expansion coefficient
        double E     = prop->E;  // Young's modulus
        double nu    = prop->nu; // Poisson's ratio

        _FORTRAN(sands23)(elm, x, y, z, E, nu, elDisp.data(), 
          (double*)elStress, (double*)elStrain, maxgus, maxstr,
          one, outerr, vmflg, strainFlg); 

	if(strInd < 7) {
          double thermalStress[4] = {0.0,0.0,0.0,0.0};
          if(strInd == 0 || strInd == 1 || strInd == 2) {
            double coef = (E*alpha)/(1.0 - 2.0*nu);
            thermalStress[0] = coef*(ndTemps[0]-Tref);
            thermalStress[1] = coef*(ndTemps[1]-Tref);
            thermalStress[2] = coef*(ndTemps[2]-Tref);
            thermalStress[3] = coef*(ndTemps[3]-Tref);
          }
          stress[0] = elStress[0][strInd] - thermalStress[0];
          stress[1] = elStress[1][strInd] - thermalStress[1];
          stress[2] = elStress[2][strInd] - thermalStress[2];
          stress[3] = elStress[3][strInd] - thermalStress[3];
	} else if(strInd < 14) {
          stress[0] = elStrain[0][strInd-7];
          stress[1] = elStrain[1][strInd-7];
          stress[2] = elStrain[2][strInd-7];
          stress[3] = elStrain[3][strInd-7];
	}
        else {
          stress[0] = 0;
          stress[1] = 0;
          stress[2] = 0;
          stress[3] = 0;
        } 
}

void
Tetrahedral::getAllStress(FullM& stress,Vector& weight,CoordSet &cs,
                         Vector& elDisp, int strInd,int surface, double* ndTemps)
{
        if(cCoefs) {
	  getAllStressAniso(stress, weight, cs,
			    elDisp, strInd, surface, ndTemps);
	  return;		    
	} 		 
        weight = 1.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);
        Node &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], z[4];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

        int vmflg,strainFlg;
        vmflg  = 0;
        strainFlg = 0;

        int maxgus = 4; // maximum gauss points
        int maxstr = 7;
        int elm    = 1;
        int outerr = 6;

        double elStress[4][7];
        double elStrain[4][7];

        const int one = 1;

        double Tref  = prop->Ta; // ambient temperature
        double alpha = prop->W;  // thermal expansion coefficient
        double E     = prop->E;  // Young's modulus
        double nu    = prop->nu; // Poisson's ratio

        _FORTRAN(sands23)(elm, x, y, z, E, nu, elDisp.data(),
          (double*)elStress, (double*)elStrain, maxgus, maxstr,
          one, outerr, vmflg, strainFlg);

// Store all Stress or all Strain as defined by strInd
        int i,j;
        if(strInd == 0) {
          double thermalStress[4] = {0.0,0.0,0.0,0.0};
	  if(ndTemps) {
            double coef = (E*alpha)/(1.0 - 2.0*nu);
            thermalStress[0] = coef*(ndTemps[0]-Tref);
            thermalStress[1] = coef*(ndTemps[1]-Tref);
            thermalStress[2] = coef*(ndTemps[2]-Tref);
            thermalStress[3] = coef*(ndTemps[3]-Tref);
	  }
          for (i=0; i<4; ++i) {
            for (j=0; j<3; ++j) {
              stress[i][j] = elStress[i][j] - thermalStress[i];
            }
            for (j=3; j<6; ++j) {
              stress[i][j] = elStress[i][j];
            }
          }
        }
        else {
          for (i=0; i<4; ++i) {
            for (j=0; j<6; ++j) {
              stress[i][j] = elStrain[i][j];
            }
          }
        }

// Get Element Principals for each node without averaging
        double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
        double pvec[3] = {0.0,0.0,0.0};

        for (i=0; i<4; ++i) {
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
Tetrahedral::getMass(CoordSet& cs)
{

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);

  double x[4], y[4], z[4], ElementMassMatrix[12][12];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

  double *gravityAcceleration = 0, *grvfor = 0, totmas = 0.0;

  int grvflg = 0, masflg = 1;

  const int numdof = 12;
  _FORTRAN(mass23)(x, y, z, prop->rho, (double*)ElementMassMatrix, numdof,
                   gravityAcceleration, grvfor, grvflg, totmas, masflg);

  return totmas;
}

void
Tetrahedral::getThermalForce(CoordSet &cs, Vector &ndTemps, 
                             Vector &elementThermalForce, 
			     int glflag,GeomState *geomState)
{
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes= 4;                                                                                           
  const int ndofs = 12;

  double X[4], Y[4], Z[4];
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
 #ifdef TETRA4_DEBUG 
   fprintf(stderr," *** DEBUG: in Tetrahedral::getThermalForce, anisotropic/orthotropic material\n");
 #endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  if(geomState) { // GEOMETRICAL NONLINEAR ANALYSIS
    if(cCoefs)
      fprintf(stderr," *** WARNING: Tetrahedral::getThermalForce: evaluation of thermal forces NOT supported for anisotropic constitutive matrix. Use isotropic coeffs.\n");
    // F = E*alpha/(1-2*nu) : de/dq_i  SUM(xi_j)*dT_j  * 0.25*Vol
    double coef = 0.5*alpha*prop->E*(1.-2.*prop->nu);
    double nGrad[4][3];
    double dedU[12][6];
    double dOmega = tetCorotator->computeShapeGrad(cs, nGrad)/24.0;
    tetCorotator->computeStrainGrad(*geomState, nGrad, dedU);

    double deltaT = 0.0;
    for (int i = 0; i < nnodes; i++)
      deltaT += ndTemps[i]-Tref;

    for (int i = 0; i < ndofs; ++i)
      elementThermalForce[i] = dOmega*coef * (dedU[i][0]+dedU[i][1]+dedU[i][2])*deltaT;

  } else { // GEOMETRICAL LINEAR ANALYSIS
    // BRUTE FORCE: NUMERICAL INETGRATION BY GAUSS PTS 
    // hard coded order 1 tetrahedral quadrature rule: {r,s,t,u(=1-r-s-t),w}
    double TetGPt1[1][5] = {{1./4.,1./4.,1./4.,1./4.,1./6.}};
    
    // integration: loop over Gauss pts
    int ngp = 1; // number of integration pts
    double m[3], Shape[4], DShape[4][3];
    double w, J;//det of jacobian
    int jSign = 0;                                                                                                                                
    for(int igp=0;igp<ngp;igp++){
      // get x, y, z  position & weight of the integration pt
      m[0] = TetGPt1[igp][0]; m[1] = TetGPt1[igp][1];  m[2] = TetGPt1[igp][2]; w = TetGPt1[igp][4];
      J = Tetra4ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
      checkJacobian(&J, &jSign, getGlNum()+1,"Tetrahedral::getThermalForce");
#endif
      w *= 0.5*fabs(J); //HB: need to check why the 0.5 is needed to recover the old implementation ...
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
    double elC[12][4];
    compute_coupling(p, elC);

    for(int i=0; i<12; ++i) {
      elementThermalForce[i] = 0.0;
      for(int j=0; j<4; ++j)
        elementThermalForce[i] += coef*elC[i][j]*(ndTemps[j] - Tref);
    }*/
  }
//  for (int kk = 0; kk < 12; kk++)
//    fprintf(stderr, "eForce[%d] = %e\n", kk, elementThermalForce[kk]);
}

void
Tetrahedral::getGravityForce(CoordSet& cs,double *gravityAcceleration, 
                             Vector& gravityForce, int gravflg, GeomState *geomState)
{
  int nnodes = 4;
  double x[4], y[4], z[4];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  // Lumped
  if (gravflg != 2) {

    double ElementMassMatrix[12][12];
    double grvfor[3], totmas = 0.0;
    int grvflg = 1, masflg = 0;
    const int numdof = 12;
    _FORTRAN(mass23)(x, y, z, prop->rho, (double*)ElementMassMatrix, numdof,
                     gravityAcceleration, grvfor, grvflg, totmas, masflg);

    grvfor[0] /= double(nnodes);
    grvfor[1] /= double(nnodes);
    grvfor[2] /= double(nnodes);

    for(int i=0; i<nnodes; ++i) {
      gravityForce[3*i+0] = grvfor[0];
      gravityForce[3*i+1] = grvfor[1];
      gravityForce[3*i+2] = grvfor[2];
    }
  }
  // Consistent
  else {

    double lforce[4];
    for(int i=0; i<nnodes; ++i) lforce[i] = 0.0;

    // integration: loop over Gauss pts
    double w;
    double m[3], Shape[8], DShape[8][3];
    double dOmega;//det of jacobian

    // hard coded order 2 tetrahedral quadrature rule: {r,s,t,u(=1-r-s-t),w}
    int ngp = 4;
    double w1 = 1./24.;
    double r1 = (5.+3.*sqrt(5.))/20.;
    double s1 = (5.-  sqrt(5.))/20.;
    double t1 = s1;
    double u1 = 1.-r1-s1-t1;
    double TetGPt2[4][5] = {{r1, s1, t1, u1, w1},
                            {s1, t1, u1, r1, w1},
                            {t1, u1, r1, s1, w1},
                            {u1, r1, s1, t1, w1}};

    for(int igp = 0; igp < ngp; igp++) {
      // get x, y, z  position & weight of the integration pt
      m[0] = TetGPt2[igp][0]; m[1] = TetGPt2[igp][1];  m[2] = TetGPt2[igp][2]; w = TetGPt2[igp][4];
      dOmega = Tetra4ShapeFct(Shape, DShape, m, x, y, z);
      w *= prop->rho*fabs(dOmega);

      for(int n = 0; n < nnodes; ++n)
        lforce[n] += w*Shape[n];
    }
    
    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = lforce[i]*gravityAcceleration[0];
      gravityForce[3*i+1] = lforce[i]*gravityAcceleration[1];
      gravityForce[3*i+2] = lforce[i]*gravityAcceleration[2];
    }

  }
}

FullSquareMatrix
Tetrahedral::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
  //int status = 0;
  const int nnodes= 4;                                                                                           
  const int ndofs = 12;

  double X[4], Y[4], Z[4];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  FullSquareMatrix M(ndofs,mel);

  if(cmflg) { //HB: consistent mass matrix
    //fprintf(stderr," *** In Tetrahedral::massMatrix: make consistent mass matrix.\n");
    M.zero();
    int ls[12] = {0,3,6,9,1,4,7,10,2,5,8,11};
 
    // integration: loop over Gauss pts
    double /*x,y,z,wx,wy,wz,*/ w;
    double m[3], Shape[8], DShape[8][3];
    double dOmega;//det of jacobian
    int jSign = 0;

    // hard coded order 2 tetrahedral quadrature rule: {r,s,t,u(=1-r-s-t),w}
    int ngp = 4;
    double w1 = 1./24.;
    double r1 = (5.+3.*sqrt(5.))/20.;
    double s1 = (5.-  sqrt(5.))/20.;
    double t1 = s1;
    double u1 = 1.-r1-s1-t1;
    double TetGPt2[4][5] = {{r1, s1, t1, u1, w1},
                            {s1, t1, u1, r1, w1},
                            {t1, u1, r1, s1, w1},
                            {u1, r1, s1, t1, w1}};
                                                                                               
    for(int igp=0;igp<ngp;igp++){
      // get x, y, z  position & weight of the integration pt
      m[0] = TetGPt2[igp][0]; m[1] = TetGPt2[igp][1];  m[2] = TetGPt2[igp][2]; w = TetGPt2[igp][4];
      dOmega = Tetra4ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
      checkJacobian(&dOmega, &jSign, getGlNum()+1,"Tetrahedral::massMatrix");;
#endif
      w *= prop->rho*fabs(dOmega);  
      addNtDNtoM3DSolid(M, Shape, w, nnodes, ls);
    }

  } else { // Lumped mass matrix
    //fprintf(stderr," *** In Tetrahedral::massMatrix: make Lumped mass matrix.\n");
    double *gravityAcceleration = 0, *grvfor = 0, totmas = 0.0;
    int grvflg = 0, masflg = 0;
    _FORTRAN(mass23)(X, Y, Z, prop->rho, (double*)mel, ndofs,
  		    gravityAcceleration, grvfor, grvflg, totmas, masflg);
  }
                                                                                               
  return(M);
}

#ifndef USE_NEW_TET4_STIFF
// HB (04/15/05): replaced by new implementation (see below) to deal with anisotropic constitutive matrix
FullSquareMatrix
Tetrahedral::stiffness(CoordSet &cs, double *d, int flg)
{
  int status = 0;
  const int nnodes= 4;                                                                                           
  const int ndofs = 12;

  double X[4], Y[4], Z[4];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

   _FORTRAN(mstf23)(X, Y, Z, prop->E, prop->nu, (double*) d, ndofs,status);

   if(status != 0) {
     fprintf(stderr, "Error in Tetra4, nodes:");
     for(int i = 0; i < 4; ++i)
       fprintf(stderr, " %d",nn[i]+1);
     fprintf(stderr, "\n");
//   cerr << 0.0 << "," << 0.0 << "," << 0.0 << "  "
//   	  << X[1]-X[0] << "," << Y[1]-Y[0] << "," << Z[1]-Z[0] << "  "
//   	  << X[2]-X[0] << "," << Y[2]-Y[0] << "," << Z[2]-Z[0] << "  "
//   	  << X[3]-X[0] << "," << Y[3]-Y[0] << "," << Z[3]-Z[0] << endl;
   }

   FullSquareMatrix ret(12,d);
   return(ret);
}
#else
//---------------------------------------------------------------------------------
//HB (04/15/05)  new implementation of the Tetra4 stiffness matrix to deal
//               with anisotropic constitutive matrix                                                                                          
FullSquareMatrix
Tetrahedral::stiffness(CoordSet &cs,double *d, int flg)
{
  //fprintf(stderr, " *** WARNING: use new Tetrahedral::stiffness method (HB)\n");
  
  int status = 0;
  const int nnodes= 4;                                                                                           
  const int ndofs = 12;

  double X[4], Y[4], Z[4];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  int ls[12] = {0,3,6,9,1,4,7,10,2,5,8,11};
  FullSquareMatrix K(ndofs,d);
  K.zero();

  // get constitutive matrix
  double C[6][6];
  if(cCoefs) {  //anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
    //cerr<<" *** DEBUG: in Tetrahedral::stiffness, anisotropic/orthotropic material\n";
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
      
  // hard coded order 1 tetrahedral quadrature rule: {r,s,t,u(=1-r-s-t),w}
  double TetGPt1[1][5] = {{1./4.,1./4.,1./4.,1./4.,1./6.}};
  
  // integration: loop over Gauss pts
  int ngp = 1; // number of integration pts
  double m[3], Shape[4], DShape[4][3];
  double w, dOmega;//det of jacobian
  int jSign = 0;                                                                                                                                
  for(int igp=0;igp<ngp;igp++){
    // get x, y, z  position & weight of the integration pt
    m[0] = TetGPt1[igp][0]; m[1] = TetGPt1[igp][1];  m[2] = TetGPt1[igp][2]; w = TetGPt1[igp][4];
    dOmega = Tetra4ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
    checkJacobian(&dOmega, &jSign, getGlNum()+1, "Tetrahedral::stiffness");;
#endif
    w *= fabs(dOmega);
    addBtCBtoK3DSolid(K, DShape, C, w, nnodes, ls);
  }
  										     
  if(status != 0) {
    fprintf(stderr, " *** FATAL ERROR in Tetrahedral::stiffness");
    fprintf(stderr, " -> Tetra4's nodes: \n"); 
    for(int i = 0; i < nnodes; ++i)
       fprintf(stderr, " %d",nn[i]+1);
    fprintf(stderr, "\n");
    exit(-1);
  }
  return(K);
}      
#endif

int
Tetrahedral::numNodes() { return(4); }

int
Tetrahedral::numDofs() { return(12); }

int
Tetrahedral::getTopNumber() { 
  return 123; 
}


int*
Tetrahedral::nodes(int *p)
{
  if(!p) p = new int[4];
  p[0] = nn[0];
  p[1] = nn[1];
  p[2] = nn[2];
  p[3] = nn[3];
  return(p);
}


int*
Tetrahedral::dofs(DofSetArray &dsa, int *p)
{
  if(!p) p = new int[12];

  dsa.number(nn[0],DofSet::XYZdisp,p);
  dsa.number(nn[1],DofSet::XYZdisp,p+3);
  dsa.number(nn[2],DofSet::XYZdisp,p+6);
  dsa.number(nn[3],DofSet::XYZdisp,p+9);

  return(p);
}

void
Tetrahedral::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, 4, DofSet::XYZdisp);
}

Corotator *
Tetrahedral::getCorotator(CoordSet &cs, double *kel, int , int )
{
  if(mat) {
#ifdef USE_EIGEN3
    MatNLElement *ele = new NLTetrahedral4(nn);
    ele->setMaterial(mat);
    ele->setGlNum(glNum);
    return new MatNLCorotator(ele);
#else
    printf("WARNING: Corotator not implemented for element %d\n", glNum+1); return 0;
#endif
  }
  else {
    tetCorotator = new TetCorotator(nn, prop->E, prop->nu, cs);
    return tetCorotator;
  }
}

//---------------------------------------------------------------------------------
//HB (05/28/05)
// Stress evaluation in case of anisotropic elastic constitutive matrix 
void
Tetrahedral::getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs,
			      Vector &elDisp, int strInd, int surface, double *ndTemps,
			      double ylayer, double zlayer, int avgnum)
{
#ifdef TETRA4_DEBUG 
  fprintf(stderr," *** DEBUG: Get in Tetrahedral::getVonMisesAniso.\n");  
#endif

  const int nnodes =  4;
  //const int ndofs  = 12;
  weight = 1.0;
  
  double X[4], Y[4], Z[4];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
     
  // Flags sands17 to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = false; //HB: to force averaging the  Von Mises stress & strain. 
                         //    I don't really know the rational behind that, but its is necessary 
			 //    if we want to recover the same result as the old (fortran based) implementation
  
  double elStress[4][7];
  double elStrain[4][7];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef TETRA4_DEBUG 
    cerr<<" *** DEBUG: Tetrahedral::getVonMisesAniso, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[4][3] = {{0.0,0.0,0.0},{1.0,0.0,0.0},{ 0.0,1.0,0.0},{0.0,0.0,1.0}};

  double Shape[4], DShape[4][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Tetra4ShapeFct(Shape, DShape, m, X, Y, Z); 
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
Tetrahedral::getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
			        Vector &elDisp, int strInd, int surface, double *ndTemps)
{
#ifdef BRICK8_DEBUG 
  fprintf(stderr," *** DEBUG: Get in Tetrahedral::getVonMisesAniso.\n");  
#endif

  const int nnodes =  4;
  //const int ndofs  = 12;
  weight = 1.0;
  
  double X[4], Y[4], Z[4];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
        
  double elStress[4][6];
  double elStrain[4][6];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef TETRA4_DEBUG 
    cerr<<" *** DEBUG: Tetrahedral::getAllStressAniso, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[4][3] = {{0.0,0.0,0.0},{1.0,0.0,0.0},{ 0.0,1.0,0.0},{0.0,0.0,1.0}};

  double Shape[4], DShape[4][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Tetra4ShapeFct(Shape, DShape, m, X, Y, Z); 
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
  for(int i=0; i<nnodes; ++i){ // get average stress/strain  
    for(int j=0; j<6; ++j) 
      svec[j] = stress[i][j];
 
    // Convert Engineering to Tensor Strains
    if(strInd != 0) { svec[3] /= 2; svec[4] /= 2; svec[5] /= 2; }
    pstress(svec,pvec); // compute principal stress (or strain) & direction
    for(int j=0; j<3; ++j) 
      stress[i][j+6] = pvec[j];
  }
}

void
Tetrahedral::setMaterial(NLMaterial *_mat)
{
  mat = _mat;
}

int
Tetrahedral::numStates()
{
  int numGaussPoints = 1;
  return (mat) ? numGaussPoints*mat->getNumStates(): 0;
}

