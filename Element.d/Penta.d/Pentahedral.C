#include <cstdio>
#include <iostream>
#include <Element.d/Penta.d/Pentahedral.h>
#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Corotational.d/PentaCorotator.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/NLPentahedral.h>
#include <Corotational.d/MatNLCorotator.h>

#define USE_NEW_PENTA6_STIFF  //HB: force using new implementation of the stiffness matrix 
                              //    that deals with ansitropic constitutive matrix
#define CHECK_JACOBIAN        //HB: force check nullity & constant sign of jacobian over el.
//#define PENTA6_DEBUG

extern "C"      {

void	_FORTRAN(mstf24)(double*, double*, double*, double&, double&, double*, 
			 const int&, int &);
void    _FORTRAN(sands24)(const int&, double*, double*, double*, double&, 
			  double&, double*, double*, double*, const int&, const int&,
			  const int&, const int&, const int&,const int&);
// HB: for getThermalForce
void  _FORTRAN(lgauss)(int &, int &, double *, double *); 
#ifdef USE_NEW_PENTA6_STIFF //HB: used in the new version of the stiffness matrix
void _FORTRAN(brkcmt)(double&, double&, double*);
#endif
}

// HB: for getThermalForce & new version of the stiffness matrix
extern double Penta6ShapeFct(double Shape[6], double DShape[6][3], double m[3], double X[6], double Y[6], double Z[6]);
#ifdef USE_NEW_PENTA6_STIFF //HB: used in the new version of the stiffness matrix
extern void addBtCBtoK3DSolid(FullSquareMatrix &K, double (*DShape)[3], double C[6][6], double alpha, int nnodes, int* ls);
#endif
#ifdef CHECK_JACOBIAN //HB: for checking zero/small & constant sign of jacobian over the el.
extern int checkJacobian(double *J, int *jSign, int elId, const char* mssg= 0, double atol = 0.0, bool stop=true, FILE* file=stderr);
#endif
extern void addNtDNtoM3DSolid(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls, double (*D)[3] = 0);
extern void rotateConstitutiveMatrix(double *Cin, double *T33, double Cout[6][6]);
//HB: for stresses & strains evaluation in case of ansitropic constitutive matrix
extern void computeStressAndEngStrain3DSolid(double Stress[6], double Strain[6], double C[6][6], double (*DShape)[3], double* U, int nnodes, int* ls=0);
extern double computeStress3DSolid(double Stress[6],double Strain[6], double C[6][6]);
extern double computeVonMisesStress(double Stress[6]);
extern double computeVonMisesStrain(double Strain[6]);

Pentahedral::Pentahedral(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nn[2] = nodenums[2];
  nn[3] = nodenums[3];
  nn[4] = nodenums[4];
  nn[5] = nodenums[5];

  pentaCorotator = 0;

  cFrame = 0;
  cCoefs = 0;
  mat = 0;
}

Element *
Pentahedral::clone()
{
 return new Pentahedral(*this);
}

void
Pentahedral::renum(int *table)
{
     nn[0] = table[nn[0]];
     nn[1] = table[nn[1]];
     nn[2] = table[nn[2]];
     nn[3] = table[nn[3]];
     nn[4] = table[nn[4]];
     nn[5] = table[nn[5]];
}

void
Pentahedral::renum(EleRenumMap& table)
{
     nn[0] = table[nn[0]];
     nn[1] = table[nn[1]];
     nn[2] = table[nn[2]];
     nn[3] = table[nn[3]];
     nn[4] = table[nn[4]];
     nn[5] = table[nn[5]];
}

void
Pentahedral::getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
			 Vector &elDisp, int strInd, int surface, double *ndTemps,
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
	Node &nd5 = cs.getNode(nn[4]);
	Node &nd6 = cs.getNode(nn[5]);

  	double x[6], y[6], z[6]; 

  	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
        x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
	x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;

        // Flags sands17 to calculate Von Mises stress.
        int vmflg = 0;
        if(strInd == 6) vmflg  = 1; 

	// Flags sands17 to calculate Von Mises strain
        int strainFlg = 0;
        if(strInd == 13) strainFlg = 1;

        int maxgus = 6; // maximum gauss points 
        int maxstr = 7; // maximum  
        int elm    = 1; // element number
	int outerr = 6;

        double elStress[6][7];
        double elStrain[6][7];

	const int one = 1;
        _FORTRAN(sands24)(elm, x, y, z, prop->E, prop->nu, elDisp.data(),
			 (double*)elStress,(double*)elStrain,maxgus,
			  maxstr,one, outerr, vmflg, strainFlg); 

	// HB (09-21-03): add thermal stress contribution 
	if(strInd < 7) {
          double thermalStress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
          if(strInd == 0 || strInd == 1 || strInd == 2) {
            // get material props.
            double Tref  = prop->Ta;
            double alpha = prop->W ;
            double E     = prop->E ;
            double nu    = prop->nu;
            double coef = (E*alpha)/(1.0 - 2.0*nu);
            thermalStress[0] = coef*(ndTemps[0]-Tref);
            thermalStress[1] = coef*(ndTemps[1]-Tref);
            thermalStress[2] = coef*(ndTemps[2]-Tref);
            thermalStress[3] = coef*(ndTemps[3]-Tref);
            thermalStress[4] = coef*(ndTemps[4]-Tref);
            thermalStress[5] = coef*(ndTemps[5]-Tref);
	  } 
          stress[0] = elStress[0][strInd] - thermalStress[0];
          stress[1] = elStress[1][strInd] - thermalStress[1];
          stress[2] = elStress[2][strInd] - thermalStress[2];
          stress[3] = elStress[3][strInd] - thermalStress[3];
          stress[4] = elStress[4][strInd] - thermalStress[4];
          stress[5] = elStress[5][strInd] - thermalStress[5];
	}
	else {
          stress[0] = elStrain[0][strInd-7];
          stress[1] = elStrain[1][strInd-7];
          stress[2] = elStrain[2][strInd-7];
          stress[3] = elStrain[3][strInd-7];
          stress[4] = elStrain[4][strInd-7];
          stress[5] = elStrain[5][strInd-7];
	}
}

void
Pentahedral::getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
			 Vector &elDisp, int strInd, int surface, double *ndTemps)
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
	Node &nd5 = cs.getNode(nn[4]);
	Node &nd6 = cs.getNode(nn[5]);

  	double x[6], y[6], z[6]; 

  	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
        x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
	x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;

        int vmflg,strainFlg;
        vmflg  = 0;
        strainFlg = 0;

        int maxgus = 6; // maximum gauss points 
        int maxstr = 7; // maximum  
        int elm    = 1; // element number
	int outerr = 6;

        double elStress[6][7];
        double elStrain[6][7];

        double Tref  = prop->Ta; // ambient temperature
        double alpha = prop->W;  // thermal expansion coefficient
        double E     = prop->E;  // Young's modulus
        double nu    = prop->nu; // Poisson's ratio

	const int one = 1;
        _FORTRAN(sands24)(elm, x, y, z, E, nu, elDisp.data(),
			 (double*)elStress,(double*)elStrain,maxgus,
			  maxstr,one, outerr, vmflg, strainFlg); 

// Store all Stress or all Strain as defined by strInd
        int i,j;
        if(strInd == 0) {
          // HB (09-23-03): add thermal stresses contribution
          double thermalStress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
	  if(ndTemps) {
            double coef = (E*alpha)/(1.0 - 2.0*nu);
            thermalStress[0] = coef*(ndTemps[0]-Tref);
            thermalStress[1] = coef*(ndTemps[1]-Tref);
            thermalStress[2] = coef*(ndTemps[2]-Tref);
            thermalStress[3] = coef*(ndTemps[3]-Tref);
            thermalStress[4] = coef*(ndTemps[4]-Tref);
            thermalStress[5] = coef*(ndTemps[5]-Tref);
	  }

          for (i=0; i<6; ++i) {
            for (j=0; j<3; ++j) {
              stress[i][j] = elStress[i][j] - thermalStress[i];
            }
            for (j=3; j<6; ++j) {
              stress[i][j] = elStress[i][j];
            }
          }
        }
        else {
          for (i=0; i<6; ++i) {
            for (j=0; j<6; ++j) {
              stress[i][j] = elStrain[i][j];
            }
          }
        }

// Get Element Principals for each node without averaging
        double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
        double pvec[3] = {0.0,0.0,0.0};

        for (i=0; i<6; ++i) {
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
Pentahedral::getMass(CoordSet& cs)
{
  const int nnodes = 6;
  double X[6], Y[6], Z[6];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
  double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                          {2./3.,1./6.,1./6.,1./6.},
                          {1./6.,2./3.,1./6.,1./6.}};

  // integration: loop over Gauss pts
  double wxy, wz;
  int ngpz = 2;  // number of (linear) Gauss pts in the (local) z direction
  int ngpxy = 3; // numbder of (triangular) integration pts (in the local x-y plane)
  double m[3], Shape[6], DShape[6][3];
  double dOmega; // det of jacobian
  double volume = 0.0;
  for(int iz = 1; iz <= ngpz; iz++){
    // get z position & weight of the Gauss pt
    _FORTRAN(lgauss)(ngpz, iz, &m[2], &wz);
    for(int ixy = 0; ixy < ngpxy; ixy++) { // triangle Gauss pts
      // get x, y  position & weight of the Gauss pt
      m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
      dOmega = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
      volume += fabs(dOmega)*wxy*wz;
    }
  }

  return volume*prop->rho;
}

void
Pentahedral::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                             Vector& gravityForce, int gravflg, GeomState *geomState)

{
  const int nnodes=  6;

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

    const int ndofs = 18;

    double lforce[20];
    for(int i=0; i<nnodes; ++i) lforce[i] = 0.0;

    double X[6], Y[6], Z[6];
    cs.getCoordinates(nn, nnodes, X, Y, Z);

    // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
    double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                            {2./3.,1./6.,1./6.,1./6.},
                            {1./6.,2./3.,1./6.,1./6.}};

    // integration: loop over Gauss pts
    double wxy, wz, w;
    int ngpz = 2; // number of (linear) Gauss pts in the (local) z direction
    int ngpxy = 3; // numbder of (triangular) integration pts (in the local x-y plane)
    double m[3], Shape[6], DShape[6][3];
    double dOmega; // det of jacobian
    for(int iz = 1; iz <= ngpz; iz++){
      // get z position & weight of the Gauss pt
      _FORTRAN(lgauss)(ngpz, iz, &m[2], &wz);
      for(int ixy = 0; ixy < ngpxy; ixy++) { // triangle Gauss pts
        // get x, y  position & weight of the Gauss pt
        m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
        dOmega = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
        w = fabs(dOmega)*wxy*wz*prop->rho;

        for (int n = 0; n < nnodes; ++n) 
          lforce[n] += w*Shape[n];
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
Pentahedral::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  const int nnodes=  6;
  const int ndofs = 18;

  double X[6], Y[6], Z[6];
  cs.getCoordinates(nn, nnodes, X, Y, Z);
  
  FullSquareMatrix M(ndofs,mel);
  
  if(cmflg) { //HB: consistent mass matrix
    //fprintf(stderr," *** In Pentahedral::massMatrix: make consistent mass matrix.\n");
    M.zero();
    int ls[18] = {0,3,6, 9,12,15,
                  1,4,7,10,13,16,
  		  2,5,8,11,14,17};

    // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
    double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                            {2./3.,1./6.,1./6.,1./6.},
                            {1./6.,2./3.,1./6.,1./6.}};
                                                                                                                                                                                                                                                                         
    // integration: loop over Gauss pts
    double wxy, wz, w;
    int ngpz = 2; // number of (linear) Gauss pts in the (local) z direction
    int ngpxy= 3; // numbder of (triangular) integration pts (in the local x-y plane)
    double m[3], Shape[6], DShape[6][3];
    double dOmega; //det of jacobian
    int jSign = 0;
    for(int iz = 1; iz <= ngpz; iz++){ 
      // get z position & weight of the Gauss pt
      _FORTRAN(lgauss)(ngpz, iz, &m[2], &wz);
      for(int ixy = 0; ixy < ngpxy; ixy++) { // triangle Gauss pts
        // get x, y  position & weight of the Gauss pt
        m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
        dOmega = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
        checkJacobian(&dOmega, &jSign, getGlNum()+1, "Pentahedral::massMatrix");
#endif
        w = fabs(dOmega)*wxy*wz*prop->rho;
        addNtDNtoM3DSolid(M, Shape, w, nnodes, ls);
      }
    }
  } 
  else { // Lumped mass matrix
    fprintf(stderr," *** In Pentahedral::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return(M);
}

#ifndef USE_NEW_PENTA6_STIFF
// HB (04/15/05): replaced by new implementation (see below) to deal with anisotropic constitutive matrix
FullSquareMatrix
Pentahedral::stiffness(CoordSet &cs,double *d, int flg)
{
  int status = 0;
  const int nnodes=  6;	
  const int ndofs = 18;
  
  double X[6], Y[6], Z[6];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  _FORTRAN(mstf24)(X, Y, Z, prop->E, prop->nu, (double*) d, ndofs, status);

  if(status != 0) {
    fprintf(stderr, "Error in Penta:");
    int i;
    for(i = 0; i < 6; ++i)
       fprintf(stderr, " %d",nn[i]+1);
    fprintf(stderr, "\n");
  }

  FullSquareMatrix ret(18,d);
      
  return(ret);
}
#else
//---------------------------------------------------------------------------------
//HB (04/15/05)  new implementation of the Penta6 stiffness matrix to deal
//               with anisotropic constitutive matrix                                                                                          
FullSquareMatrix
Pentahedral::stiffness(CoordSet &cs,double *d, int flg)
{
  //fprintf(stderr, " *** WARNING: use new Pentahedral::stiffness method (HB)\n");
  
  int status = 0;
  const int nnodes=  6;	
  //const int ndofs = 18;

  double X[6], Y[6], Z[6];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
                                                                                            
  int ls[18] = {0,3,6, 9,12,15,
                1,4,7,10,13,16,
		2,5,8,11,14,17};
  FullSquareMatrix K(18,d);
  K.zero();

  // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
    //cerr<<" *** DEBUG: in Pentahedral::stiffness, anisotropic/orthotropic material\n";
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
  double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                          {2./3.,1./6.,1./6.,1./6.},
                          {1./6.,2./3.,1./6.,1./6.}};
                                                                                                                                  
  // integration: loop over Gauss pts
  double wxy,wz,w;
  int ngpz = 2; // number of (linear) Gauss pts in the (local) z direction
  int ngpxy= 3; // numbder of (triangular) integration pts (in the local x-y plane)
  double m[3], Shape[6], DShape[6][3];
  double dOmega;//det of jacobian
  int jSign = 0;                                                                                                                                
  for(int iz=1;iz<=ngpz;iz++){ 
    // get z position & weight of the Gauss pt
   _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz);
    for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
      // get x, y  position & weight of the Gauss pt
      m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
      dOmega = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
      checkJacobian(&dOmega, &jSign, getGlNum()+1, "Pentahedral::stiffness");
#endif
      w = fabs(dOmega)*wxy*wz;
      addBtCBtoK3DSolid(K, DShape, C, w, nnodes, ls);
    }
  }  
  										     
  if(status != 0) {
    fprintf(stderr, " *** FATAL ERROR in Pentahedral::stiffness");
    fprintf(stderr, " -> Pent6's nodes: \n"); 
    for(int i = 0; i < 6; ++i)
       fprintf(stderr, " %d",nn[i]+1);
    fprintf(stderr, "\n");
    exit(-1);
  }
  
  return(K);
}      
#endif

int
Pentahedral::numNodes()
{
 	return 6;
}

int*
Pentahedral::nodes(int *p)
{
  if(!p) p = new int[6];
  p[0] = nn[0];
  p[1] = nn[1];
  p[2] = nn[2];
  p[3] = nn[3];
  p[4] = nn[4];
  p[5] = nn[5];
  return(p);
}

int
Pentahedral::numDofs()
{
  return(18);
}

int*
Pentahedral::dofs(DofSetArray &dsa, int *p)
{
  if(!p) p = new int[18];

  dsa.number(nn[0],DofSet::XYZdisp,p);
  dsa.number(nn[1],DofSet::XYZdisp,p+3);
  dsa.number(nn[2],DofSet::XYZdisp,p+6);
  dsa.number(nn[3],DofSet::XYZdisp,p+9);
  dsa.number(nn[4],DofSet::XYZdisp,p+12);
  dsa.number(nn[5],DofSet::XYZdisp,p+15);

  return(p);
}

void
Pentahedral::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, 6, DofSet::XYZdisp);
}

int
Pentahedral::getTopNumber()
{
  return 124;
}

// HB (09-21-03): implementation of getThermalForce for 6 nodes pentahedral element
// Compute nodal thermal forces
void
Pentahedral::getThermalForce(CoordSet &cs, Vector &ndTemps,
                             Vector &elementThermalForce, int glflag, GeomState *geomState)
{
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes=  6;	
  const int ndofs = 18;
                                                                                                                   
  // extract nodes coordinates
  double X[6], Y[6], Z[6];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // initialize nodal thermal forces
  for(int i=0; i<ndofs; i++) elementThermalForce[i] = 0.0;

  // get material props & constitutive matrix
  double Tref  = prop->Ta;
  double alpha = prop->W ;
  double coef  = prop->E*(1.-2.*prop->nu);
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
  double Shape[6], DShape[6][3], m[3];
  double wxy,wz,w,J;
  int ngpz  = 2; // number of (linear) Gauss pts in the (local) z direction
  int ngpxy = 3; // numbder of (triangular) integration pts (in the local x-y plane)
  // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
  double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                          {2./3.,1./6.,1./6.,1./6.},
                          {1./6.,2./3.,1./6.,1./6.}};

  if (geomState)  { // GEOMETRICAL NONLINEAR ANALYSIS
    if(cCoefs)
      fprintf(stderr," *** WARNING: Pentahedral::getThermalForce: evaluation of thermal forces NOT supported for anisotropic constitutive matrix. Use isotropic coeffs.\n");
    double dedU[18][6];
    // integration: loop over Gauss pts
    for(int iz=1;iz<=ngpz;iz++){        // z Gauss pts
      // get z position & weight of the Gauss pt
      _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz);
      for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
        // get x, y  position & weight of the Gauss pt
        m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
        // compute shape fcts & their derivatives at the Gauss pt
        J = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
        w = wxy*wz*fabs(J);
        // compute strain gradient (total Lagrangian formulation)
        pentaCorotator->computeStrainGrad(*geomState, cs, dedU, m);
        // compute theta
        double theta = 0.0;
        for(int inode=0; inode<nnodes; inode++) theta += Shape[inode]*(ndTemps[inode] - Tref);
        // sum contribution
        theta *= coef*alpha*w;
        for(int i=0; i<ndofs; i++)
          elementThermalForce[i] += theta*(dedU[i][0]+dedU[i][1]+dedU[i][2]);
      }
    }
  } else { // GEOMETRICAL LINEAR ANALYSIS
    // integration: loop over Gauss pts
    for(int iz=1;iz<=ngpz;iz++){        // z Gauss pts
      // get z position & weight of the Gauss pt
      _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz);
      for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
        // get x, y  position & weight of the Gauss pt
        m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
        // compute shape fcts & their derivatives at the Gauss pt
        J = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
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
  //for(int i=0; i<18; i++){
  //  cerr<<" elementThermalForce["<<i<<"] = "<<elementThermalForce[i]<<endl;
  //}
}

//---------------------------------------------------------------------------------
Corotator*
Pentahedral::getCorotator(CoordSet &cs, double *kel, int , int )
{
  if(mat) {
#ifdef USE_EIGEN3
    MatNLElement *ele = new NLPentahedral6(nn);
    ele->setMaterial(mat);
    ele->setGlNum(glNum);
    return new MatNLCorotator(ele);
#else
    printf("WARNING: Corotator not implemented for element %d\n", glNum+1); return 0;
#endif
  }
  else {
    pentaCorotator = new PentaCorotator(nn, prop->E, prop->nu, cs);
    return pentaCorotator;
  }
}

//---------------------------------------------------------------------------------
//HB (05/28/05)
// Stress evaluation in case of anisotropic elastic constitutive matrix                                                                                                                                                                                                             
void
Pentahedral::getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs,
			      Vector &elDisp, int strInd, int surface, double *ndTemps,
			      double ylayer, double zlayer, int avgnum)
{
#ifdef PENTA6_DEBUG 
  fprintf(stderr," *** DEBUG: Get in Pentahedral::getVonMisesAniso.\n");  
#endif

  const int nnodes =  6;
  //const int ndofs  = 18;
  weight = 1.0;
  
  double X[6], Y[6], Z[6];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
     
  // Flags sands17 to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = false; //HB: to force averaging the  Von Mises stress & strain. 
                         //    I don't really know the rational behind that, but its is necessary 
			 //    if we want to recover the same result as the old (fortran based) implementation
   
  double elStress[6][7];
  double elStrain[6][7];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
    cerr<<" *** DEBUG: Pentahedral::getVonMises, anisotropic/orthotropic material\n";
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 

  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[6][3] = {{0.0,0.0,-1.0},{1.0,0.0,-1.0},{0.0,1.0,-1.0},
                               {0.0,0.0, 1.0},{1.0,0.0, 1.0},{0.0,1.0, 1.0}};

  double Shape[6], DShape[6][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta6ShapeFct(Shape, DShape, m, X, Y, Z); 
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
Pentahedral::getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
			       Vector &elDisp, int strInd, int surface, double *ndTemps)
{
#ifdef PENTA6_DEBUG 
  fprintf(stderr," *** DEBUG: Get in Pentahedral::getAllStressAniso.\n");  
#endif

  const int nnodes =  6;
  //const int ndofs  = 18;
  weight = 1.0;
  
  double X[6], Y[6], Z[6];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
  
  double elStress[6][6];
  double elStrain[6][6];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
    cerr<<" *** DEBUG: Pentahedral::getVonMises, anisotropic/orthotropic material\n";
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[6][3] = {{0.0,0.0,-1.0},{1.0,0.0,-1.0},{0.0,1.0,-1.0},
                               {0.0,0.0, 1.0},{1.0,0.0, 1.0},{0.0,1.0, 1.0}};

  double Shape[6], DShape[6][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta6ShapeFct(Shape, DShape, m, X, Y, Z); 
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

void
Pentahedral::setMaterial(NLMaterial *_mat)
{
  mat = _mat;
}

int
Pentahedral::numStates()
{
  int numGaussPoints = 6;
  return (mat) ? numGaussPoints*mat->getNumStates(): 0;
}

