#include <stdio.h>
#include <iostream>
#include <Element.d/Brick.d/EightNodeBrick.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Corotational.d/BrickCorotator.h>

//#define BRICK8_DEBUG

extern "C" {

void  _FORTRAN(brkcmt)(double&, double&, double*);

void  _FORTRAN(brik8v)(double*, double*, double*,double*,const int&,
                       double*, int &);
		       
void  _FORTRAN(vol17)(const int&,double*,double*,double*,double &);
void  _FORTRAN(br8mas)(const int&,double&,double*,double*,double*,
                        double*,double*,double*,const int&, double &);
			
void  _FORTRAN(sands17)(const int&,double*,double*,double*,double*,
                        double*,double*,double*,const int&,const int&,
			const int&,const int&,const int&);

void  _FORTRAN(sands17c)(const int&,double*,double*,double*,double*,
                         DComplex*,DComplex*,DComplex*,const int&,const int&,
                         const int&,const int&,const int&);

void _FORTRAN(hxgaus)(int &, int &, int &, int &, int &, int &,
                      double &,  double &, double &, double &);

void _FORTRAN(h8shpe)(double &, double &, double &, double *, double *,
                      double *, double *, double *, double *, double *,
                      double &);

void _FORTRAN(lgauss)(const int &, int &, double *, double *);
}

//HB: for anisotropic elastic constitutive matrix
void rotateConstitutiveMatrix(double *Cin, double *T33, double Cout[6][6]);

//HB: for consistent mass matrix 
double Hexa8ShapeFct(double Shape[8], double DShape[8][3], double m[3], double X[8], double Y[8], double Z[8]);
void addNtDNtoM3DSolid(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls, double (*D)[3] = 0);

//HB: for stresses & strains evaluation in case of ansitropic constitutive matrix
extern void computeStressAndEngStrain3DSolid(double Stress[6], double Strain[6], double C[6][6], double (*DShape)[3], double* U, int nnodes, int* ls=0);
extern double computeStress3DSolid(double Stress[6],double Strain[6], double C[6][6]);
extern double computeVonMisesStress(double Stress[6]);
extern double computeVonMisesStrain(double Strain[6]);

EightNodeBrick::EightNodeBrick(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nn[2] = nodenums[2];
  nn[3] = nodenums[3];
  nn[4] = nodenums[4];
  nn[5] = nodenums[5];
  nn[6] = nodenums[6];
  nn[7] = nodenums[7];

  brickCorotator = 0;

  cFrame = 0; 
  cCoefs = 0;
}

Element *
EightNodeBrick::clone()
{
  return(new EightNodeBrick(*this));
}

void
EightNodeBrick::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
  nn[6] = table[nn[6]];
  nn[7] = table[nn[7]];
}

void
EightNodeBrick::getVonMises(Vector& stress,Vector& weight,CoordSet &cs,
                            Vector& elDisp, int strInd, int surface, 
			    double *ndTemps, double ylayer, double zlayer, int avgnum)
{
  if(strInd == 17) { stress.zero(); return; }

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
  	Node &nd7 = cs.getNode(nn[6]);
  	Node &nd8 = cs.getNode(nn[7]);

  	double x[8], y[8], z[8], c[6][6];

  	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
  	x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
  	x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
  	x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
  	x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;

       if(cCoefs) {  // PJSA 3-30-05: orthotropic material
          // transform local constitutive matrix to global
          rotateConstitutiveMatrix(cCoefs, cFrame, c);
        }
        else { // isotropic material
          _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)c);
        }

        int vmflg = 0,strainFlg = 0;
	// Flags sands17 to calculate Von Mises stress.
        if(strInd == 6)  vmflg  = 1; 
	// Flags sands17 to calculate Von Mises strain
	if(strInd == 13) strainFlg = 1;

        int maxgus = 8; // maximum gauss points 
        int maxstr = 7; // maximum  
        int elm    = 1; // element number
        int maxsze = 1;

        double elStress[8][7], elStrain[8][7];

       _FORTRAN(sands17)(elm,x,y,z,(double*)c,elDisp.data(),
			(double*)elStress,(double*)elStrain,
                        maxgus,maxstr,maxsze,vmflg,strainFlg); 

        // HB (09-21-03): add thermal stresses contribution 
        // WARNING: thermal stress only supported for isotropic
        if(strInd <= 6) {
          if(cCoefs)
            fprintf(stderr," *** WARNING: EightNodeBrick::getVonMises: thermal stress evaluation NOT supported for anisotropic constitutive matrix. Use isotropic coeffs.\n");
          double thermalStress[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
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
            thermalStress[6] = coef*(ndTemps[6]-Tref);
            thermalStress[7] = coef*(ndTemps[7]-Tref);
          }
          stress[0] = elStress[0][strInd] - thermalStress[0];
          stress[1] = elStress[1][strInd] - thermalStress[1];
          stress[2] = elStress[2][strInd] - thermalStress[2];
          stress[3] = elStress[3][strInd] - thermalStress[3];
          stress[4] = elStress[4][strInd] - thermalStress[4];
          stress[5] = elStress[5][strInd] - thermalStress[5];
          stress[6] = elStress[6][strInd] - thermalStress[6];
          stress[7] = elStress[7][strInd] - thermalStress[7];
        }
        else {
          stress[0] = elStrain[0][strInd-7];
          stress[1] = elStrain[1][strInd-7];
          stress[2] = elStrain[2][strInd-7];
          stress[3] = elStrain[3][strInd-7];
          stress[4] = elStrain[4][strInd-7];
          stress[5] = elStrain[5][strInd-7];
          stress[6] = elStrain[6][strInd-7];
          stress[7] = elStrain[7][strInd-7];
        }
}

void
EightNodeBrick::getAllStress(FullM& stress,Vector& weight,CoordSet &cs,
Vector& elDisp, int strInd,int surface, double *ndTemps)
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
        Node &nd7 = cs.getNode(nn[6]);
        Node &nd8 = cs.getNode(nn[7]);

        double x[8], y[8], z[8], c[6][6];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
        x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
        x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
        x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
        x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;

       if(cCoefs) {  // PJSA 3-30-05: orthotropic material
          // transform local constitutive matrix to global
          rotateConstitutiveMatrix(cCoefs, cFrame, c);
        }
        else { // isotropic material
          _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)c);
        }

        int vmflg,strainFlg;
        vmflg  = 0;
        strainFlg = 0;

        int maxgus = 8; // maximum gauss points
        int maxstr = 7; // maximum
        int elm    = 1; // element number
        int maxsze = 1;

        double elStress[8][7], elStrain[8][7];

       _FORTRAN(sands17)(elm,x,y,z,(double*)c,elDisp.data(),
                        (double*)elStress,(double*)elStrain,
                        maxgus,maxstr,maxsze,vmflg,strainFlg);

// Store all Stress or all Strain as defined by strInd
        int i,j;
        if(strInd == 0) {
          // HB (09-21-03): add thermal stresses contribution
          double thermalStress[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	  if(ndTemps) {
            if(cCoefs)
              fprintf(stderr," *** WARNING: EightNodeBrick::getAllStress: thermal stress evaluation NOT supported for anisotropic constitutive matrix. Use isotropic coeffs.\n");
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
            thermalStress[6] = coef*(ndTemps[6]-Tref);
            thermalStress[7] = coef*(ndTemps[7]-Tref);
	  }
          for (i=0; i<8; ++i) {
            for (j=0; j<3; ++j) {
              stress[i][j] = elStress[i][j] - thermalStress[i];
            }
            for (j=3; j<6; ++j) {
              stress[i][j] = elStress[i][j];
            }
          }
        }
        else {
          for (i=0; i<8; ++i) {
            for (j=0; j<6; ++j) {
              stress[i][j] = elStrain[i][j];
            }
          }
        }

// Get Element Principals
        double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
        double pvec[3] = {0.0,0.0,0.0};
        for (j=0; j<6; ++j) {
          for (i=0; i<8; ++i) {
            svec[j] += stress[i][j];
          }
          svec[j] /= 8;
        }
// Convert Engineering to Tensor Strains
        if(strInd != 0) {
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
        }
        pstress(svec,pvec);
        for (i=0; i<8; ++i) {
          for (j=0; j<3; ++j) {
            stress[i][j+6] = pvec[j];
          }
        }
}

double
EightNodeBrick::getMass(CoordSet& cs)
{
  const int nnodes=  8;
  const int numgauss = 2;
 
  double X[8], Y[8], Z[8];                                                                                                                             
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  double volume=0.0;
  _FORTRAN(vol17)(numgauss,X,Y,Z,volume);
  double totmas = prop->rho*volume;

  return(totmas);
}

void
EightNodeBrick::getGravityForce(CoordSet& cs,double *gravityAcceleration, 
                                Vector& gravityForce, int gravflg, GeomState *geomState)
{
  const int nnodes=  8;
  int numgauss = 2;
 
  double X[8], Y[8], Z[8];                                                                                                                             
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  double grvfor[3];
  double ElementMassMatrix[24][24];
  
// Lumped
  if (gravflg == 1) {

    int grvflg = 1;
    double totmas = 0.0;

    _FORTRAN(br8mas)(numgauss,prop->rho,(double*)ElementMassMatrix,X, Y, Z,
                    gravityAcceleration,grvfor,grvflg,totmas);

// ... DISTRIBUTE GRAVITY FORCE AMONG NODES

    grvfor[0] *= 0.125;
    grvfor[1] *= 0.125;
    grvfor[2] *= 0.125;

    for(int i=0; i<nnodes; ++i) {
      gravityForce[3*i+0] = grvfor[0];
      gravityForce[3*i+1] = grvfor[1];
      gravityForce[3*i+2] = grvfor[2];
    }
  }
// Consistent
  else if (gravflg == 2) {
    double lforce[8];
    for(int i=0; i<nnodes; ++i) lforce[i] = 0.0;

    int fortran = 1;  // fortran routines start from index 1
    int pt1, pt2, pt3;
    for (pt1 = 0 + fortran; pt1 < numgauss + fortran; pt1++)  {
      for (pt2 = 0 + fortran; pt2 < numgauss + fortran; pt2++)  {
        for (pt3 = 0 + fortran; pt3 < numgauss + fortran; pt3++)  {
          // get gauss point
          double xi, eta, mu, wt;
          _FORTRAN(hxgaus)(numgauss, pt1, numgauss, pt2, numgauss, pt3, xi, eta, mu, wt);

          //compute shape functions
          double shapeFunc[8], shapeGradX[8], shapeGradY[8], shapeGradZ[8];
          double detJ;  //det of jacobian

          _FORTRAN(h8shpe)(xi, eta, mu, X, Y, Z,
                           shapeFunc, shapeGradX, shapeGradY, shapeGradZ, detJ);

          for (int i = 0; i < 8; ++i)
            lforce[i] += wt*shapeFunc[i]*detJ;
        }
      }
    }
    double rho = prop->rho;
    for(int i=0; i<nnodes; ++i) {
      gravityForce[3*i+0] = lforce[i]*rho*gravityAcceleration[0];
      gravityForce[3*i+1] = lforce[i]*rho*gravityAcceleration[1];
      gravityForce[3*i+2] = lforce[i]*rho*gravityAcceleration[2];
    }
  }
  else {
    for(int i=0; i<nnodes; ++i) {
      gravityForce[3*i+0] = 0.0;
      gravityForce[3*i+1] = 0.0;
      gravityForce[3*i+2] = 0.0;
    }
  }
}

FullSquareMatrix
EightNodeBrick::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
  //int status = 0;
  const int nnodes= 8;
  const int ndofs = 24;
  const int numgauss= 2;

  double X[8], Y[8], Z[8];                                                                                                                             
  cs.getCoordinates(nn, nnodes, X, Y, Z);
  
  FullSquareMatrix M(ndofs,mel);
  double *gravityAcceleration = 0, *grvfor = 0;

  int grvflg = 0;
  double totmas = 0.0;

  if(cmflg) { //HB: consistent mass matrix
    //fprintf(stderr," *** In EightNodeBrick::massMatrix: make consistent mass matrix.\n");
    M.zero();
    int ls[24] = {0,3,6,9,12,15,18,21,
                  1,4,7,10,13,16,19,22,
                  2,5,8,11,14,17,20,23};
                                                                                                                                       
    // integration: loop over Gauss pts
    double m[3], Shape[8], DShape[8][3];
    double wx,wy,wz,w, dOmega;//det of jacobian
#ifdef CHECK_JACOBIAN
    int jSign = 0; 
#endif
    for(int i=1;i<=numgauss;i++){
      _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
      for(int j=1;j<=numgauss;j++){
        _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
        for(int k=1;k<=numgauss;k++){
          _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
          dOmega = Hexa8ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
          checkJacobian(&dOmega, &jSign,"EightNodeBrick::massMatrix");
#endif
          w = fabs(dOmega)*wx*wy*wz*prop->rho;
          addNtDNtoM3DSolid(M, Shape, w, nnodes, ls);
        }
      }
    }
  } else { // Lumped mass matrix
    //fprintf(stderr," *** In EightNodeBrick::massMatrix: make Lumped mass matrix.\n");
    _FORTRAN(br8mas)(numgauss,prop->rho,(double*)mel,
                     X,Y,Z,gravityAcceleration,grvfor,grvflg,totmas);
  }

  return(M);
}


FullSquareMatrix
EightNodeBrick::stiffness(CoordSet &cs, double *d, int flg)
{
  int status = 0;
  const int nnodes= 8;
  const int ndofs = 24;
  const int numgauss= 2;

  double X[8], Y[8], Z[8];                                                                                                                             
  cs.getCoordinates(nn, nnodes, X, Y, Z);
 
  FullSquareMatrix K(ndofs,d);

  double C[6][6];
  if(cCoefs) {  // PJSA 3-30-05: orthotropic material
    // transform local constitutive matrix to global
    rotateConstitutiveMatrix(cCoefs, cFrame, C); 
  }
  else { // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C); 
  }
  _FORTRAN(brik8v)(X, Y, Z, (double *)C, numgauss, (double *)d, status);

  return(K);
}

int
EightNodeBrick::numNodes() { return(8); }

int
EightNodeBrick::numDofs() { return(24); }

int
EightNodeBrick::getTopNumber() { 
  return 117; 
}

int*
EightNodeBrick::nodes(int *p)
{
  if(!p) p = new int[numNodes()];
  p[0] = nn[0];
  p[1] = nn[1];
  p[2] = nn[2];
  p[3] = nn[3];
  p[4] = nn[4];
  p[5] = nn[5];
  p[6] = nn[6];
  p[7] = nn[7];
  return(p);
}


int*
EightNodeBrick::dofs(DofSetArray &dsa, int *p)
{
  if(!p) p = new int[numDofs()];

  dsa.number(nn[0],DofSet::XYZdisp, p  );
  dsa.number(nn[1],DofSet::XYZdisp, p+3);
  dsa.number(nn[2],DofSet::XYZdisp, p+6);
  dsa.number(nn[3],DofSet::XYZdisp, p+9);
  dsa.number(nn[4],DofSet::XYZdisp, p+12);
  dsa.number(nn[5],DofSet::XYZdisp, p+15);
  dsa.number(nn[6],DofSet::XYZdisp, p+18);
  dsa.number(nn[7],DofSet::XYZdisp, p+21);

  return(p);
}

void
EightNodeBrick::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

//-------------------------------------------------------------------
Corotator*
EightNodeBrick::getCorotator(CoordSet &cs, double *kel, int , int )
{
  return(new BrickCorotator(nn, prop->E, prop->nu, cs));
}

void EightNodeBrick::buildCorotator(CoordSet &cs)  
{
  if (!brickCorotator)
    brickCorotator = new BrickCorotator(nn, prop->E, prop->nu, cs);
}


//  HB (09-23-03): implement thermal force 
void
EightNodeBrick::getThermalForce(CoordSet &cs, Vector &ndTemps,
                                Vector &elementThermalForce, int glflag, 
				GeomState *geomState)
{
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes= 8;
  const int ndofs = 24;
  const int numgauss= 2;

  // extract nodes coordinates
  double X[8], Y[8], Z[8];                                                                                                                             
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
  double Shape[8], DShape[8][3], m[3];
  double wx,wy,wz,w,J; 
  if (geomState)  {
    if(cCoefs)
      fprintf(stderr," *** WARNING: Pentahedral::getThermalForce: evaluation of thermal forces NOT supported for anisotropic constitutive matrix. Use isotropic coeffs.\n");
     double dedU[24][6];
     for(int i=1;i<=numgauss;i++){
      _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
      for(int j=1;j<=numgauss;j++){
        _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
        for(int k=1;k<=numgauss;k++){
          _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
          brickCorotator->computeStrainGrad(*geomState, cs, dedU, i, j, k);
          J = Hexa8ShapeFct(Shape, DShape, m, X, Y, Z);
          w = fabs(J)*wx*wy*wz;
          double theta = 0.0;
          for(int inode=0; inode<nnodes; inode++) theta += Shape[inode]*(ndTemps[inode] - Tref);
          theta *= coef*J*w;

          for(int i=0; i<ndofs; i++)
            elementThermalForce[i] += theta*(dedU[i][0]+dedU[i][1]+dedU[i][2]);
        }
      }
    }
  }
  else {
    // integration: loop over Gauss pts
    for(int i=1;i<=numgauss;i++){
      _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
      for(int j=1;j<=numgauss;j++){
        _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
        for(int k=1;k<=numgauss;k++){
          _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
          J = Hexa8ShapeFct(Shape, DShape, m, X, Y, Z);
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
  //for(int i=0; i<24; i++) 
   //cerr<<" eelementThermalForce["<<i+1<<"] = "<<elementThermalForce[i]<<endl;
}

//---------------------------------------------------------------------------------
//HB (05/28/05)
// Stress evaluation in case of anisotropic elastic constitutive matrix                                                                                                                                                                                                             
void
EightNodeBrick::getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs,
			         Vector &elDisp, int strInd, int surface, double *ndTemps,
			         double ylayer, double zlayer, int avgnum)
{
#ifdef BRICK8_DEBUG 
  fprintf(stderr," *** DEBUG: Get in EightNodeBrick::getVonMisesAniso.\n");  
#endif

  const int nnodes =  8;
  //const int ndofs  = 24;
  weight = 1.0;
  
  double X[8], Y[8], Z[8];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
     
  // Flags sands17 to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = true; //HB: to force averaging the  Von Mises stress & strain. 
                         //    I don't really know the rational behind that, but its is necessary 
			 //    if we want to recover the same result as the old (fortran based) implementation
  
   
  double elStress[8][7];
  double elStrain[8][7];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef BRICK8_DEBUG 
    cerr<<" *** DEBUG: EightNodeBrick::getVonMisesAniso, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[8][3] = {{-1.0,-1.0,-1.0},{1.0,-1.0,-1.0},{ 1.0,1.0,-1.0},{-1.0,1.0,-1.0},
                               {-1.0,-1.0, 1.0},{1.0,-1.0, 1.0},{ 1.0,1.0, 1.0},{-1.0,1.0, 1.0}};

  double Shape[8], DShape[8][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Hexa8ShapeFct(Shape, DShape, m, X, Y, Z); 
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
  
  // compute average Von Mises stress and/or Von Mises strain: to match old Fortran code 
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
EightNodeBrick::getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
			          Vector &elDisp, int strInd, int surface, double *ndTemps)
{
#ifdef BRICK8_DEBUG 
  fprintf(stderr," *** DEBUG: Get in EightNodeBrick::getVonMisesAniso.\n");  
#endif

  const int nnodes =  8;
  //const int ndofs  = 24;
  weight = 1.0;
  
  double X[8], Y[8], Z[8];                                                                                                                              
  cs.getCoordinates(nn, nnodes, X, Y, Z);
        
  double elStress[8][6];
  double elStrain[8][6];
 
   // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic/orthotropic material
    // transform local constitutive matrix to global frame
#ifdef BRICK8_DEBUG 
    cerr<<" *** DEBUG: EightNodeBrick::getAllStressAniso, anisotropic/orthotropic material\n";
#endif
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else  // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
 
  //Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[8][3] = {{-1.0,-1.0,-1.0},{1.0,-1.0,-1.0},{ 1.0,1.0,-1.0},{-1.0,1.0,-1.0},
                               {-1.0,-1.0, 1.0},{1.0,-1.0, 1.0},{ 1.0,1.0, 1.0},{-1.0,1.0, 1.0}};

  double Shape[8], DShape[8][3];
  for(int inode=0; inode<nnodes; inode++){
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Hexa8ShapeFct(Shape, DShape, m, X, Y, Z); 
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
                                                                                     
                    

