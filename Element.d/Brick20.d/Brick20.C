#include <Element.d/Brick20.d/Brick20.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>

#define CHECK_JACOBIAN //HB: force check nullity & constant sign of jacobian over el.
//#define BRICK20_DEBUG

extern bool useFull;

extern "C"      {

void  _FORTRAN(brkcmt)(double&, double&, double*);

void  _FORTRAN(brik20v)(double*, double*, double*,double*,const int&,
                        double*, int &);

void   _FORTRAN(sands20)(const int&,double*,double*,double*,double*,
                         double*,double*,double*,const int&,const int&,
	 		 const int&,const int&,const int&);

void  _FORTRAN(br20vmint)(double*, double*, double*, double*, double*, 
	                  const int&, double&, double&, double&, double&, 
                          int&, int &);

void _FORTRAN(lgauss)(const int &, int &, double *, double *);
}

//HB: for anisotropic elastic constitutive matrix
void rotateConstitutiveMatrix(double *Cin, double *T33, double Cout[6][6]);

//HB: for consistent mass matrix
double Hexa20ShapeFct(double Shape[20], double DShape[20][3], double m[3], double X[20], double Y[20], double Z[20]);
void addNtDNtoM3DSolid(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls, double (*D)[3] = 0);
                                                                                               
#ifdef CHECK_JACOBIAN //HB: for checking zero/small & constant sign of jacobian over the el.
extern int checkJacobian(double *J, int *jSign, int elId, const char* mssg= 0, double atol = 0.0, bool stop=true, FILE* file=stderr);
#endif

//HB: for stresses & strains evaluation in case of ansitropic constitutive matrix
extern void computeStressAndEngStrain3DSolid(double Stress[6], double Strain[6], double C[6][6], double (*DShape)[3], double* U, int nnodes, int* ls=0);
extern double computeStress3DSolid(double Stress[6],double Strain[6], double C[6][6]);
extern double computeVonMisesStress(double Stress[6]);
extern double computeVonMisesStrain(double Strain[6]);

//-----------------------------------------------------------------------
//    Brick20 FEM input nodes numbering       Brick20 class nodes numbering
//                     16                                 20 
// 	       5+-------+-------+8		  5+-------+-------+8
// 	       /|	       /|		  /|		  /|		 
// 	      / |	      / |		 / |		 / |		 
// 	   13+  |	   15+  |	      17+  |	      19+  |		   
// 	    / 17+	    /	+20	       / 13+	       /   +16  	   
// 	   /	| 14	   /	|	      /    | 18       /    |		  
// 	 6+-------+-------+7	|	    6+-------+-------+7    |		  
// 	  |	|      12 |	|	     |     |	  12 |     |		  
// 	  |    1+-------+-|-----+ 4   ===>   |    1+-------+-|-----+ 4       
// 	  |    /	  |    /	     |    /	     |    /		  
// 	18+   / 	19+   / 	   14+   /	   15+   /		 
// 	  | 9+  	  |  +11	     | 9+	     |  +11	    
// 	  | /		  | /		     | /	     | /	    
// 	  |/		  |/		     |/ 	     |/ 		  
// 	 2+-------+-------+3		    2+-------+-------+3 		 
//   		 10		                    10
//-----------------------------------------------------------------------
Brick20::Brick20(int* nodenums)
{
  int i;
  for(i=0; i<12; ++i)
    nn[i] = nodenums[i];
  for(i =0; i < 4; ++i) {
    nn[i+16] = nodenums[i+12];
    nn[i+12] = nodenums[i+16];
  }

  cFrame = 0;
  cCoefs = 0;
}

Element *
Brick20::clone()
{
   return(new Brick20(*this));
}

void
Brick20::renum(int *table)
{
  for(int i=0; i<20; ++i)
    nn[i] = table[nn[i]];
}

void
Brick20::getVonMises(Vector& stress,Vector& weight,CoordSet &cs,
                     Vector& elDisp, int strInd, int, double*,
		     double ylayer, double zlayer, int avgnum)
{
	weight = 1.0;

  	double x[20], y[20], z[20], c[6][6];

        cs.getCoordinates(nn, numNodes(), x, y, z);

       if(cCoefs) {  // PJSA 3-30-05: orthotropic material
          // transform local constitutive matrix to global
          rotateConstitutiveMatrix(cCoefs, cFrame, c);
        }
        else { // isotropic material
          _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)c);
        }

        int vmflg     = 0;
        int strainFlg = 0;

	// Flags sands17 to calculate Von Mises stress.
        if(strInd == 6)  vmflg  = 1; 
	// Flags sands17 to calculate Von Mises strain
	if(strInd == 13) strainFlg = 1;

        int maxgus = 20; // number of nodes: 20 
        int maxstr = 7;  // maximum stress types  
        int elm    = 1;  // element number
        int maxsze = 1;

        double elStress[20][7], elStrain[20][7];

       _FORTRAN(sands20)(elm,x,y,z,(double*)c,elDisp.data(),
			(double*)elStress,(double*)elStrain,
                        maxgus,maxstr,maxsze,vmflg,strainFlg); 

        if(strInd <= 6) {
          stress[0]  = elStress[ 0][strInd];
          stress[1]  = elStress[ 1][strInd];
          stress[2]  = elStress[ 2][strInd];
          stress[3]  = elStress[ 3][strInd];
          stress[4]  = elStress[ 4][strInd];
          stress[5]  = elStress[ 5][strInd];
          stress[6]  = elStress[ 6][strInd];
          stress[7]  = elStress[ 7][strInd];
          stress[8]  = elStress[ 8][strInd];
          stress[9]  = elStress[ 9][strInd];
          stress[10] = elStress[10][strInd];
          stress[11] = elStress[11][strInd];
          stress[12] = elStress[12][strInd];
          stress[13] = elStress[13][strInd];
          stress[14] = elStress[14][strInd];
          stress[15] = elStress[15][strInd];
          stress[16] = elStress[16][strInd];
          stress[17] = elStress[17][strInd];
          stress[18] = elStress[18][strInd];
          stress[19] = elStress[19][strInd];
        }
        else {
          stress[0]  = elStrain[ 0][strInd-7];
          stress[1]  = elStrain[ 1][strInd-7];
          stress[2]  = elStrain[ 2][strInd-7];
          stress[3]  = elStrain[ 3][strInd-7];
          stress[4]  = elStrain[ 4][strInd-7];
          stress[5]  = elStrain[ 5][strInd-7];
          stress[6]  = elStrain[ 6][strInd-7];
          stress[7]  = elStrain[ 7][strInd-7];
          stress[8]  = elStrain[ 8][strInd-7];
          stress[9]  = elStrain[ 9][strInd-7];
          stress[10] = elStrain[10][strInd-7];
          stress[11] = elStrain[11][strInd-7];
          stress[12] = elStrain[12][strInd-7];
          stress[13] = elStrain[13][strInd-7];
          stress[14] = elStrain[14][strInd-7];
          stress[15] = elStrain[15][strInd-7];
          stress[16] = elStrain[16][strInd-7];
          stress[17] = elStrain[17][strInd-7];
          stress[18] = elStrain[18][strInd-7];
          stress[19] = elStrain[19][strInd-7];
        }

}

void
Brick20::getAllStress(FullM& stress,Vector& weight,CoordSet &cs,
Vector& elDisp, int strInd,int,double* )
{
	weight = 1.0;

  	double x[20], y[20], z[20], c[6][6];

        cs.getCoordinates(nn, numNodes(), x, y, z);

       if(cCoefs) {  // PJSA 3-30-05: orthotropic material
          // transform local constitutive matrix to global
          rotateConstitutiveMatrix(cCoefs, cFrame, c);
        }
        else { // isotropic material
          _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)c);
        }
                                                                                                                                        
        int vmflg     = 0;
        int strainFlg = 0;

        int maxgus = 8; // maximum gauss points 
        int maxstr = 7; // maximum  
        int elm    = 1; // element number
        int maxsze = 1;

        double elStress[8][7], elStrain[8][7];

       _FORTRAN(sands20)(elm,x,y,z,(double*)c,elDisp.data(),
			(double*)elStress,(double*)elStrain,
                        maxgus,maxstr,maxsze,vmflg,strainFlg); 

// Store all Stress or all Strain as defined by strInd
        int i,j;
        if(strInd == 0) {
          for (i=0; i<8; ++i) {
            for (j=0; j<6; ++j) {
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
        pstress(svec,pvec);
        for (i=0; i<8; ++i) {
          for (j=0; j<3; ++j) {
            stress[i][j+6] = pvec[j];
          }
        }
}

void
Brick20::getVonMisesInt(CoordSet &cs, Vector &d, double& sigbar, 
                        double& fac, int areaFlag,
                        double& vmint, double& vol)
{
        Node &nd1  = cs.getNode(nn[0]);
        Node &nd2  = cs.getNode(nn[1]);
        Node &nd3  = cs.getNode(nn[2]);
        Node &nd4  = cs.getNode(nn[3]);
        Node &nd5  = cs.getNode(nn[4]);
        Node &nd6  = cs.getNode(nn[5]);
        Node &nd7  = cs.getNode(nn[6]);
        Node &nd8  = cs.getNode(nn[7]);
        Node &nd9  = cs.getNode(nn[8]);
        Node &nd10 = cs.getNode(nn[9]);
        Node &nd11 = cs.getNode(nn[10]);
        Node &nd12 = cs.getNode(nn[11]);
        Node &nd13 = cs.getNode(nn[12]);
        Node &nd14 = cs.getNode(nn[13]);
        Node &nd15 = cs.getNode(nn[14]);
        Node &nd16 = cs.getNode(nn[15]);
        Node &nd17 = cs.getNode(nn[16]);
        Node &nd18 = cs.getNode(nn[17]);
        Node &nd19 = cs.getNode(nn[18]);
        Node &nd20 = cs.getNode(nn[19]);

	double x[20], y[20], z[20], c[6][6];

        x[0]  = nd1.x;   y[0]  = nd1.y;   z[0]  = nd1.z;
        x[1]  = nd2.x;   y[1]  = nd2.y;   z[1]  = nd2.z;
        x[2]  = nd3.x;   y[2]  = nd3.y;   z[2]  = nd3.z;
        x[3]  = nd4.x;   y[3]  = nd4.y;   z[3]  = nd4.z;
        x[4]  = nd5.x;   y[4]  = nd5.y;   z[4]  = nd5.z;
        x[5]  = nd6.x;   y[5]  = nd6.y;   z[5]  = nd6.z;
        x[6]  = nd7.x;   y[6]  = nd7.y;   z[6]  = nd7.z;
        x[7]  = nd8.x;   y[7]  = nd8.y;   z[7]  = nd8.z;
        x[8]  = nd9.x;   y[8]  = nd9.y;   z[8]  = nd9.z;
        x[9]  = nd10.x;  y[9]  = nd10.y;  z[9]  = nd10.z;
        x[10] = nd11.x;  y[10] = nd11.y;  z[10] = nd11.z;
        x[11] = nd12.x;  y[11] = nd12.y;  z[11] = nd12.z;
        x[12] = nd13.x;  y[12] = nd13.y;  z[12] = nd13.z;
        x[13] = nd14.x;  y[13] = nd14.y;  z[13] = nd14.z;
        x[14] = nd15.x;  y[14] = nd15.y;  z[14] = nd15.z;
        x[15] = nd16.x;  y[15] = nd16.y;  z[15] = nd16.z;
        x[16] = nd17.x;  y[16] = nd17.y;  z[16] = nd17.z;
        x[17] = nd18.x;  y[17] = nd18.y;  z[17] = nd18.z;
        x[18] = nd19.x;  y[18] = nd19.y;  z[18] = nd19.z;
        x[19] = nd20.x;  y[19] = nd20.y;  z[19] = nd20.z;

	_FORTRAN(brkcmt)(prop->E, prop->nu, (double*)c);

  	const int numgauss = 3;
	int status;

        _FORTRAN(br20vmint)(x, y, z, d.data(), (double *)c, numgauss, 
	                  vmint, vol, sigbar, fac, areaFlag, status);
}

double
Brick20::getMass(CoordSet& cs)
{
  double x[20], y[20], z[20];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  const int numgauss = 3;
  double wx, wy, wz;
  double m[3], Shape[20], DShape[20][3];
  double dOmega; // det of jacobian
  double volume = 0.0;
  for(int i = 1; i <= numgauss; i++) {
    _FORTRAN(lgauss)(numgauss, i, &m[0], &wx);
    for(int j = 1; j <= numgauss; j++) {
      _FORTRAN(lgauss)(numgauss, j, &m[1], &wy);
      for(int k = 1; k <= numgauss; k++) {
        _FORTRAN(lgauss)(numgauss, k, &m[2], &wz);
        dOmega = Hexa20ShapeFct(Shape, DShape, m, x, y, z);
        volume += fabs(dOmega)*wx*wy*wz;
      }
    }
  }
  double totmas = prop->rho*volume;

  return totmas;
}

void
Brick20::getGravityForce(CoordSet& cs,double *gravityAcceleration, 
                         Vector& gravityForce, int gravflg, GeomState *geomState)
{
  int nnodes = 20;

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

    double lforce[20];
    for(int i=0; i<nnodes; ++i) lforce[i] = 0.0;

    double x[20], y[20], z[20];
    cs.getCoordinates(nn, numNodes(), x, y, z);

    int numgauss = 3;
    // integration: loop over Gauss pts
    double wx, wy, wz, w;
    double m[3], Shape[20], DShape[20][3];
    double dOmega; // det of jacobian
    for(int i = 1; i <= numgauss; i++) {
      _FORTRAN(lgauss)(numgauss, i, &m[0], &wx);
      for(int j = 1; j <= numgauss; j++) {
        _FORTRAN(lgauss)(numgauss, j, &m[1], &wy);
        for(int k = 1; k <= numgauss; k++) {
          _FORTRAN(lgauss)(numgauss, k, &m[2], &wz);

          dOmega = Hexa20ShapeFct(Shape, DShape, m, x, y, z);
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
Brick20::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  double X[20], Y[20], Z[20];
  cs.getCoordinates(nn, numNodes(), X, Y, Z);

  const int numgauss = 3;
  const int nnodes = 20;
  const int ndofs = 60;
  FullSquareMatrix M(ndofs,mel);

  if(cmflg) { //HB: consistent mass matrix
    //fprintf(stderr," *** In Brick20::massMatrix: make consistent mass matrix.\n");
    M.zero();
    int ls[60] = {0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,
                  1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,
                  2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59}; 
                                                                                               
    // integration: loop over Gauss pts
    double wx,wy,wz,w;
    double m[3], Shape[20], DShape[20][3];
    double dOmega;//det of jacobian
    int jSign = 0;                                                                                                
    for(int i=1;i<=numgauss;i++){
      _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
      for(int j=1;j<=numgauss;j++){
        _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
        for(int k=1;k<=numgauss;k++){
          _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
          dOmega = Hexa20ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
          checkJacobian(&dOmega, &jSign, getGlNum()+1, "Brick20::massMatrix");
#endif
          w = fabs(dOmega)*wx*wy*wz*prop->rho;
          addNtDNtoM3DSolid(M, Shape, w, nnodes, ls);
        }
      }
    }
  } else { // Lumped mass matrix
    fprintf(stderr," *** In Brick20::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return(M);
}


FullSquareMatrix
Brick20::stiffness(CoordSet &cs, double *d, int flg)
{
  double X[20], Y[20], Z[20], C[6][6];
  cs.getCoordinates(nn, numNodes(), X, Y, Z);

  if(cCoefs) {  // HB 06-19-05: anisotropic constitutive matrix 
    // transform local constitutive matrix to global
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  }
  else { // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
  }

  const int numgauss = 3;
  int status;

  _FORTRAN(brik20v)(X, Y, Z, (double *)C, numgauss, (double *)d, status);

  FullSquareMatrix K(60,d);

  return(K);
}

int
Brick20::numNodes(){ 
  if(useFull)
    return(20); 
  else
    return 8;
}

int
Brick20::numDofs(){ return(60); }

int
Brick20::getTopNumber()
{
  return(172);//9; // Cube Geometry, watch: it is returned as 8 node brick
}

int*
Brick20::nodes(int *p)
{
  if(useFull)
    {
      if(!p) p = new int[20];
      
      for(int i=0; i<12; ++i)
	p[i] = nn[i];
      
      for(int i =0; i < 4; ++i) {
	p[i+12] = nn[i+16];
	p[i+16] = nn[i+12];
      }
      
      return(p);
    }
  else
    {
        if(!p) p = new int[8];
      
      for(int i=0; i<8; ++i)
	p[i] = nn[i];
      
      
      return(p);
    }
}


int*
Brick20::dofs(DofSetArray &dsa, int *p)
{
  if(!p) p = new int[60];

  for(int i=0; i<20; ++i)
    dsa.number(nn[i],DofSet::XYZdisp, p+3*i);

  return(p);
}

// 3D element with x,y,z degrees of freedom per node
void
Brick20::markDofs(DofSetArray &dsa)
{
  for(int i=0; i<20; ++i)
    dsa.mark(nn[i],DofSet::XYZdisp);
}


// HB (09-21-03): implement WARNING message
void
Brick20::getThermalForce(CoordSet &cs, Vector &ndTemps,
                         Vector &elementThermalForce, int glflag, GeomState *geomState)
{
  fprintf(stderr," *** WARNING: Brick20::getThermalForce NOT implemented. Return NULL force.\n");
  elementThermalForce.zero();
  return;

  //HB: NOT TESTED, VALIDATED YET ...
#ifdef BRICK20_DEBUG 
  fprintf(stderr," *** DEBUG: Get in Brick20::getThermalForce.\n");  
#endif
 
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes= 20;
  const int ndofs = 60;
  const int numgauss= 2;

  // extract nodes coordinates
  double X[20], Y[20], Z[20];                                                                                                                             
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
#ifdef BRICK20_DEBUG 
    cerr<<" *** DEBUG: Brick20::getThermalForc, anisotropic/orthotropic material\n";
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
  double Shape[20], DShape[20][3], m[3];
  double wx,wy,wz,w,J; 
  if (geomState)  {
    fprintf(stderr," *** ERROR: Brick20::getThermalForce NOT supported for geometric nonlinear analysis. Abort.\n");
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
          J = Hexa20ShapeFct(Shape, DShape, m, X, Y, Z);
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
  //for(int i=0; i<60; i++) 
   //cerr<<" eelementThermalForce["<<i+1<<"] = "<<elementThermalForce[i]<<endl;
}

