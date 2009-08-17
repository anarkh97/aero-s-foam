#include	<stdio.h>
#include        <math.h>
#include	<stdlib.h>
#include	<string.h>

#include	<Element.d/Shell.d/ThreeNodeShell.h>
#include        <Driver.d/PolygonSet.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Math.d/Vector.h>
#include        <Corotational.d/Shell3Corotator.h>
#include        <Corotational.d/utilities.h>
#include        <Corotational.d/GeomState.h>
#include	<Element.d/State.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>
#include        <Utils.d/pstress.h>
#include        <Hetero.d/InterpPoint.h>

// tria3d   - three node shell stiffness routine
// mass8    - three node shell mass routine
// sands8   - three node shell stress routine
// trithmfr - three node shell thermal load routine 

extern "C"      {
void _FORTRAN(tria3d)(int&, double* ,double* ,double* ,double& , double& ,
                      double* ,double*);
void _FORTRAN(mass8)(double* ,double* ,double* ,double* , double& ,
                     double* ,const int&, double* ,double*,const int&,
                     double&,const int&);
void _FORTRAN(sands8)(double*,double*,double*,double&,double&,double*,
                      double*,double*, const int&, const int&, const int&, 
                      const int&, const int&,const int&, double&);
void  _FORTRAN(trithmfr)(double*, double*, double*, double&, double&, double&,
 			 double&, double&, double&, double*, int&);
}

ThreeNodeShell::ThreeNodeShell(int* nodenums, double _w)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	w = _w;
}

Element *
ThreeNodeShell::clone()
{
	return new ThreeNodeShell(*this);
}

void
ThreeNodeShell::renum(int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
}

void
ThreeNodeShell::getVonMises(Vector& stress,Vector& weight,CoordSet &cs,
		       	    Vector& elDisp, int strInd,int surface,
                            double *ndTemps, double ylayer, double zlayer, int avgnum)
{ 
	weight = 1.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        // Thickness
        h[0] = h[1] = h[2] = prop->eh;
	
	//determine average nodal temperature difference relative to ambient
        double thermalStrain = 0.0;
        if(ndTemps) {
	  double dt = (ndTemps[0]+ndTemps[1]+ndTemps[2])/3 - prop->Ta;
	  thermalStrain = (prop->W)*dt;
        }

        int maxsze = 1;
        int maxstr = 7;
        int maxgus = 3;
        int elm    = 1;

        double elStress[3][7];

	// SET STRAIN FLAG IF USER WANTS STRAIN OUTPUT
	int strainFlg = 0;
	if( strInd > 6) strainFlg = 1;

       _FORTRAN(sands8)(x,y,z,prop->E,prop->nu,h,elDisp.data(),
                      (double*)elStress,
                      strainFlg, maxsze,maxstr,maxgus,elm,surface,thermalStrain);

        if(strInd < 7) {
          stress[0] = elStress[0][strInd];
          stress[1] = elStress[1][strInd];
          stress[2] = elStress[2][strInd];
        } else {
          stress[0] = elStress[0][strInd-7];
          stress[1] = elStress[1][strInd-7];
          stress[2] = elStress[2][strInd-7];
        }
}

void
ThreeNodeShell::getAllStress(FullM& stress,Vector& weight,CoordSet &cs,
                            Vector& elDisp, int strInd,int surface,
                            double *ndTemps)
{
        weight = 1.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        // Thickness
        h[0] = h[1] = h[2] = prop->eh;
	
	//determine average nodal temperature difference relative to ambient
	double dt = 0.0;
	if(ndTemps) // PJSA
          dt = (ndTemps[0]+ndTemps[1]+ndTemps[2])/3 - prop->Ta;
        double thermalStrain = (prop->W)*dt;

        int maxsze = 1;
        int maxstr = 7;
        int maxgus = 3;
        int elm    = 1;

        double elStress[3][7];

       _FORTRAN(sands8)(x,y,z,prop->E,prop->nu,h,elDisp.data(),
                      (double*)elStress,
                      strInd, maxsze,maxstr,maxgus,elm,surface,thermalStrain);

// Store all Stress or all Strain as defined by strInd
        int i,j;
        for (i=0; i<3; ++i) {
          for (j=0; j<6; ++j) {
            stress[i][j] = elStress[i][j];
          }
        }

// Get Element Principals
        double svec[6], pvec[3];
        for (j=0; j<6; ++j) {
          svec[j] = stress[0][j];
        }
// Convert Engineering to Tensor Strains
        if(strInd != 0) {
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
        }
        pstress(svec,pvec);
        for (i=0; i<3; ++i) {
          for (j=0; j<3; ++j) {
            stress[i][j+6] = pvec[j];
          }
        }
}

double
ThreeNodeShell::getMass(CoordSet& cs)
{
        if (prop == NULL) return 0.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double r1[3], r2[3], r3[3], v1[3], v2[3], v3[3];

        r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
        r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
        r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

        v1[0] = r3[0] - r1[0];
        v1[1] = r3[1] - r1[1];
        v1[2] = r3[2] - r1[2];

        v2[0] = r2[0] - r1[0];
        v2[1] = r2[1] - r1[1];
        v2[2] = r2[2] - r1[2];

        crossprod(v1, v2, v3);

        double area = 0.5*sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
        double mass = area*prop->rho*prop->eh;

        return mass;

}

void
ThreeNodeShell::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                                Vector& gravityForce, int gravflg, GeomState *geomState)
{
        if (prop == NULL) {
	   gravityForce.zero();
	   return;
        }

        double mass = getMass(cs);
        double massPerNode = mass/3.0;

        double fx = massPerNode*gravityAcceleration[0];
        double fy = massPerNode*gravityAcceleration[1];
        double fz = massPerNode*gravityAcceleration[2];
        double mx[3],my[3],mz[3];
        int i;   
        for(i=0; i<3; ++i) {
          mx[i]=0.0;
          my[i]=0.0;
          mz[i]=0.0;
        }

// Lumped
        if (gravflg == 1) {

// Consistent
//  Compute treating shell as 3 beams.
        }
        else if (gravflg == 2) {
          
          Node &nd1 = cs.getNode(nn[0]);
          Node &nd2 = cs.getNode(nn[1]);
          Node &nd3 = cs.getNode(nn[2]);
          double x[3], y[3], z[3];
          x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
          x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
          x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

          double T1[3],T2[3],T3[3];
          // Vector 1 from Node 1->2
          T1[0] = x[1] - x[0];
          T1[1] = y[1] - y[0];
          T1[2] = z[1] - z[0];
          normalize( T1 );
          // Vector 2 from Node 1->3
          T2[0] = x[2] - x[0];
          T2[1] = y[2] - y[0];
          T2[2] = z[2] - z[0];
          normalize( T2 );
          // Local Z-axis as cross between V1 and V2
          crossprod( T1, T2, T3 );
          normalize( T3);

          int beam, beamnode[3][2];
          beamnode[0][0] = 0;
          beamnode[0][1] = 1;
          beamnode[1][0] = 0;
          beamnode[1][1] = 2;
          beamnode[2][0] = 1;
          beamnode[2][1] = 2;

          for(beam=0; beam<3; ++beam) {
            double length, dx, dy, dz, localg[3];
            int n1, n2;
            n1 = beamnode[beam][0];
            n2 = beamnode[beam][1];
            dx = x[n2] - x[n1];
            dy = y[n2] - y[n1];
            dz = z[n2] - z[n1];
            length = sqrt(dx*dx + dy*dy + dz*dz);
            // Local X-axis from Node 1->2
            T1[0] = x[n2] - x[n1];
            T1[1] = y[n2] - y[n1];
            T1[2] = z[n2] - z[n1];
            normalize( T1 );
            // Local Y-axis as cross between Z and X
            crossprod( T3, T1, T2 );
            normalize( T2);

            for(i=0; i<3; ++i)
              localg[i] = 0.0;
            for(i=0; i<3; ++i) {
              localg[0] += T1[i]*gravityAcceleration[i];
              localg[1] += T2[i]*gravityAcceleration[i];
              localg[2] += T3[i]*gravityAcceleration[i];
            }
            double lmy,lmz;
            lmy = massPerNode*localg[2]*length/12.0;
            lmz = massPerNode*localg[1]*length/12.0;
            mx[n1] += ((T2[0]*lmy) + (T3[0]*lmz));
            my[n1] += ((T2[1]*lmy) + (T3[1]*lmz));
            mz[n1] += ((T2[2]*lmy) + (T3[2]*lmz));
            mx[n2] -= ((T2[0]*lmy) + (T3[0]*lmz));
            my[n2] -= ((T2[1]*lmy) + (T3[1]*lmz));
            mz[n2] -= ((T2[2]*lmy) + (T3[2]*lmz));
          }
        }
        else {
          fx = 0.0;
          fy = 0.0;
          fz = 0.0;
        }
        
        gravityForce[0]  =  fx;
        gravityForce[1]  =  fy;
        gravityForce[2]  =  fz;
        gravityForce[3]  =  mx[0];
        gravityForce[4]  =  my[0];
        gravityForce[5]  =  mz[0];
        gravityForce[6]  =  fx;
        gravityForce[7]  =  fy;
        gravityForce[8]  =  fz;
        gravityForce[9]  =  mx[1];
        gravityForce[10] =  my[1];
        gravityForce[11] =  mz[1];
        gravityForce[12] =  fx;
        gravityForce[13] =  fy;
        gravityForce[14] =  fz;
        gravityForce[15] =  mx[2];
        gravityForce[16] =  my[2];
        gravityForce[17] =  mz[2];

       //cerr << "shell gravityForce = "; for(int i=0; i<18; ++i) cerr << gravityForce[i] << " "; cerr << endl;
}

FullSquareMatrix
ThreeNodeShell::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
        // Check for phantom element, which has no stiffness
        if(prop == NULL) {
           FullSquareMatrix ret(18,mel);
	   ret.zero();
           return ret;
        }

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];
        double *gravityAcceleration = NULL, *grvfor = NULL,totmas = 0.0;

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        h[0] = h[1] = h[2] = prop->eh;

        int grvflg = 0, masflg = 0;

	const int numdof = 18;

       _FORTRAN(mass8)(x,y,z,h,prop->rho,(double *)mel,numdof,
               gravityAcceleration,grvfor,grvflg,totmas,masflg);

        FullSquareMatrix ret(18,mel);

        //cerr << "shell mass matrix = \n"; ret.print();
        //for(int i=0; i<18; ++i) if(ret[i][i] == 0.0) cerr << "oops-a-daisy\n";
 
        return ret;

}

FullSquareMatrix
ThreeNodeShell::stiffness(CoordSet &cs,double *d, int flg)
{
        // Check for phantom element, which has no stiffness
        if (prop == NULL) {
           FullSquareMatrix ret(18,d);
	   ret.zero();
           return ret;
        }

	Node &nd1 = cs.getNode(nn[0]);
	Node &nd2 = cs.getNode(nn[1]);
	Node &nd3 = cs.getNode(nn[2]);

	double x[3], y[3], z[3], h[3];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

	h[0] = h[1] = h[2] = prop->eh;

        // Check for a zero thickness
	if(h[0] <= 0.0) {
          fprintf(stderr," *** ERROR: Zero shell thickness "
                         "(ThreeNodeShell.C) %d %d %d\n",
                         nn[0]+1, nn[1]+1, nn[2]+1);
          fprintf(stderr," ... exiting fem program ...\n");
          exit(-1);
        }

        _FORTRAN(tria3d)(flg, x, y, z, prop->E, prop->nu, h, (double *)d);

        FullSquareMatrix ret(18,d);
       
        return ret;
}

int
ThreeNodeShell::numNodes()
{
 	return 3;
}

int*
ThreeNodeShell::nodes(int *p)
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
ThreeNodeShell::numDofs()
{
 	return 18;
}

int*
ThreeNodeShell::dofs(DofSetArray &dsa, int *p)
{
 	if(p == 0) p = new int[18];

        dsa.number(nn[0],DofSet::XYZdisp | DofSet::XYZrot, p  );
        dsa.number(nn[1],DofSet::XYZdisp | DofSet::XYZrot, p+6);
        dsa.number(nn[2],DofSet::XYZdisp | DofSet::XYZrot, p+12);

	return p;
}

void
ThreeNodeShell::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 3,  DofSet::XYZdisp | DofSet::XYZrot);
}

Corotator *
ThreeNodeShell::getCorotator(CoordSet &cs, double *kel, int fitAlgShell, int )
{
 int flag = 0; // signals stiffness routine to keep local matrix 
 FullSquareMatrix myStiff = stiffness(cs, kel, flag);
 return new Shell3Corotator(nn[0], nn[1], nn[2], myStiff, fitAlgShell);
}

void
ThreeNodeShell::computeDisp(CoordSet&, State &state, const InterpPoint &ip,
                            double *res, GeomState *gs)
{
 const double *gp = ip.xy;
 double xyz[3][6];
 state.getDV(nn[0], xyz[0], xyz[0]+3);
 state.getDV(nn[1], xyz[1], xyz[1]+3);
 state.getDV(nn[2], xyz[2], xyz[2]+3);

 int j;
 for(j=0; j<6; ++j)
    res[j] = (1.0-gp[0]-gp[1]) * xyz[0][j] + gp[0]*xyz[1][j] + gp[1]*xyz[2][j];

}

void
ThreeNodeShell::printInfo(CoordSet&, State &state, double gp[2])
{
 double xyz[3][6];
 state.getDV(nn[0], xyz[0], xyz[0]+3);
 state.getDV(nn[1], xyz[1], xyz[1]+3);
 state.getDV(nn[2], xyz[2], xyz[2]+3);

 fprintf(stderr, "INfo %e %e %e %e %e %e %e %e %e\n",
     xyz[0][0], xyz[0][1], xyz[0][2],
     xyz[1][0], xyz[1][1], xyz[1][2],
     xyz[2][0], xyz[2][1], xyz[2][2]);
}

void
ThreeNodeShell::getFlLoad(CoordSet &, const InterpPoint &ip, double *flF, 
                          double *resF, GeomState *gs)
{
 const double *gp = ip.xy;
 int i;
 for(i = 0; i < 3; ++i) {
   resF[i+3]  = resF[i+9] = resF[i+15] = 0.0;
   resF[i]    = (1.0-gp[0]-gp[1]) * flF[i];
   resF[6+i]  = gp[0] * flF[i];
   resF[12+i] = gp[1] * flF[i];
  }
}

int
ThreeNodeShell::getTopNumber()
{
  return 108;
}

void
ThreeNodeShell::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                     GeomState *geomState, int cflg)
{
     double px = 0.0;
     double py = 0.0;
     double pz = 0.0;

     double mx[3],my[3],mz[3];
     int i;
     for(i=0; i<3; ++i) {
       mx[i]=0.0;
       my[i]=0.0;
       mz[i]=0.0;
     }

     if (geomState) {

       // Get Nodes current coordinates
       NodeState &ns1 = (*geomState)[nn[0]];
       NodeState &ns2 = (*geomState)[nn[1]];
       NodeState &ns3 = (*geomState)[nn[2]];

       double  t0[3][3];
       double  xn[3][3];

       xn[0][0]  = ns1.x; // x coordinate of node state 1
       xn[0][1]  = ns1.y; // y coordinate of node state 1
       xn[0][2]  = ns1.z; // z coordinate of node state 1

       xn[1][0]  = ns2.x; // x coordinate of node state 2
       xn[1][1]  = ns2.y; // y coordinate of node state 2
       xn[1][2]  = ns2.z; // z coordinate of node state 2
      
       xn[2][0]  = ns3.x; // x coordinate of node state 3
       xn[2][1]  = ns3.y; // y coordinate of node state 3
       xn[2][2]  = ns3.z; // z coordinate of node state 3

       // Determine the area of the triangle
       //TM there is a problem here
       double xij[3][3];
       double yij[3][3];
       double zij[3][3];

       // Compute nodal delta coordinates
       int inod, jnod;
       for(inod=0; inod<3; inod++ ) {
         for(jnod=0; jnod<3; jnod++ ) {
           xij[inod][jnod] = xn[inod][0] - xn[jnod][0];
           yij[inod][jnod] = xn[inod][1] - xn[jnod][1];
           zij[inod][jnod] = xn[inod][2] - xn[jnod][2];
         }
       }

       /* TRIANGLE IN SPACE : WE COMPUTE THE LENGTH
        *.... OF ONE SIDE AND THE DISTANCE OF THE
        *.... OPPOSING NODE TO THAT SIDE TO COMPUTE THE AREA*/

       double rlr = sqrt( xij[1][0]*xij[1][0]+yij[1][0]*yij[1][0]+zij[1][0]*zij[1][0] );
       double rlb = sqrt( xij[2][1]*xij[2][1]+yij[2][1]*yij[2][1]+zij[2][1]*zij[2][1] );
       double bpr = sqrt((xij[1][0]*xij[2][1]+yij[1][0]*yij[2][1]+zij[1][0]*zij[2][1] ) 
                         *(xij[1][0]*xij[2][1]+yij[1][0]*yij[2][1]+zij[1][0]*zij[2][1] ))/rlr;

       double area= rlr*(sqrt(rlb*rlb-bpr*bpr))/2.0;

       // determine the normal to the plane of the element
       /* Compute t0 transformation matrix with x axis along side 1-2 */
       int i;
       for(i=0; i<3; i++ ) t0[0][i] = xn[1][i] - xn[0][i];
       normalize( t0[0] );

       // local z axis
       for( i=0; i<3; i++ ) t0[1][i] = xn[2][i] - xn[0][i];
       crossprod( t0[0], t0[1], t0[2] );
       normalize( t0[2] );

       double normal[3];

       normal[0] = t0[2][0];
       normal[1] = t0[2][1];
       normal[2] = t0[2][2];
 
       // compute pressure force per node
       double pressureForce = pressure * area / 3.0;

       px = pressureForce*normal[0];
       py = pressureForce*normal[1];
       pz = pressureForce*normal[2];

       if (cflg == 1) {

         int beam, beamnode[3][2];
         beamnode[0][0] = 0;
         beamnode[0][1] = 1;
         beamnode[1][0] = 0;
         beamnode[1][1] = 2;
         beamnode[2][0] = 1;
         beamnode[2][1] = 2;

         for(beam=0; beam<3; ++beam) {
           double length, dx, dy, dz;
           int n1, n2;
           n1 = beamnode[beam][0];
           n2 = beamnode[beam][1];
           dx = xn[n2][0] - xn[n1][0];
           dy = xn[n2][1] - xn[n1][1];
           dz = xn[n2][2] - xn[n1][2];
           length = sqrt(dx*dx + dy*dy + dz*dz);
           // Local X-axis from Node 1->2
           for(i=0; i<3; i++ ) t0[0][i] = xn[n2][i] - xn[n1][i];
           normalize( t0[0] );
           // Local Y-axis as cross between Z and X
           crossprod( t0[2], t0[0], t0[1] );
           normalize( t0[1] );
 
           double lmy = pressureForce*length/12.0;
           mx[n1] += (t0[1][0]*lmy);
           my[n1] += (t0[1][1]*lmy);
           mz[n1] += (t0[1][2]*lmy);
           mx[n2] -= (t0[1][0]*lmy);
           my[n2] -= (t0[1][1]*lmy);
           mz[n2] -= (t0[1][2]*lmy);
         }
       }

     }
     else {

       // Compute area of shell
       Node &nd1 = cs.getNode(nn[0]);
       Node &nd2 = cs.getNode(nn[1]);
       Node &nd3 = cs.getNode(nn[2]);

       double r1[3], r2[3], r3[3], v1[3], v2[3], normal[3];

       r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
       r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
       r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

       v1[0] = r3[0] - r1[0];
       v1[1] = r3[1] - r1[1];
       v1[2] = r3[2] - r1[2];

       v2[0] = r2[0] - r1[0];
       v2[1] = r2[1] - r1[1];
       v2[2] = r2[2] - r1[2];

       // Compute normal to shell using vector cross product rule
       crossprod(v2, v1, normal);

       double magnitude = sqrt(normal[0]*normal[0] + normal[1]*normal[1] 
                                                    + normal[2]*normal[2]);
       double area = 0.5*magnitude;

       // compute pressure force per node
       double pressureForce = pressure * area / 3.0;

       // compute unit normal to shell surface

       normal[0] /= magnitude;
       normal[1] /= magnitude;
       normal[2] /= magnitude;

       px = pressureForce*normal[0];
       py = pressureForce*normal[1];
       pz = pressureForce*normal[2];

       if (cflg == 1) {

         int beam, beamnode[3][2];
         beamnode[0][0] = 0;
         beamnode[0][1] = 1;
         beamnode[1][0] = 0;
         beamnode[1][1] = 2;
         beamnode[2][0] = 1;
         beamnode[2][1] = 2;

         double xn[3][3];
         for(i=0; i<3; ++i) {
           xn[0][i] = r1[i];
           xn[1][i] = r2[i];
           xn[2][i] = r3[i];
         }

         for(beam=0; beam<3; ++beam) {
           double length, dx, dy, dz;
           int n1, n2;
           n1 = beamnode[beam][0];
           n2 = beamnode[beam][1];
           dx = xn[n2][0] - xn[n1][0];
           dy = xn[n2][1] - xn[n1][1];
           dz = xn[n2][2] - xn[n1][2];
           length = sqrt(dx*dx + dy*dy + dz*dz);
           // Local X-axis from Node 1->2
           for(i=0; i<3; i++ ) v1[i] = xn[n2][i] - xn[n1][i];
           normalize( v1 );
           // Local Y-axis as cross between Z and X
           crossprod( normal, v1, v2 );
           normalize( v2 );

           double lmy = pressureForce*length/12.0;
           mx[n1] += (v2[0]*lmy);
           my[n1] += (v2[1]*lmy);
           mz[n1] += (v2[2]*lmy);
           mx[n2] -= (v2[0]*lmy);
           my[n2] -= (v2[1]*lmy);
           mz[n2] -= (v2[2]*lmy);
         }
       }
     }

     elPressureForce[0]  = px;
     elPressureForce[1]  = py;
     elPressureForce[2]  = pz;
     elPressureForce[3]  = mx[0];
     elPressureForce[4]  = my[0];
     elPressureForce[5]  = mz[0];

     elPressureForce[6]  = px;
     elPressureForce[7]  = py;
     elPressureForce[8]  = pz;
     elPressureForce[9]  = mx[1];
     elPressureForce[10] = my[1];
     elPressureForce[11] = mz[1];

     elPressureForce[12] = px;
     elPressureForce[13] = py;
     elPressureForce[14] = pz;
     elPressureForce[15] = mx[2];
     elPressureForce[16] = my[2];
     elPressureForce[17] = mz[2];
}

double
ThreeNodeShell::getMoment(Vector& force, CoordSet& cs, int node, int idir)
{
        if (prop == NULL) return 0.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        Node &nod = cs.getNode(node);

        double arm[3], ff[3], fm[3], mom[3];

        arm[0] = (nd1.x + nd2.x + nd3.x) / 3.0 - nod.x;
        arm[1] = (nd1.y + nd2.y + nd3.y) / 3.0 - nod.y;
	arm[2] = (nd1.z + nd2.z + nd3.z) / 3.0 - nod.z;
	
	ff[0] = (force[0] + force[6]  + force[12])/3.0;
	ff[1] = (force[1] + force[7]  + force[13])/3.0;
	ff[2] = (force[2] + force[8]  + force[14])/3.0;
	
	fm[0] = (force[3] + force[9]  + force[15])/3.0;
	fm[1] = (force[4] + force[10] + force[16])/3.0;
	fm[2] = (force[5] + force[11] + force[17])/3.0;
	

        crossprod(arm, ff, mom);

        double moment;
	
	switch (idir) {
	
	  case 3:
	    moment = mom[0] + fm[0];
	    break;
	  case 4:
	    moment = mom[1] + fm[1];
	    break;
	  case 5:
	    moment = mom[2] + fm[2];
	    break;
	  default:
	    moment = 0.0;
	}
        return moment;
}



/*/ dec
int
ThreeNodeShell::facelist(PolygonSet &pgs, int *flist)
{
 if(flist != 0) {
    flist[0] = pgs.addLine2( nn[0], nn[1] );
    flist[1] = pgs.addLine2( nn[1], nn[2] );
    flist[2] = pgs.addLine2( nn[2], nn[0] );
 }
 return 3 ;
}



// end dec */


//------------------------------------------------------------------------

void
ThreeNodeShell::getThermalForce(CoordSet& cs, Vector& ndTemps,
				Vector &elThermalForce,int glflag, 
				GeomState *geomState)
 {
    
   int i;
   
   double Tref = prop->Ta;
   
   double x[3], y[3], z[3];
   
   //get nodal temperatures and relate to reference temperature
   double meant = 0.0 ; //determine the average nodal temperature
   for (i=0;i<3;i++)
      meant += ndTemps[i]/3;
      meant -= Tref;


     if (geomState) {

       // Get Nodes current coordinates
       NodeState &ns1 = (*geomState)[nn[0]];
       NodeState &ns2 = (*geomState)[nn[1]];
       NodeState &ns3 = (*geomState)[nn[2]];   
       
       x[0] = ns1.x; y[0] = ns1.y; z[0] = ns1.z;
       x[1] = ns2.x; y[1] = ns2.y; z[1] = ns2.z;
       x[2] = ns3.x; y[2] = ns3.y; z[2] = ns3.z;
       
     } else {
             
       // Compute Node's coordinates
       Node &nd1 = cs.getNode(nn[0]);
       Node &nd2 = cs.getNode(nn[1]);
       Node &nd3 = cs.getNode(nn[2]);
      
       x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
       x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
       x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;      
    
     }
  
   
      double alpha = 1.5;
     
     
   //Call FORTRAN routine to determine elemental thermal force vector from
   //membranre effects -- returned in global coordinates if glflag = 0
   _FORTRAN(trithmfr)(x, y, z, meant, prop->W, prop->E, prop->nu,
		      prop->eh, alpha, (double *)elThermalForce.data(), 
                      glflag);
}

