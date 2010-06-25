#include <stdio.h>
#include <math.h>

#include <Element.d/Membrane.d/Membrane.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Corotational.d/utilities.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/Vector.h>
#include <Corotational.d/Shell3Corotator.h>

extern "C"      {
void _FORTRAN(trimem)(int&, double* ,double* ,double* ,double& , double& ,
                      double* ,double* );
void _FORTRAN(sands19)(double*,double*,double*,double&,double&,double*,double*,
        double*, const int&, const int&, const int&, const int&, const int&);
}

Membrane::Membrane(int* nodeNumbers)
{
	nn[0] = nodeNumbers[0];
	nn[1] = nodeNumbers[1];
	nn[2] = nodeNumbers[2];
}

Element *
Membrane::clone()
{
	return new Membrane(*this);
}

void
Membrane::renum(int *renumberingTable)
{
	nn[0] = renumberingTable[nn[0]];
	nn[1] = renumberingTable[nn[1]];
	nn[2] = renumberingTable[nn[2]];
}

void
Membrane::getVonMises(Vector& stress,Vector& weight,CoordSet &cs, Vector& elDisp, 
                      int strInd,int,double*,double ylayer, double zlayer, int avgnum)
{
	weight = 1.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        h[0] = h[1] = h[2] = prop->eh;

        double elStress[3][7];

        int  msize = 1;
        int maxstr = 7;
        int maxgus = 3;
        int    elm = 1; 

	// SET STRAIN FLAG IF USER WANTS STRAIN OUTPUT 
	int strainFlg;
	if(strInd > 6) strainFlg = 1; 

       _FORTRAN(sands19)(x,y,z,prop->E,prop->nu,h,elDisp.data(),(double*)elStress,
                         msize,maxstr,maxgus,elm,strainFlg);

	if(strInd < 7) {
          stress[0] = elStress[0][strInd];
          stress[1] = elStress[1][strInd];
          stress[2] = elStress[2][strInd];
	}
	else {
          stress[0] = elStress[0][strInd-7];
          stress[1] = elStress[1][strInd-7];
          stress[2] = elStress[2][strInd-7];
	}

}

void
Membrane::getAllStress(FullM& stress,Vector& weight,CoordSet &cs, Vector& elDisp, int strInd,int,double*)
{
	weight = 1.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        h[0] = h[1] = h[2] = prop->eh;

        double elStress[3][7];

        int  msize = 1;
        int maxstr = 7;
        int maxgus = 3;
        int    elm = 1; 

       _FORTRAN(sands19)(x,y,z,prop->E,prop->nu,h,elDisp.data(),(double*)elStress,
                         msize,maxstr,maxgus,elm,strInd);

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
Membrane::getMass(CoordSet& cs)
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

        v3[0] =  v1[1]*v2[2] - v2[1]*v1[2];
        v3[1] = -v1[0]*v2[2] + v2[0]*v1[2];
        v3[2] =  v1[0]*v2[1] - v2[0]*v1[1];

	double area = 0.5*sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
        double mass = area*prop->rho*prop->eh;

        return mass;

}

void
Membrane::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                          Vector& gravityForce, int gravflg, GeomState *geomState)
{
        double mass = getMass(cs);
        double massPerNode = mass/3.0;
        double fx, fy, fz;

        // Lumped
        if(gravflg != 2) {

          fx = massPerNode*gravityAcceleration[0];
          fy = massPerNode*gravityAcceleration[1];
          fz = massPerNode*gravityAcceleration[2];

        }
        // Consistent
        else {
          int i;
          Node &nd1 = cs.getNode(nn[0]);
          Node &nd2 = cs.getNode(nn[1]);
          Node &nd3 = cs.getNode(nn[2]);
          double x[3], y[3], z[3], localg[3];
          double T1[3],T2[3],T3[3];

          // Set the coordinates
          x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
          x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
          x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

          // Local X-axis
          T1[0] = x[1] - x[0];
          T1[1] = y[1] - y[0];
          T1[2] = z[1] - z[0];
          normalize( T1 );
          // 2nd Vector In Plane
          T2[0] = x[2] - x[0];
          T2[1] = y[2] - y[0];
          T2[2] = z[2] - z[0];
          normalize( T2 );
          // Local Z-axis as cross product of x-axis and in-plane vector
          crossprod( T1, T2, T3 );
          normalize( T3 );
          // Local Y-axis as cross product of x-axis and z-axis
          crossprod( T3, T1, T2 );
          normalize( T2 );

          for(i=0; i<3; ++i)
            localg[i] = 0.0;

          for(i=0; i<3; ++i) {
            localg[0] += T1[i]*gravityAcceleration[i];
            localg[1] += T2[i]*gravityAcceleration[i];
            localg[2] += T3[i]*gravityAcceleration[i];
          }
          double localf[3];
          localf[0] = massPerNode*localg[0];
          localf[1] = massPerNode*localg[1];
          localf[2] = massPerNode*localg[2];

          fx = (T1[0]*localf[0]) + (T2[0]*localf[1]) + (T3[0]*localf[2]);
          fy = (T1[1]*localf[0]) + (T2[1]*localf[1]) + (T3[1]*localf[2]);
          fz = (T1[2]*localf[0]) + (T2[2]*localf[1]) + (T3[2]*localf[2]);
        }

        gravityForce[0] =  fx;
        gravityForce[1] =  fy;
        gravityForce[2] =  fz;
        gravityForce[3] =  fx;
        gravityForce[4] =  fy;
        gravityForce[5] =  fz;
        gravityForce[6] =  fx;
        gravityForce[7] =  fy;
        gravityForce[8] =  fz;
}

FullSquareMatrix
Membrane::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        double x21   = x[1] - x[0];
        double y21   = y[1] - y[0];
        double z21   = z[1] - z[0];
        double x32   = x[2] - x[1];
        double y32   = y[2] - y[1];
        double z32   = z[2] - z[1];
        double x13   = x[0] - x[2];
        double y13   = y[0] - y[2];
        double z13   = z[0] - z[2];

        double rl[3];
        rl[0] = sqrt(x21*x21 + y21*y21 + z21*z21);
        rl[1] = sqrt(x32*x32 + y32*y32 + z32*z32);
        rl[2] = sqrt(x13*x13 + y13*y13 + z13*z13);

        double rmas = getMass(cs)/3.0; 

        FullSquareMatrix ret(9,mel);

        ret.zero();

        ret[0][0] = rmas;
        ret[1][1] = rmas;
        ret[2][2] = (rl[0]*rl[0]+rl[2]*rl[2])/420.0*rmas;
        ret[3][3] = rmas;
        ret[4][4] = rmas;
        ret[5][5] = (rl[1]*rl[1]+rl[0]*rl[0])/420.0*rmas;
        ret[6][6] = rmas;
        ret[7][7] = rmas;
        ret[8][8] = (rl[2]*rl[2]+rl[1]*rl[1])/420.0*rmas;

        return ret;
}

FullSquareMatrix
Membrane::stiffness(CoordSet &cs, double *d, int flg)
{
	Node &nd1 = cs.getNode(nn[0]);
	Node &nd2 = cs.getNode(nn[1]);
	Node &nd3 = cs.getNode(nn[2]);

	double x[3], y[3], z[3], h[3];

        // Set the coordinates
	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        // Set the thickness
	h[0] = h[1] = h[2] = prop->eh;

        if(h[0] <= 0.0)
          fprintf(stderr,"ERROR: Zero shell thickness (ThreeNodeShell.C) %d %d %d\n",
                nn[0], nn[1], nn[2]);

        _FORTRAN(trimem)(flg, x, y, z, prop->E, prop->nu, h, (double *)d);

        FullSquareMatrix ret(18,d);

        //cerr << "here in Membrane::stiffness\n"; ret.print();

        return ret;
}

int
Membrane::numNodes()
{
 	return 3;
}

int*
Membrane::nodes(int *p)
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
Membrane::numDofs()
{
 	return 18;
}

int*
Membrane::dofs(DofSetArray &dsa, int *p)
{
 	if(p == 0) p = new int[18];
        dsa.number(nn[0], DofSet::XYZdisp | DofSet::XYZrot, p  );
        dsa.number(nn[1], DofSet::XYZdisp | DofSet::XYZrot, p+6);
        dsa.number(nn[2], DofSet::XYZdisp | DofSet::XYZrot, p+12);
	return p;
}

void
Membrane::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 3, DofSet::XYZdisp | DofSet::XYZrot );
}

int
Membrane::getTopNumber()
{
  return 119;//4;
}

Corotator *
Membrane::getCorotator(CoordSet &cs, double *kel, int fitAlg, int)
{
 int flag = 0; // signals stiffness routine to keep local matrix 
 FullSquareMatrix myStiff = stiffness(cs, kel, flag);
 return new Shell3Corotator(nn[0], nn[1], nn[2], myStiff, fitAlg);
}

