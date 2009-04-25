#include	<Element.d/Shell.d/ThreeNodeShell.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Math.d/Vector.h>
#include        <Corotational.d/Shell3Corotator.h>
#include        <Corotational.d/utilities.h>
#include	<Element.d/State.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>
#include        <math.h>
#include	<stdlib.h>

// tria3d - three node shell stiffness routine
// mass8  - three node shell mass routine
// sands8 - three node shell stress routine

extern "C"      {
void _FORTRAN(gxtria3d)(int&, double* ,double* ,double* ,double* ,double* ,double* ,
                        double& ,double& ,double& ,
                        double* ,double* ,double* ,double*);

void _FORTRAN(gxmass8)(double* ,double* ,double* ,double* ,double* ,double* ,
                       double* ,double* ,double& ,double& ,
                       double* ,double* ,const int&, double* ,
		       double* ,const int& ,double& ,const int&);

void _FORTRAN(gxsands8)(double*,double*,double*,double*,double*,double*,
                        double&,double&,double&,double*,double*,
			double*,double*,
			double*,double*,
			const int&, const int&, const int&, 
                        const int&, const int&, const int&);
}

//------------------------------------------------------------------------------

void
ThreeNodeShell::getGradVonMises(Vector &dstress, Vector &weight, 
	                        CoordSet &cs, CoordSet &dcs, 
				Vector &elDisp, Vector &elGrad,
				int strInd, int surface)
{ 	
	weight = 1.0;

        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);

        Node dnd1 = dcs.getNode(nn[0]);
        Node dnd2 = dcs.getNode(nn[1]);
        Node dnd3 = dcs.getNode(nn[2]);

        double  x[3],  y[3],  z[3],  h[3];
        double dx[3], dy[3], dz[3], dh[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
	
        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
        dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;
	
         h[0] =  h[1] =  h[2] = prop->eh;
        dh[0] = dh[1] = dh[2] = gradprop->eh;

	double emod  = prop->E;
	double demod = gradprop->E;

	double nu  = prop->nu;

        int maxsze = 1;
        int maxstr = 7;
        int maxgus = 3;
        int elm    = 1;

        double  elStress[3][7];
        double delStress[3][7];

	// SET STRAIN FLAG IF USER WANTS STRAIN OUTPUT
	int strainFlg = 0;
	if( strInd > 6) strainFlg = 1;
               
	_FORTRAN( gxsands8)(x,dx,y,dy,z,dz,
                            emod,demod,nu,h,dh,
                            elDisp.data(),elGrad.data(),
		            (double*)elStress,(double*)delStress,
                            strainFlg,maxsze,maxstr,maxgus,elm,surface);

        if(strInd < 7) {
          dstress[0] = delStress[0][strInd];
          dstress[1] = delStress[1][strInd];
          dstress[2] = delStress[2][strInd];
        }
        else {
          dstress[0] = delStress[0][strInd-7];
          dstress[1] = delStress[1][strInd-7];
          dstress[2] = delStress[2][strInd-7];
        }
}

//------------------------------------------------------------------------------

double
ThreeNodeShell::getGradMass(CoordSet& cs, CoordSet& dcs)
{
        if (prop == NULL) return 0.0;

        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);

        Node dnd1 = dcs.getNode(nn[0]);
        Node dnd2 = dcs.getNode(nn[1]);
        Node dnd3 = dcs.getNode(nn[2]);
        
	double rho = prop->rho;
	double eh  = prop->eh;
        
	double drho = gradprop->rho;
	double deh  = gradprop->eh;

	double  r1[3],  r2[3],  r3[3],  v1[3],  v2[3],  v3[3];
	double dr1[3], dr2[3], dr3[3], dv1[3], dv2[3], dv3[3];

        r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
        r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
        r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

        dr1[0] = dnd1.x; dr1[1] = dnd1.y; dr1[2] = dnd1.z;
        dr2[0] = dnd2.x; dr2[1] = dnd2.y; dr2[2] = dnd2.z;
        dr3[0] = dnd3.x; dr3[1] = dnd3.y; dr3[2] = dnd3.z;

        v1[0] = r3[0] - r1[0];
        v1[1] = r3[1] - r1[1];
        v1[2] = r3[2] - r1[2];

        dv1[0] = dr3[0] - dr1[0];
        dv1[1] = dr3[1] - dr1[1];
        dv1[2] = dr3[2] - dr1[2];

        v2[0] = r2[0] - r1[0];
        v2[1] = r2[1] - r1[1];
        v2[2] = r2[2] - r1[2];

        dv2[0] = dr2[0] - dr1[0];
        dv2[1] = dr2[1] - dr1[1];
        dv2[2] = dr2[2] - dr1[2];

        dcrossprod(v1, v2, v3, dv1, dv2, dv3);

        double  area = 0.5*sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
	double darea = 0.25 *(v3[0]*dv3[0] + v3[1]*dv3[1] + v3[2]*dv3[2])
	             / area;

        double dmass = darea*rho*eh + area*drho*eh  + area*rho*deh;

        return dmass;

}

//------------------------------------------------------------------------------

void
ThreeNodeShell::gradMassMatrix(CoordSet &cs,CoordSet &dcs,
                               FullSquareMatrix &dmel)
{
        if (prop == 0) {
	   dmel.zero();
           return;
        }

       // gradient orientated data

        Node dnd1 = dcs.getNode(nn[0]);
        Node dnd2 = dcs.getNode(nn[1]);
        Node dnd3 = dcs.getNode(nn[2]);

        double dx[3], dy[3], dz[3], dh[3];

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
        dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;

        dh[0] = dh[1] = dh[2] = gradprop->eh;
	
	double drho = gradprop->rho;

        // non-gradient orientated data

        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        h[0] = h[1] = h[2] = prop->eh;
	
	double rho = prop->rho;

	const int numdof = 18;
	
        double * mel = new double[numdof*numdof];

        int grvflg = 0, masflg = 0;

        double *gravityAcceleration = NULL, *grvfor = NULL,totmas = 0.0;

       _FORTRAN(gxmass8)(x,dx,y,dy,z,dz,h,dh,rho,drho,
                         mel,dmel.data(),numdof,gravityAcceleration,
			 grvfor,grvflg,totmas,masflg);

        delete [] mel;
}

//------------------------------------------------------------------------------

void
ThreeNodeShell::gradstiffness(CoordSet &cs,CoordSet &dcs,FullSquareMatrix &dd,
                              int flg)
{
        if (prop == 0) {
	   dd.zero();
	   return;
        }

       // gradient orientated data

        Node dnd1 = dcs.getNode(nn[0]);
        Node dnd2 = dcs.getNode(nn[1]);
        Node dnd3 = dcs.getNode(nn[2]);

        double dx[3], dy[3], dz[3], dh[3];

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
        dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;

        dh[0] = dh[1] = dh[2] = gradprop->eh;
	
	double demod = gradprop->E;

        // non-gradient orientated data

 	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);

	double x[3], y[3], z[3], h[3];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

	h[0] = h[1] = h[2] = prop->eh;

	double emod = prop->E;
	double nu   = prop->nu;

        double * d = new double[18*18];

        _FORTRAN(gxtria3d)(flg,x,dx,y,dy,z,dz,emod,demod,nu,h,dh,d,dd.data());

        delete [] d;
}

//------------------------------------------------------------------------------

int ThreeNodeShell::chkOptInf(CoordSet &dcs)
{
        if (prop == 0) return 0;

        Node dnd1 = dcs.getNode(nn[0]);
        Node dnd2 = dcs.getNode(nn[1]);
        Node dnd3 = dcs.getNode(nn[2]);

        double dx[3], dy[3], dz[3], dh[3];

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
        dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;

        dh[0] = dh[1] = dh[2] = gradprop->eh;
	
	double drho = gradprop->rho;

	double demod = gradprop->E;

	double dvd = 0.0;
	
	dvd += dx[0]*dx[0] + dy[0]*dy[0] + dz[0]*dz[0];
	dvd += dx[1]*dx[1] + dy[1]*dy[1] + dz[1]*dz[1];
	dvd += dx[2]*dx[2] + dy[2]*dy[2] + dz[2]*dz[2];
		
	dvd += dh[0]*dh[0] + dh[1]*dh[1] + dh[2]*dh[2];
	dvd += drho*drho;
	dvd += demod;
	
	return   (dvd > 0.0) ? 1 : 0;
}
