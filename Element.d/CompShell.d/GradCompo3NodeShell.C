#include	<stdio.h>
#include	<stdlib.h>
#include 	<Utils.d/dbg_alloca.h>
#include 	<math.h>

#include	<Element.d/CompShell.d/Compo3NodeShell.h>
#include        <Corotational.d/Shell3Corotator.h>
#include        <Corotational.d/utilities.h>
#include        <Corotational.d/GeomState.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

extern "C"      {
void _FORTRAN(gxcompvms)(const int&, const int&, const int&, const int&, const int&,
                            double&,    double&,    double&,    double*,    double*,
                            double*,    double*,    double*,    double*,    double*,
			    double*,    double*,    double*,    double*,    double*,
			    double*, const int&, const int&, const int&,    double*,
			    double*,       int*,    double*,    double*,    double*,
			 const int&,       int*, const int&, const int&, const int&,
			 const int&, const int&, const int&, const int&, const int&);

void _FORTRAN(gxcompms1)(double*,    double*,    double*,       double*,       double*,
                         double*,    double*,    double*, const double&, const double&,
		         double*, const int&, const int&,    const int&,    const int&,
		            int*,    double*,    double*,       double*,    const int&,
		      const int&, const int&, const int&,       double*,       double*, 
		      const int&,    double&,    double&,    const int&);

void _FORTRAN(gxcompms2)(double*,    double*,    double*,       double*,       double*,
                         double*,    double*,    double*, const double&, const double&,
		         double*,    double*, const int&,    const int&,    const int&,    
	              const int&,       int*,    double*,       double*,       double*,
		      const int&, const int&, const int&,    const int&,       double*,
		         double*, const int&,    double&,    const int&);

void _FORTRAN(gxcompst)(const double&, const double&, const int&,       double*,    double*, 
                              double*,       double*, const int&, const double&,    double*,
			      double*,       double*,    double*,       double*,    double*,
			   const int&,    const int&, const int&,       double*,    double*,
			         int*,       double*,    double*,       double*, const int&,
	                   const int&,    const int&, const int&,          int&);
}

//------------------------------------------------------------------------------

void
Compo3NodeShell::getGradVonMises(Vector &dstress, Vector &weight, 
	                        CoordSet &cs, CoordSet &dcs, 
				Vector &elDisp, Vector &elGrad,
				int strInd, int surface)
{
	weight = 1.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);
        Node &dnd3 = dcs.getNode(nn[2]);

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
	
	double  emod = prop->E;
	double demod = gradprop->E;
	
	double    nu = prop->nu;

        const int maxstr =  7;
        const int maxgus =  3;

        double  elStress[3][7];
        double delStress[3][7];

        double *laysigm = 0;// new double[numLayers*maxstr*maxgus];

        int *laysgid = (int*)dbg_alloca(sizeof(int)*5*numLayers); 

        int i;
        for(i=0; i<5*numLayers; ++i)
          laysgid[i] = 0;

	int strainFlg;
	if(strInd > 6 ) strainFlg = 1;

	double * dcCoefs = new double[36];

        for (i=0;i<36;i++) dcCoefs[i]=0.0;	

       _FORTRAN(gxcompvms)(1,1,1,maxstr,maxgus,emod,demod,nu,h,dh,x,dx,y,dy,z,dz,
                           elDisp.data(),elGrad.data(),(double*)elStress,(double*)delStress,
                           (double*)laysigm,1,numLayers,1,(double*)cCoefs,dcCoefs,(int*)idlay,
			   layData,layGrad,cFrame,0,laysgid,1,type,1,1,numLayers,1,strainFlg,
			   surface);

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

        delete [] dcCoefs;
}

//------------------------------------------------------------------------------

double
Compo3NodeShell::getGradMass(CoordSet& cs, CoordSet& dcs)
{ 
       if (prop == 0) return 0.0;

        // gradient oriented data

        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);
        Node &dnd3 = dcs.getNode(nn[2]);

	double dx[3], dy[3], dz[3],  dh[3];

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
        dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;

        dh[0] = dh[1] = dh[2] = gradprop->eh;

	double drho = gradprop->rho;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double  x[3],  y[3],  z[3],   h[3];

        double ElementMassMatrix[18][18];

        double *gravityAcceleration=0, *grvfor=0;

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        h[0] =  h[1] =  h[2] = prop->eh;
      
	double rho  = prop->rho;

        int grvflg = 0, masflg = 1;

        double  totmas = 0.0;
        double dtotmas = 0.0;

       _FORTRAN(gxcompms1)(x,dx,y,dy,z,dz,h,dh,rho,drho,(double *)ElementMassMatrix,
                           18,numLayers,1,1,(int *)idlay, layData,layGrad,cFrame,
                           1,type,1,1,gravityAcceleration,grvfor, grvflg,
                           totmas,dtotmas,masflg);

	return dtotmas;
}

//------------------------------------------------------------------------------

void
Compo3NodeShell::gradMassMatrix(CoordSet &cs,CoordSet &dcs,
                               FullSquareMatrix &dmel)
{
        if (prop == 0) {
          dmel.zero();
          return;
        }

       // gradient orientated data

        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);
        Node &dnd3 = dcs.getNode(nn[2]);

        double dx[3], dy[3], dz[3], dh[3];

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
        dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;

        dh[0] = dh[1] = dh[2] = gradprop->eh;
	
	double drho = gradprop->rho;

        // non-gradient orientated data

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        h[0] = h[1] = h[2] = prop->eh;

	double rho  = prop->rho;

        double * mel = new double[18*18];
        
	int grvflg = 0, masflg = 0;

        double totmas = 0;

        double *gravityAcceleration=0, *grvfor=0;

       _FORTRAN(gxcompms2)(x,dx,y,dy,z,dz,h,dh,rho,drho,mel,dmel.data(),
                           18,numLayers,1,1,(int *)idlay, layData,layGrad,cFrame,
                           1,type,1,1,gravityAcceleration,grvfor, grvflg,
                           totmas,masflg);

        delete [] mel;
}

//------------------------------------------------------------------------------

void
Compo3NodeShell::gradstiffness(CoordSet &cs,CoordSet &dcs,FullSquareMatrix &dd,
                               int flg)
{
        if (prop == 0) {
	   dd.zero();
	   return;
        }

        // gradient orientated data

        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);
        Node &dnd3 = dcs.getNode(nn[2]);

        double dx[3], dy[3], dz[3], dh[3];

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
        dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;

        dh[0] = dh[1] = dh[2] = gradprop->eh;
	
	double demod = gradprop->E;

        // non-gradient orientated data

	Node &nd1 = cs.getNode(nn[0]);
	Node &nd2 = cs.getNode(nn[1]);
	Node &nd3 = cs.getNode(nn[2]);

	double x[3], y[3], z[3], h[3];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

	h[0] = h[1] = h[2] = prop->eh;

	double emod = prop->E;
	double nu   = prop->nu;

        double *       d = new double[18*18];
	double * dcCoefs = new double[36];

        int i;	
        for (i=0;i<36;i++) dcCoefs[i]=0.0;	

        _FORTRAN(gxcompst)(emod,demod,1,h,dh,(double *)d,dd.data(),18,nu,
                           x,dx,y,dy,z,dz,1,numLayers,1,(double*) cCoefs,
			   dcCoefs,(int *)idlay,layData,layGrad,
			   cFrame,1,type,1,1,flg);

        delete [] d;
	delete [] dcCoefs;
}

//------------------------------------------------------------------------------

int Compo3NodeShell::chkOptInf(CoordSet &dcs)
{
        if (prop == 0) return 0;

        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);
        Node &dnd3 = dcs.getNode(nn[2]);

        double dx[3], dy[3], dz[3], dh[3];

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
        dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;

        dh[0] = dh[1] = dh[2] = gradprop->eh;
	
	double demod = gradprop->E;

	double drho = gradprop->rho;

	double dvd = 0.0;
	
	dvd += dx[0]*dx[0] + dy[0]*dy[0] + dz[0]*dz[0];
	dvd += dx[1]*dx[1] + dy[1]*dy[1] + dz[1]*dz[1];
	dvd += dx[2]*dx[2] + dy[2]*dy[2] + dz[2]*dz[2];
		
	dvd += dh[0]*dh[0] + dh[1]*dh[1] + dh[2]*dh[2];
	dvd += demod*demod;
	dvd += drho*drho;
	
        int i;	
	for (i=0;i<numLayers;i++) {
	   dvd += layGrad[i*9]   * layGrad[i*9] 
	        + layGrad[i*9+1] * layGrad[i*9+1]
	        + layGrad[i*9+2] * layGrad[i*9+2]
	        + layGrad[i*9+3] * layGrad[i*9+3]
	        + layGrad[i*9+4] * layGrad[i*9+4]
	        + layGrad[i*9+5] * layGrad[i*9+5]
	        + layGrad[i*9+6] * layGrad[i*9+6]
	        + layGrad[i*9+7] * layGrad[i*9+7]
	        + layGrad[i*9+8] * layGrad[i*9+8];
	}

	return   (dvd > 0.0) ? 1 : 0;
}
