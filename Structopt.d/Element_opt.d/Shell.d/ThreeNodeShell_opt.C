#ifdef STRUCTOPT

#include <Structopt.d/Optvar.h>
#include <Structopt.d/Element_opt.d/Shell.d/ThreeNodeShell_opt.h>
#include <Corotational.d/utilities.h>

//------------------------------------------------------------------------------
int ThreeNodeShell_opt::chkOptInf(CoordSet &dcs)
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
  
  double drho  = gradprop->rho;
  double demod = gradprop->E;
  double dnu   = gradprop->nu;
  double dcte  = gradprop->W;
  
  double dvd = 0.0;
  double dat = 0.0;
  
  dvd += dx[0]*dx[0] + dy[0]*dy[0] + dz[0]*dz[0];
  dvd += dx[1]*dx[1] + dy[1]*dy[1] + dz[1]*dz[1];
  dvd += dx[2]*dx[2] + dy[2]*dy[2] + dz[2]*dz[2];
  
  dat += dh[0]*dh[0] + dh[1]*dh[1] + dh[2]*dh[2];
  dat += drho*drho;
  dat += demod*demod;
  dat += dnu*dnu;
  dat += dcte*dcte;
  
  if ( dvd > 0 && dat >0 ) return Optvar::mixedvar;
  if ( dvd > 0 )           return Optvar::coordinate;
  if ( dat > 0 )           return Optvar::attribute;
  
  return 0;
}


extern "C"      {
void _FORTRAN(gxtria3d)(int&, double* ,double* ,double* ,double* ,double* ,double* ,
                        double& ,double& ,double& ,
                        double* ,double* ,double* ,double*, double& ,double&);

void _FORTRAN(gxmass8)(double* ,double* ,double* ,double* ,double* ,double* ,
                       double* ,double* ,double& ,double& ,
                       double* ,double* ,const int&, double* ,
		       double* ,const int& ,double& ,const int&);

void _FORTRAN(gxsands8)(double*,double*,double*,double*,double*,double*,
                        double&,double&,double&,double*,double*,
			double*,double*, double*,double*,
			int&, int&, int&, int&, int&, int&, double& ,double&,
			double&, double&);

void  _FORTRAN(gxshvmint)(double*, double*, double*, double*, double*, double*, 
                          double&, double&, double&, double*, double*, 
	                  double*, double*, double&, double&, double&, double&,
                          double&, double&, int&,    double& ,double&);
			  
}

//------------------------------------------------------------------------------
void
ThreeNodeShell_opt::getGradVonMises(Vector &dstress, Vector &weight, 
	                        CoordSet &cs, CoordSet &dcs, 
				Vector &elDisp, Vector &elGrad,
				int strInd, int surface, 
                                double* ndTemp, double* dndTemp)
{ 	
        if (prop == NULL || prop->E <= 0) {
	  dstress.zero();
	  return;
        }
	
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

        //determine average nodal temperature difference relative to ambient
	double dt = (ndTemp[0]+ndTemp[1]+ndTemp[2])/3 - prop->Ta;
	double thermalStrain  = (prop->W)*dt;
	
	
	//the sensitivities of w the thermal strain portion taken out has yet to
	//be checked
	double ddt = (dndTemp[0]+dndTemp[1]+dndTemp[2])/3;	
        double dthermalStrain  = (prop->W)*ddt + (gradprop->W)*dt;

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

	// quick-fix
	_FORTRAN(gxsands8)(x,dx,y,dy,z,dz,
			   emod,demod,nu,h,dh,
			   elDisp.data(),elGrad.data(),
			   reinterpret_cast<double*>(elStress),
			   reinterpret_cast<double*>(delStress),
			   strainFlg,maxsze,maxstr,maxgus,elm,
			   surface,alpha,beta,thermalStrain, dthermalStrain);
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

void
ThreeNodeShell_opt::getGradVonMisesInt(CoordSet &cs, CoordSet &dcs, Vector &d,
                                   Vector &dd, double& sigbar, 
                                   double& fac, int areaFlag,
                                   double& vmint, double& dvmint, double& vol,
                                   double& dvol, 
                                   double* ndTemp, double* dndTemp)
{
        if (prop == NULL || prop->E <= 0) {
	  vmint=0; dvmint=0; vol=0; dvol=0;
	  return;
        }

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

	double emod  = prop->E;
	double demod = gradprop->E;
	double nu  = prop->nu;

        _FORTRAN(gxshvmint)(x, y, z, dx, dy, dz, emod, demod, nu, h, dh, 
                           d.data(), dd.data(), vmint, dvmint, vol, dvol, 
                           sigbar, fac, areaFlag,alpha,beta);

}

//------------------------------------------------------------------------------

double
ThreeNodeShell_opt::getGradMass(CoordSet& cs, CoordSet& dcs)
{
  if (prop == 0) { return 0.0; }

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);
        Node &dnd3 = dcs.getNode(nn[2]);
        
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

        double dmass = darea*rho*eh + area*drho*eh + area*rho*deh;

        return dmass;

}

double
ThreeNodeShell_opt::getGradDCmass(CoordSet& cs, CoordSet& dcs, Vector& d, Vector &gd, double& powFac)
{
        if (prop == NULL) return 0.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);
        Node &dnd3 = dcs.getNode(nn[2]);
        
	double rho = prop->rho;
	double eh  = prop->eh;
        
	double drho = gradprop->rho;
	double deh  = gradprop->eh;

	double  r1[3],  r2[3],  r3[3],  v1[3],  v2[3],  v3[3];
	double dr1[3], dr2[3], dr3[3], dv1[3], dv2[3], dv3[3];

        r1[0] = nd1.x + d[0];  r1[1] = nd1.y + d[1];  r1[2] = nd1.z + d[2];
        r2[0] = nd2.x + d[6];  r2[1] = nd2.y + d[7];  r2[2] = nd2.z + d[8];
        r3[0] = nd3.x + d[12]; r3[1] = nd3.y + d[13]; r3[2] = nd3.z + d[14];

        dr1[0] = dnd1.x + gd[0];  dr1[1] = dnd1.y + gd[1]; ; dr1[2] = dnd1.z + gd[2];
        dr2[0] = dnd2.x + gd[6];  dr2[1] = dnd2.y + gd[7]; ; dr2[2] = dnd2.z + gd[8];
        dr3[0] = dnd3.x + gd[12]; dr3[1] = dnd3.y + gd[13];; dr3[2] = dnd3.z + gd[14];

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

        double   scl = pow(0.5,powFac);
	
	double  v3nm = v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2];
	
	double dv3nm = 2.0*(v3[0]*dv3[0] + v3[1]*dv3[1] + v3[2]*dv3[2]);

        double  area = scl*pow(v3nm,(powFac/2.0));

	double darea = scl*(powFac/2.0)*dv3nm*pow(v3nm,(powFac/2.0-1.0));

        double dmass = darea*rho*eh + area*drho*eh  + area*rho*deh;

        return dmass;

}

//------------------------------------------------------------------------------

void
ThreeNodeShell_opt::gradMassMatrix(CoordSet &cs,CoordSet &dcs,
				   FullSquareMatrix &dmel, double mratio)
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
  
  double rho = prop->rho;
  
  const int numdof = 18;
  
  double * mel = new double[numdof*numdof];
  
  int grvflg = 0, masflg = 0;
  
  double *gravityAcceleration = NULL, *grvfor = NULL,totmas = 0.0;
  
  _FORTRAN(gxmass8)(x,dx,y,dy,z,dz,h,dh,rho,drho,
		    mel,dmel.data(),numdof,gravityAcceleration,
		    grvfor,grvflg,totmas,masflg);
  
  delete [] mel;
  lumpMatrix(dmel, mratio);
}

//------------------------------------------------------------------------------

void
ThreeNodeShell_opt::gradstiffness(CoordSet &cs,CoordSet &dcs,FullSquareMatrix &dd,
                              int flg)
{
        if (prop == 0 || prop->E <= 0 ) {
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

        double * d = new double[18*18];

        _FORTRAN(gxtria3d)(flg,x,dx,y,dy,z,dz,emod,demod,nu,h,dh,d,
                           dd.data(),alpha,beta);

        delete [] d;
}


//------------------------------------------------------------------------------

void
ThreeNodeShell_opt::computeDfDsPressure(CoordSet& cs, CoordSet& dcs, Vector& eleDfDs, 
                                    GeomState* gs, int cflg)
{
     double px = 0.0;
     double py = 0.0;
     double pz = 0.0;
     
     double dpx = 0.0;
     double dpy = 0.0;
     double dpz = 0.0;

     double dmx[3],dmy[3],dmz[3];
     double mx[3],my[3],mz[3];

     int i;
     for(i=0; i<3; ++i) {
       mx[i]=0.0;
       my[i]=0.0;
       mz[i]=0.0;
     
       dmx[i]=0.0;
       dmy[i]=0.0;
       dmz[i]=0.0;
     }
     
     // Compute area of shell

     Node &nd1 = cs.getNode(nn[0]);
     Node &nd2 = cs.getNode(nn[1]);
     Node &nd3 = cs.getNode(nn[2]);
    
     //gradnodes = Coordset class type
     
     Node &dnd1 = dcs.getNode(nn[0]);
     Node &dnd2 = dcs.getNode(nn[1]);
     Node &dnd3 = dcs.getNode(nn[2]);

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
    
     double dr1[3], dr2[3], dr3[3], dv1[3], dv2[3], dnormal[3];

     dr1[0] = dnd1.x; dr1[1] = dnd1.y; dr1[2] = dnd1.z;
     dr2[0] = dnd2.x; dr2[1] = dnd2.y; dr2[2] = dnd2.z;
     dr3[0] = dnd3.x; dr3[1] = dnd3.y; dr3[2] = dnd3.z;
    
     dv1[0] = dr3[0] - dr1[0];
     dv1[1] = dr3[1] - dr1[1];
     dv1[2] = dr3[2] - dr1[2];
    
     dv2[0] = dr2[0] - dr1[0];
     dv2[1] = dr2[1] - dr1[1];
     dv2[2] = dr2[2] - dr1[2];

     // Compute normal to shell using vector cross product rule

     dcrossprod(v2, v1, normal, dv2, dv1, dnormal);
    
     double magnitude = sqrt(normal[0]*normal[0] + normal[1]*normal[1] 
     						 + normal[2]*normal[2]);

     double dmagnitude = (normal[0]*dnormal[0] + normal[1]*dnormal[1] 
     			 + normal[2]*dnormal[2])/magnitude;

     double  area = 0.5*magnitude;
     double darea = 0.5*dmagnitude;
     
     double pressureForce  = pressure * area / 3.0;
     double dpressureForce = ( pressure * darea ) / 3.0;
     
     double xn[3][3], dxn[3][3];
       for(i=0; i<3; ++i) {
           xn[0][i] =  r1[i];
           xn[1][i] =  r2[i];
           xn[2][i] =  r3[i];
	  dxn[0][i] = dr1[i];
          dxn[1][i] = dr2[i];
          dxn[2][i] = dr3[i];
         }
    
     if (! gs ){
       normal[0] /= magnitude;
       normal[1] /= magnitude;
       normal[2] /= magnitude;
      dnormal[0] = -dmagnitude/magnitude*normal[0]+dnormal[0]/magnitude;
      dnormal[1] = -dmagnitude/magnitude*normal[1]+dnormal[1]/magnitude;
      dnormal[2] = -dmagnitude/magnitude*normal[2]+dnormal[2]/magnitude;
      
      dpx = dpressureForce*normal[0] + pressureForce*dnormal[0];
      dpy = dpressureForce*normal[1] + pressureForce*dnormal[1];
      dpz = dpressureForce*normal[2] + pressureForce*dnormal[2];
    
       px = pressureForce*normal[0];
       py = pressureForce*normal[1];
       pz = pressureForce*normal[2];
       
       
       if (cflg) 
       {
         int beam, beamnode[3][2];
         beamnode[0][0] = 0;
         beamnode[0][1] = 1;
         beamnode[1][0] = 0;
         beamnode[1][1] = 2;
         beamnode[2][0] = 1;
         beamnode[2][1] = 2;
	 
         
	 
           double length, dx, dy, dz;
	   double dlength, ddx, ddy, ddz;
           int n1, n2;
           
	  for(beam=0; beam<3; ++beam) {
	   n1 = beamnode[beam][0];
           n2 = beamnode[beam][1];
           dx =  xn[n2][0] -  xn[n1][0];
           dy =  xn[n2][1] -  xn[n1][1];
           dz =  xn[n2][2] -  xn[n1][2];
	  ddx = dxn[n2][0] - dxn[n1][0];
          ddy = dxn[n2][1] - dxn[n1][1];
          ddz = dxn[n2][2] - dxn[n1][2];
           length = sqrt(dx*dx + dy*dy + dz*dz);
	  dlength = (ddx*dx + ddy*dy + ddz*dz)/length; 
           for(i=0; i<3; i++ ) {
	   	v1[i] = xn[n2][i] - xn[n1][i];
	       dv1[i] = dxn[n2][i] - dxn[n1][i];
	      }	
	   dnormalize( v1 , dv1);
	   // Local Y-axis as cross between Z and X
           dcrossprod( normal, v1, v2, dnormal , dv1, dv2 );
           dnormalize( v2 , dv2);
           double  lmy = pressureForce*length/12.0;
	   double dlmy = (dpressureForce*length + pressureForce*dlength)/12.0;
	   
	   mx[n1] += (v2[0]*lmy);
           my[n1] += (v2[1]*lmy);
           mz[n1] += (v2[2]*lmy);
           mx[n2] -= (v2[0]*lmy);
           my[n2] -= (v2[1]*lmy);
           mz[n2] -= (v2[2]*lmy);
	  
	   
	   dmx[n1] += (dv2[0]*lmy) + (v2[0]*dlmy);
           dmy[n1] += (dv2[1]*lmy) + (v2[1]*dlmy);
           dmz[n1] += (dv2[2]*lmy) + (v2[2]*dlmy);
           dmx[n2] -= (dv2[0]*lmy) + (v2[0]*dlmy);
           dmy[n2] -= (dv2[1]*lmy) + (v2[1]*dlmy);
           dmz[n2] -= (dv2[2]*lmy) + (v2[2]*dlmy);
       }
    }
         
    
    eleDfDs[0]  = dpx;
    eleDfDs[1]  = dpy;
    eleDfDs[2]  = dpz;
    eleDfDs[3]  = dmx[0];
    eleDfDs[4]  = dmy[0];
    eleDfDs[5]  = dmz[0];

    eleDfDs[6]  = dpx;
    eleDfDs[7]  = dpy;
    eleDfDs[8]  = dpz;
    eleDfDs[9]  = dmx[1];
    eleDfDs[10] = dmy[1];
    eleDfDs[11] = dmz[1];

    eleDfDs[12] = dpx;
    eleDfDs[13] = dpy;
    eleDfDs[14] = dpz;
    eleDfDs[15] = dmx[2];
    eleDfDs[16] = dmy[2];
    eleDfDs[17] = dmz[2];
    
    }
    
    //if local corrdinates are needed for nonlinear analysis
     if (gs){
     //compute centroid
     double xc0[3], dxc0[3];
     double t0[3][3],xl0[3][3];
     double dt0[3][3],dxl0[3][3];
     int nod;
     for( i=0; i<3; ++i ) {
         xc0[i] = (  xn[0][i] +   xn[1][i] +  xn[2][i] )/3.0;
        dxc0[i]= ( dxn[0][i] +  dxn[1][i] + dxn[2][i] )/3.0;
      }
     
     // Compute t0 transformation matrix with x axis along side 1-2 
     for( i=0; i<3; i++ ) {
         t0[0][i] = xn[1][i] - xn[0][i];
        dt0[0][i] = dxn[1][i] - dxn[0][i];
      }
     dnormalize( t0[0], dt0[0] );
     
     
     // local y axis
     for( i=0; i<3; i++ ){
                       t0[1][i] = xn[2][i] - xn[0][i];
		       dt0[1][i] = dxn[2][i] - dxn[0][i];
         }
     dcrossprod( t0[0], t0[1], t0[2], dt0[0], dt0[1], dt0[2]);
     dnormalize( t0[2], dt0[2] );
     
     // local z axis
     dcrossprod( t0[2], t0[0], t0[1], dt0[2], dt0[0], dt0[1] );
     dnormalize( t0[1], dt0[1] );
     
     // Compute local coordinates of undeformed element
     for( nod=0; nod<3; nod++) {
        for( i=0; i<3; i++ ) {
           xl0[nod][i] = t0[i][0]*(xn[nod][0] - xc0[0])
                        +t0[i][1]*(xn[nod][1] - xc0[1])
                        +t0[i][2]*(xn[nod][2] - xc0[2]);
	  dxl0[nod][i] = dt0[i][0]*(xn[nod][0] - xc0[0])
	                +t0[i][0]*(dxn[nod][0] - dxc0[0])
		        +dt0[i][1]*(xn[nod][1] - xc0[1])
		        +t0[i][1]*(dxn[nod][1] - dxc0[1])
		        +dt0[i][2]*(xn[nod][2] - xc0[2])
		        +t0[i][2]*(dxn[nod][2] - dxc0[2]);       
             }
         }
     
     
        eleDfDs[0]  = 0;
	eleDfDs[1]  = 0;
	eleDfDs[2]  = dpressureForce;
	eleDfDs[3]  = cflg*( dpressureForce*( xl0[2][1] - xl0[0][1]
	                                     +xl0[1][1] - xl0[0][1])
		            + pressureForce*( dxl0[2][1] - dxl0[0][1]
	                                     +dxl0[1][1] - dxl0[0][1]))/8;
	eleDfDs[4]  = cflg*( dpressureForce*( xl0[0][0] - xl0[2][0]
	                                     +xl0[0][0] - xl0[1][0])
			    + pressureForce*( dxl0[0][0] - dxl0[2][0]
	                                     +dxl0[0][0] - dxl0[1][0]))/8;
	eleDfDs[5]  = 0;

	eleDfDs[6]  = 0;
	eleDfDs[7]  = 0;
	eleDfDs[8]  = dpressureForce;
	eleDfDs[9]  = cflg*( dpressureForce*( xl0[0][1] - xl0[1][1]
	                                     +xl0[2][1] - xl0[1][1])
			    + pressureForce*( dxl0[0][1] - dxl0[1][1]
	                                     +dxl0[2][1] - dxl0[1][1]))/8;
	eleDfDs[10] = cflg*( dpressureForce*( xl0[1][0] - xl0[0][0]
	                                     +xl0[1][0] - xl0[2][0])
			    + pressureForce*( dxl0[1][0] - dxl0[0][0]
	                                     +dxl0[1][0] - dxl0[2][0]))/8;
	eleDfDs[11] = 0;

	eleDfDs[12] = 0;
	eleDfDs[13] = 0;
	eleDfDs[14] = dpressureForce;
	eleDfDs[15] = cflg*( dpressureForce*( xl0[1][1] - xl0[2][1]
	                                     +xl0[0][1] - xl0[2][1])
			    + pressureForce*( dxl0[1][1] - dxl0[2][1]
	                                     +dxl0[0][1] - dxl0[2][1]))/8;
	eleDfDs[16] = cflg*( dpressureForce*( xl0[2][0] - xl0[1][0]
	                                     +xl0[2][0] - xl0[0][0])
			    + pressureForce*( dxl0[2][0] - dxl0[1][0]
	                                     +dxl0[2][0] - dxl0[0][0]))/8;
	eleDfDs[17] = 0;
	
      }

}

//------------------------------------------------------------------------------

void ThreeNodeShell_opt::computeGradDisp(double gp[2],double *res) 
{

  // this extrapolates gradient of displacement wrt each dof 
  // of the element to matched nodes (p)
  // only one shape function is used for each node  b/c
  // values that shape functions are multiplied by are only non-zero for 
  // the x,y and z (X) associated with the same node

  res[0] = (1.0 - gp[0] - gp[1]);  // dpX/duX1(e)
  res[1] = gp[0];                  // dpX/duX2(e)
  res[2] = gp[1];                  // dpX/duX3(e)
  res[3] = 0.0;                    // dummy space holder

}


//------------------------------------------------------------------------------
double* ThreeNodeShell_opt::getMidPoint(CoordSet &cs)
{ 
  double * midPoint = new double[3];
  ThreeNodeShell_opt::getMidPoint(cs, midPoint);
  return midPoint;
}

//------------------------------------------------------------------------------
void ThreeNodeShell_opt::getMidPoint(CoordSet &cs, double* midPoint)
{ 
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  
  midPoint[0] = ( nd1.x + nd2.x + nd3.x ) / 3.0; 
  midPoint[1] = ( nd1.y + nd2.y + nd3.y ) / 3.0; 
  midPoint[2] = ( nd1.z + nd2.z + nd3.z ) / 3.0; 
  return;
}


#endif
