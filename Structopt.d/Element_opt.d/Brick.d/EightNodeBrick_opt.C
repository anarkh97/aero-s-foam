#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/Brick.d/EightNodeBrick_opt.h>
#include <Structopt.d/Optvar.h>

extern "C"      
{
void  _FORTRAN(gxbrkcmt)(double&, double&, double&, double&, double*, double*);

void  _FORTRAN(gxbrik8v)(double*, double*, double*, double*, double*, 
			 double*, double*, double*, const int&,    
			 double*, int&);

void  _FORTRAN(gxbr8mas)(const int&p,double&rho,double&g_rho,double*g_elmass,
			 double*x,double*y,double*z,
			 double*g_x,double*g_y,double*g_z,double*gamma,
			 double*grvfor,double*g_grvfor,const int&grvflag, 
			 double &g_totmas, const int&masflg, double& mratio);
  
void  _FORTRAN(gxsands17)(const int&,double*,double*,double*,double*,double*,
                        double*,double*,double*,double*,double*,
                        double*,double*,const int&,const int&,
			const int&,const int&,const int&,const int&);       

void  _FORTRAN(gxvol17)(const int&,double*,double*,double*,
                         double*,double*,double*,double &,double &);

void  _FORTRAN(gxbrvmint)(double*, double*, double*,double*, double*, double*, 
                          double*, double*, double*, double*, const int&,  
	                  double&, double&, double&, double&, double&, 
                          double&, int&,    int&);
}

//-------------------------------------------------------------------
double* EightNodeBrick_opt::getMidPoint(CoordSet &cs)
{ 
  double * midPoint = new double[3];
  EightNodeBrick_opt::getMidPoint(cs, midPoint);
  return midPoint;
}

//-------------------------------------------------------------------
void EightNodeBrick_opt::getMidPoint(CoordSet &cs, double* midPoint)
{ 
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);
  Node &nd5 = cs.getNode(nn[4]);
  Node &nd6 = cs.getNode(nn[5]);
  Node &nd7 = cs.getNode(nn[6]);
  Node &nd8 = cs.getNode(nn[7]);
  
  midPoint[0] = ( nd1.x + nd2.x + nd3.x + nd4.x 
		  + nd5.x + nd6.x + nd7.x + nd8.x ) / 8.0; 
  midPoint[1] = ( nd1.y + nd2.y + nd3.y + nd4.y 
		  + nd5.y + nd6.y + nd7.y + nd8.y ) / 8.0; 
  midPoint[2] = ( nd1.z + nd2.z + nd3.z + nd4.z 
		  + nd5.z + nd6.z + nd7.z + nd8.z ) / 8.0; 
  return;
}

//------------------------------------------------------------------------------
int EightNodeBrick_opt::chkOptInf(CoordSet &dcs)
{
  if (prop == 0) return 0;

  Node &dnd1 = dcs.getNode(nn[0]);
  Node &dnd2 = dcs.getNode(nn[1]);
  Node &dnd3 = dcs.getNode(nn[2]);
  Node &dnd4 = dcs.getNode(nn[3]);
  Node &dnd5 = dcs.getNode(nn[4]);
  Node &dnd6 = dcs.getNode(nn[5]);
  Node &dnd7 = dcs.getNode(nn[6]);
  Node &dnd8 = dcs.getNode(nn[7]);

  double dx[8], dy[8], dz[8];

  dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
  dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
  dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;
  dx[3] = dnd4.x; dy[3] = dnd4.y; dz[3] = dnd4.z;
  dx[4] = dnd5.x; dy[4] = dnd5.y; dz[4] = dnd5.z;
  dx[5] = dnd6.x; dy[5] = dnd6.y; dz[5] = dnd6.z;
  dx[6] = dnd7.x; dy[6] = dnd7.y; dz[6] = dnd7.z;
  dx[7] = dnd8.x; dy[7] = dnd8.y; dz[7] = dnd8.z;
	
  double drho = gradprop->rho;
  double demod = gradprop->E;
  double dnu   = gradprop->nu;

  double dvd = 0.0;
  double dat = 0.0;
	
  dvd += dx[0]*dx[0] + dy[0]*dy[0] + dz[0]*dz[0];
  dvd += dx[1]*dx[1] + dy[1]*dy[1] + dz[1]*dz[1];
  dvd += dx[2]*dx[2] + dy[2]*dy[2] + dz[2]*dz[2];
  dvd += dx[3]*dx[3] + dy[3]*dy[3] + dz[3]*dz[3];
  dvd += dx[4]*dx[4] + dy[4]*dy[4] + dz[4]*dz[4];
  dvd += dx[5]*dx[5] + dy[5]*dy[5] + dz[5]*dz[5];
  dvd += dx[6]*dx[6] + dy[6]*dy[6] + dz[6]*dz[6];
  dvd += dx[7]*dx[7] + dy[7]*dy[7] + dz[7]*dz[7];
		
  dat += drho*drho;
  dat += demod*demod;
  dat += dnu*dnu;

  if ( dvd > 0 && dat >0 ) return Optvar::mixedvar;
  if ( dvd > 0 )           return Optvar::coordinate;
  if ( dat > 0 )           return Optvar::attribute;

  return 0;
}

//------------------------------------------------------------------------------
void EightNodeBrick_opt::gradstiffness(CoordSet &cs, CoordSet &dcs, FullSquareMatrix &dd, int flg)
{
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

  Node &dnd1 = dcs.getNode(nn[0]);
  Node &dnd2 = dcs.getNode(nn[1]);
  Node &dnd3 = dcs.getNode(nn[2]);
  Node &dnd4 = dcs.getNode(nn[3]);
  Node &dnd5 = dcs.getNode(nn[4]);
  Node &dnd6 = dcs.getNode(nn[5]);
  Node &dnd7 = dcs.getNode(nn[6]);
  Node &dnd8 = dcs.getNode(nn[7]);

  double dx[8], dy[8], dz[8], dc[6][6];

  dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
  dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
  dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;
  dx[3] = dnd4.x; dy[3] = dnd4.y; dz[3] = dnd4.z;
  dx[4] = dnd5.x; dy[4] = dnd5.y; dz[4] = dnd5.z;
  dx[5] = dnd6.x; dy[5] = dnd6.y; dz[5] = dnd6.z;
  dx[6] = dnd7.x; dy[6] = dnd7.y; dz[6] = dnd7.z;
  dx[7] = dnd8.x; dy[7] = dnd8.y; dz[7] = dnd8.z;
	
  _FORTRAN(gxbrkcmt)(prop->E, gradprop->E, prop->nu,
		     gradprop->nu,
		     reinterpret_cast<double*>(c), 
		     reinterpret_cast<double*>(dc));

  const int numgauss = 2;

  int status;

  _FORTRAN(gxbrik8v)(x, dx, y, dy, z, dz, reinterpret_cast<double *>(c), 
		     reinterpret_cast<double *>(dc), numgauss, dd.data(), status);
}

//------------------------------------------------------------------------------
void EightNodeBrick_opt::gradMassMatrix(CoordSet &cs, CoordSet &dcs, FullSquareMatrix &dmel, double mratio)
{
  if (prop == 0) {
     dmel.zero();
     return;
  }
	
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);
  Node &nd5 = cs.getNode(nn[4]);
  Node &nd6 = cs.getNode(nn[5]);
  Node &nd7 = cs.getNode(nn[6]);
  Node &nd8 = cs.getNode(nn[7]);

  double x[8], y[8], z[8];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
  x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
  x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
  x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
  x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;

  Node &dnd1 = dcs.getNode(nn[0]);
  Node &dnd2 = dcs.getNode(nn[1]);
  Node &dnd3 = dcs.getNode(nn[2]);
  Node &dnd4 = dcs.getNode(nn[3]);
  Node &dnd5 = dcs.getNode(nn[4]);
  Node &dnd6 = dcs.getNode(nn[5]);
  Node &dnd7 = dcs.getNode(nn[6]);
  Node &dnd8 = dcs.getNode(nn[7]);

  double dx[8], dy[8], dz[8];

  dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
  dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
  dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;
  dx[3] = dnd4.x; dy[3] = dnd4.y; dz[3] = dnd4.z;
  dx[4] = dnd5.x; dy[4] = dnd5.y; dz[4] = dnd5.z;
  dx[5] = dnd6.x; dy[5] = dnd6.y; dz[5] = dnd6.z;
  dx[6] = dnd7.x; dy[6] = dnd7.y; dz[6] = dnd7.z;
  dx[7] = dnd8.x; dy[7] = dnd8.y; dz[7] = dnd8.z;

  double *gravityAcceleration = 0, *grvfor = 0, *g_grvfor=0;

  int grvflg = 0, masflg=0;

  double dtotmas = 0.0;
 
  const int numgauss = 2;
  mratio=1.0-mratio;
  _FORTRAN(gxbr8mas)(numgauss,prop->rho,gradprop->rho,dmel.data(),
		     x,y,z,dx,dy,dz,gravityAcceleration,
		     grvfor,g_grvfor,grvflg,dtotmas,masflg,mratio);
}

//-----------------------------------------------------------------------

void
EightNodeBrick_opt::getGradVonMises(Vector& dstress, Vector& weight, CoordSet &cs,
                CoordSet &dcs, Vector& elDisp, Vector &elGrad, int strInd,
		int surface, double* ndTemp, double*  dndTemp)
{
	weight = 1.0;

	Node &nd1 = cs.getNode(nn[0]);
	Node &nd2 = cs.getNode(nn[1]);
	Node &nd3 = cs.getNode(nn[2]);
	Node &nd4 = cs.getNode(nn[3]);
  	Node &nd5 = cs.getNode(nn[4]);
  	Node &nd6 = cs.getNode(nn[5]);
  	Node &nd7 = cs.getNode(nn[6]);
  	Node &nd8 = cs.getNode(nn[7]);
	
	Node &dnd1 = dcs.getNode(nn[0]);
	Node &dnd2 = dcs.getNode(nn[1]);
	Node &dnd3 = dcs.getNode(nn[2]);
	Node &dnd4 = dcs.getNode(nn[3]);
  	Node &dnd5 = dcs.getNode(nn[4]);
  	Node &dnd6 = dcs.getNode(nn[5]);
  	Node &dnd7 = dcs.getNode(nn[6]);
  	Node &dnd8 = dcs.getNode(nn[7]);

  	double x[8], y[8], z[8], c[6][6];
  	double dx[8], dy[8], dz[8], dc[6][6];

  	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
  	x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
  	x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
  	x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
  	x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;

  	dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
  	dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
  	dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;
  	dx[3] = dnd4.x; dy[3] = dnd4.y; dz[3] = dnd4.z;
  	dx[4] = dnd5.x; dy[4] = dnd5.y; dz[4] = dnd5.z;
  	dx[5] = dnd6.x; dy[5] = dnd6.y; dz[5] = dnd6.z;
  	dx[6] = dnd7.x; dy[6] = dnd7.y; dz[6] = dnd7.z;
  	dx[7] = dnd8.x; dy[7] = dnd8.y; dz[7] = dnd8.z;
	
	_FORTRAN(gxbrkcmt)(prop->E, gradprop->E, prop->nu, gradprop->nu,
			   reinterpret_cast<double*>(c), 
			   reinterpret_cast<double*>(dc));

        int maxgus = 8; // maximum gauss points 
        int maxstr = 7; // maximum  
        int elm    = 1; // element number
        int maxsze = 1;

	int vmflg     = 0;
        int strainFlg = 0;
	
	// sands17 to calculate stresses.
        if(strInd <= 6){
	      // Flags sands17 to calculate Von Mises stress.
	      if(strInd == 6) vmflg  = 1;
	 
              double elStress[8][7];  
              double delStress[8][7]; 
	
              _FORTRAN(gxsands17)(elm,x,y,z,dx,dy,dz,
				  reinterpret_cast<double*>(c),
				  reinterpret_cast<double*>(dc),
				  elDisp.data(),elGrad.data(),
				  reinterpret_cast<double*>(elStress),
				  reinterpret_cast<double*>(delStress),
				  maxgus,maxstr,maxsze,vmflg,
				  strainFlg,strInd);
	     		 
              dstress[0] = delStress[0][strInd];
              dstress[1] = delStress[1][strInd];
              dstress[2] = delStress[2][strInd];
              dstress[3] = delStress[3][strInd];
              dstress[4] = delStress[4][strInd];
              dstress[5] = delStress[5][strInd];
              dstress[6] = delStress[6][strInd];
              dstress[7] = delStress[7][strInd];
	} 
	else {
	      // sands17 to calculate strains.
	
	      // Flags sands17 to calculate Von Mises strain
	      if(strInd == 13) strainFlg = 1;

  	      double elStrain[8][7]; 
  	      double delStrain[8][7];  
	
              _FORTRAN(gxsands17)(elm,x,y,z,dx,dy,dz,
				  reinterpret_cast<double*>(c),
				  reinterpret_cast<double*>(dc),
				  elDisp.data(),elGrad.data(),
				  reinterpret_cast<double*>(elStrain),
				  reinterpret_cast<double*>(delStrain),
				  maxgus,maxstr,maxsze,vmflg,
				  strainFlg,strInd); 
	        	   
	      dstress[0] = delStrain[0][strInd-7];
	      dstress[1] = delStrain[1][strInd-7];
	      dstress[2] = delStrain[2][strInd-7];
	      dstress[3] = delStrain[3][strInd-7];
	      dstress[4] = delStrain[4][strInd-7];
              dstress[5] = delStrain[5][strInd-7];
              dstress[6] = delStrain[6][strInd-7];
              dstress[7] = delStrain[7][strInd-7];
        }
	
}

//-----------------------------------------------------------------------

void
EightNodeBrick_opt::getGradVonMisesInt(CoordSet &cs, CoordSet &dcs, Vector &d, Vector &dd, 
                            double& sigbar, double& fac, int areaFlag,
                            double& vmint, double& dvmint, double& vol, double& dvol, 
                            double* ndTemp, double*  dndTemp)
{
	Node &nd1 = cs.getNode(nn[0]);
	Node &nd2 = cs.getNode(nn[1]);
	Node &nd3 = cs.getNode(nn[2]);
	Node &nd4 = cs.getNode(nn[3]);
	Node &nd5 = cs.getNode(nn[4]);
	Node &nd6 = cs.getNode(nn[5]);
	Node &nd7 = cs.getNode(nn[6]);
	Node &nd8 = cs.getNode(nn[7]);
	
	Node &dnd1 = dcs.getNode(nn[0]);
	Node &dnd2 = dcs.getNode(nn[1]);
	Node &dnd3 = dcs.getNode(nn[2]);
	Node &dnd4 = dcs.getNode(nn[3]);
  	Node &dnd5 = dcs.getNode(nn[4]);
  	Node &dnd6 = dcs.getNode(nn[5]);
  	Node &dnd7 = dcs.getNode(nn[6]);
  	Node &dnd8 = dcs.getNode(nn[7]);

  	double x[8], y[8], z[8], c[6][6];
  	double dx[8], dy[8], dz[8], dc[6][6];

  	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
  	x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
  	x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
  	x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
  	x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;

  	dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
  	dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
  	dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;
  	dx[3] = dnd4.x; dy[3] = dnd4.y; dz[3] = dnd4.z;
  	dx[4] = dnd5.x; dy[4] = dnd5.y; dz[4] = dnd5.z;
  	dx[5] = dnd6.x; dy[5] = dnd6.y; dz[5] = dnd6.z;
  	dx[6] = dnd7.x; dy[6] = dnd7.y; dz[6] = dnd7.z;
  	dx[7] = dnd8.x; dy[7] = dnd8.y; dz[7] = dnd8.z;
	
	_FORTRAN(gxbrkcmt)(prop->E, gradprop->E, prop->nu, gradprop->nu, 
			   reinterpret_cast<double*>(c), 
			   reinterpret_cast<double*>(dc));

  	const int numgauss = 2;
	int status;

        _FORTRAN(gxbrvmint)(x, y, z, dx, dy, dz, d.data(), dd.data(),
			    reinterpret_cast<double *>(c), 
			    reinterpret_cast<double *>(dc), 
			    numgauss, vmint, 
                            dvmint, vol, dvol, sigbar, fac, areaFlag, status);

}

//-----------------------------------------------------------------------

double
EightNodeBrick_opt::getGradMass(CoordSet& cs, CoordSet& dcs)
{
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);
  Node &nd5 = cs.getNode(nn[4]);
  Node &nd6 = cs.getNode(nn[5]);
  Node &nd7 = cs.getNode(nn[6]);
  Node &nd8 = cs.getNode(nn[7]);
	
  Node &dnd1 = dcs.getNode(nn[0]);
  Node &dnd2 = dcs.getNode(nn[1]);
  Node &dnd3 = dcs.getNode(nn[2]);
  Node &dnd4 = dcs.getNode(nn[3]);
  Node &dnd5 = dcs.getNode(nn[4]);
  Node &dnd6 = dcs.getNode(nn[5]);
  Node &dnd7 = dcs.getNode(nn[6]);
  Node &dnd8 = dcs.getNode(nn[7]);

  double x[8], y[8], z[8];
  double dx[8], dy[8], dz[8];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
  x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
  x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
  x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
  x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;

  dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
  dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;
  dx[2] = dnd3.x; dy[2] = dnd3.y; dz[2] = dnd3.z;
  dx[3] = dnd4.x; dy[3] = dnd4.y; dz[3] = dnd4.z;
  dx[4] = dnd5.x; dy[4] = dnd5.y; dz[4] = dnd5.z;
  dx[5] = dnd6.x; dy[5] = dnd6.y; dz[5] = dnd6.z;
  dx[6] = dnd7.x; dy[6] = dnd7.y; dz[6] = dnd7.z;
  dx[7] = dnd8.x; dy[7] = dnd8.y; dz[7] = dnd8.z;

  const int numgauss = 2;

  double volume =0.0;
  double dvolume=0.0;
  
  _FORTRAN(gxvol17)(numgauss,x,y,z,dx,dy,dz,volume,dvolume);
  
  double dtotmas = gradprop->rho*volume + prop->rho*dvolume;

  return dtotmas;

}

#endif
