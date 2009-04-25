#include	<stdio.h>
#include	<stdlib.h>
#include 	<math.h>

#include 	<Element.d/Beam.d/TimoshenkoBeam.h>
#include	<Math.d/FullSquareMatrix.h>
#include        <Corotational.d/BeamCorotator.h>
#include        <Corotational.d/GeomState.h>
#include        <Corotational.d/utilities.h>
#include 	<Utils.d/dofset.h>
#include 	<Utils.d/linkfc.h>


extern "C"      {
void    _FORTRAN(gxmodmstif7)(double*, double&, double&, double&, double&, 
                              double*, double*, 
			      double&, double&, double&, double&, double&, double&,
			      double&, double&, double&, double&, double&, double&, 
			      double&, double&, double*, double*, double*, double*, 
			      double*, double*, int&);
			      
void    _FORTRAN(gxmass7)(int &,   double*, double&,   double&,   double&,   double&,
        		  double*, double*, double*,   double*,   double*,   double*, 
	   		  double*, double*, const int&, double &, const int& );

void    _FORTRAN(gxsands7)(const int&, double&, double&, double&, double&, double*, double*, 
                           double&,    double&, double&, double&, double&, double&, double&, double&, 
			   double&,    double&, double&, double&, double&, double&,    
			   double*,    double*, double*, double*, double*, double*, double*,
                           double*,    double*, 
			   const int&, const int&, const int&, const int&);
}

//------------------------------------------------------------------------------

void
TimoshenkoBeam::getGradIntrnForce(Vector& delForce,
                                  CoordSet& cs, CoordSet& dcs,
			          double *elDisp, double *elGrad, 
				  int forceIndex, double *ndTemps)
{
        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);

        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);

        double  x[2],  y[2],  z[2];
        double dx[2], dy[2], dz[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;

        int maxsze = 1;
        int maxstr = 7;
        int maxgus = 2;
        int elm    = 1;
        int numel  = 1;

        double delStress[2][7];
	
	double      A = prop->A;        double     dA  = gradprop->A;
	double      E = prop->E;        double     dE  = gradprop->E;
        double    Ixx = prop->Ixx;      double    dIxx = gradprop->Ixx;
	double    Iyy = prop->Iyy;	double    dIyy = gradprop->Iyy;
	double    Izz = prop->Izz;	double    dIzz = gradprop->Izz;
	double alphaY = prop->alphaY;	double dalphaY = gradprop->alphaY;
	double alphaZ = prop->alphaZ;	double dalphaZ = gradprop->alphaZ;
	double      c = prop->c;	double      dc = gradprop->c;
	double     nu = prop->nu;	double     dnu = gradprop->nu;
	
        // variation of frame
	int varFrame = gradFrame(cs,dcs);

	double *  eleFrame = (double*) *elemframe;
	double * deleFrame = (double*)*delemframe;

       _FORTRAN(gxsands7)(elm,A,dA,E,dE,eleFrame,deleFrame,
                          Ixx,dIxx,Iyy,dIyy,Izz,dIzz,alphaY,dalphaY,
                          alphaZ,dalphaZ,c,dc,nu,dnu,x,dx,y,dy,z,dz,
			  elDisp,elGrad,(double*)delStress,
			  numel,maxgus,maxstr,maxsze);

// forceIndex    = 0 =  Nodal stresses along longitudinal axis of the beam
//               = 1 =  Nodal  strains along longitudinal axis of the beam
//               = 2 =  Nodal curvatures in local x-y plane
//               = 3 =  Nodal moments in local x-y plane
//               = 4 =  Nodal curvatures in local x-z plane
//               = 5 =  Nodal moments in local x-z plane

        delForce[0] = delStress[0][forceIndex];
        delForce[1] = delStress[1][forceIndex];
}

//------------------------------------------------------------------------------

double
TimoshenkoBeam::getGradMass(CoordSet& cs, CoordSet& dcs)
{
        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);

        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);

        double x[2], y[2], z[2];
        double dx[2], dy[2], dz[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;

        double ax = x[1] - x[0];
        double ay = y[1] - y[0];
        double az = z[1] - z[0];

        double dax = dx[1] - dx[0];
        double day = dy[1] - dy[0];
        double daz = dz[1] - dz[0];

        double  length = sqrt(ax*ax + ay*ay + az*az);
        double dlength = (ax*dax + ay*day + az*daz)/length;

        double dmass =  dlength*(prop->rho)   *(prop->A)
	             +  length*(gradprop->rho)*(prop->A)
		     +  length*(prop->rho)    *(gradprop->A);

	return dmass;
}

//------------------------------------------------------------------------------

void
TimoshenkoBeam::getGradGravityForce(CoordSet& cs, CoordSet& dcs,
                                    double *gravityAcceleration, Vector& dgravityForce)
{
	double dmassPerNode = 0.5*getGradMass(cs,dcs);

        double dfx = dmassPerNode*gravityAcceleration[0];
        double dfy = dmassPerNode*gravityAcceleration[1];
        double dfz = dmassPerNode*gravityAcceleration[2];

        dgravityForce[0]  =  dfx;
        dgravityForce[1]  =  dfy;
        dgravityForce[2]  =  dfz;
        dgravityForce[3]  = 0.0;
        dgravityForce[4]  = 0.0;
        dgravityForce[5]  = 0.0;
        dgravityForce[6]  =  dfx;
        dgravityForce[7]  =  dfy;
        dgravityForce[8]  =  dfz;
        dgravityForce[9]  = 0.0;
        dgravityForce[10] = 0.0;
        dgravityForce[11] = 0.0;
}

//------------------------------------------------------------------------------

void
TimoshenkoBeam::gradMassMatrix(CoordSet &cs, CoordSet &dcs, FullSquareMatrix& dmel)
{
        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);

        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);

        double  x[2],  y[2],  z[2];
        double dx[2], dy[2], dz[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;

        double *gravityAcceleration = 0, *dgrvfor = 0, totmas = 0.0;

  	int elm = 0, grvflg = 0, masflg = 0;
	
        double      A = prop->A;        
	double     dA = gradprop->A;
	double    rho = prop->rho;	
	double   drho = gradprop->rho;

       _FORTRAN(gxmass7)(elm,dmel.data(),A,dA,rho,drho,x,dx,y,dy,z,dz,
                         gravityAcceleration,dgrvfor,grvflg,totmas,masflg);

}

//------------------------------------------------------------------------------

void
TimoshenkoBeam::gradstiffness(CoordSet &cs,  CoordSet& dcs, FullSquareMatrix& d, int flg)
{
        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);

        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);

	double  x[2],  y[2],  z[2];
        double dx[2], dy[2], dz[2];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;

	double      A = prop->A;        double     dA  = gradprop->A;
	double      E = prop->E;        double     dE  = gradprop->E;
        double    Ixx = prop->Ixx;      double    dIxx = gradprop->Ixx;
	double    Iyy = prop->Iyy;	double    dIyy = gradprop->Iyy;
	double    Izz = prop->Izz;	double    dIzz = gradprop->Izz;
	double alphaY = prop->alphaY;	double dalphaY = gradprop->alphaY;
	double alphaZ = prop->alphaZ;	double dalphaZ = gradprop->alphaZ;
	double     C1 = prop->C1;	double     dC1 = gradprop->C1;
	double     nu = prop->nu;	double     dnu = gradprop->nu;

        // variation of frame
	double varFrame = gradFrame(cs,dcs);
	
	//check whether anything is changing
	
	double chg = dx[0]*dx[0] + dx[1]*dx[1]
	           + dy[0]*dy[0] + dy[1]*dy[1]
	           + dz[0]*dz[0] + dz[1]*dz[1]
		   +      dA*dA 
		   +      dE*dE   
                   +    dIxx*dIxx
                   +    dIyy*dIyy
                   +    dIzz*dIzz
                   + dalphaY*dalphaY
                   + dalphaZ*dalphaZ
                   +     dC1*dC1
                   +     dnu*dnu
		   + varFrame;

	if ( chg == 0 ) {
	  d.zero();
	  return;
	}
	
	double *  eleFrame = (double*) *elemframe;
	double * deleFrame = (double*)*delemframe;

	_FORTRAN(gxmodmstif7)(d.data(), A, dA, E, dE,
			      eleFrame, deleFrame,
                              Ixx, dIxx, Iyy, dIyy, Izz, dIzz,
                              alphaY, dalphaY, alphaZ, dalphaZ, C1, dC1, 
			      nu, dnu, x, dx, y, dy, z, dz, flg);
}

//------------------------------------------------------------------------------

void TimoshenkoBeam::updFrame(CoordSet& cs)
{ 
    double len;

    double * e1    = new double[3];
    double * e2    = new double[3];
    double * e3    = new double[3];

    double * xloc  = new double[3];
    double * curOr = new double[3];
    
    EFrame &theFrame = *elemframe;

    double limit = 1.0e-18;

    // get current orientation

    Node &nd1 = cs.getNode(nn[0]);
    Node &nd2 = cs.getNode(nn[1]);

    curOr[0] =  nd2.x-nd1.x;
    curOr[1] =  nd2.y-nd1.y;
    curOr[2] =  nd2.z-nd1.z;

    len = curOr[0]*curOr[0] + curOr[1]*curOr[1] + curOr[2]*curOr[2];
    len = sqrt(len);

    curOr[0] =  curOr[0]/len;
    curOr[1] =  curOr[1]/len;
    curOr[2] =  curOr[2]/len;	   

    // get rotation angle and rotation axis      

    crossprod(iniOr, curOr, e3);

    len = e3[0]*e3[0] + e3[1]*e3[1] + e3[2]*e3[2];
    len = sqrt(len);

    if (len < limit) return;
      
    e3[0] = e3[0]/len;
    e3[1] = e3[1]/len;
    e3[2] = e3[2]/len;

    double xsin = len;
    double xcos = iniOr[0]*curOr[0] + iniOr[1]*curOr[1] + iniOr[2]*curOr[2];
 
    // loop for all frame axis 

    int i;
    double dummy;
     	
    for (i=0;i<3;i++) {

      // projection onto rotation plane

      dummy = e3[0]*theFrame[i][0] 
            + e3[1]*theFrame[i][1] 
            + e3[2]*theFrame[i][2];	    
      
      e1[0] = theFrame[i][0] - dummy * e3[0];
      e1[1] = theFrame[i][1] - dummy * e3[1];
      e1[2] = theFrame[i][2] - dummy * e3[2];

      len = e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2];
      len = sqrt(len);

      // check if rotation axis and frame-axis are parallel     
 
      if (len < limit) continue;

      e1[0] = e1[0]/len;
      e1[1] = e1[1]/len;
      e1[2] = e1[2]/len;

      dummy = e1[0]*theFrame[i][0] 
            + e1[1]*theFrame[i][1] 
            + e1[2]*theFrame[i][2];	
      
      // local coordinates of rotated frame-axis

      xloc[0] = dummy * xcos;
      xloc[1] = dummy * xsin;

      xloc[2] = e3[0]*theFrame[i][0] 
              + e3[1]*theFrame[i][1] 
              + e3[2]*theFrame[i][2];
  
      // check if rotation axis and frame-axis are parallel     

      crossprod(e3, e1, e2);

      len = e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2];
      len = sqrt(len);

      if (len < limit) continue;
      
      e2[0] = e2[0]/len;
      e2[1] = e2[1]/len;
      e2[2] = e2[2]/len;

      // transformation 
      
      theFrame[i][0] = e1[0]*xloc[0] + e2[0]*xloc[1] + e3[0]*xloc[2];
      theFrame[i][1] = e1[1]*xloc[0] + e2[1]*xloc[1] + e3[1]*xloc[2];
      theFrame[i][2] = e1[2]*xloc[0] + e2[2]*xloc[1] + e3[2]*xloc[2];

   }
   
   // update reference frame
   
   iniOr[0] = curOr[0];
   iniOr[1] = curOr[1];
   iniOr[2] = curOr[2];

   // cleanup
   
   delete [] e1;  
   delete [] e2;  
   delete [] e3;  
   delete [] xloc;
   delete [] curOr;

}   

//------------------------------------------------------------------------------

double TimoshenkoBeam::gradFrame(CoordSet& cs, CoordSet& dcs)
{ 
    double len,dlen;

    double * curOr = new double[3];
    double * dcurOr = new double[3];
    
    EFrame &theFrame  = *elemframe;
    EFrame &gradFrame = *delemframe;

    double limit = 1.0e-18;

    // get current orientation

    Node &nd1 = cs.getNode(nn[0]);
    Node &nd2 = cs.getNode(nn[1]);

    Node &dnd1 = dcs.getNode(nn[0]);
    Node &dnd2 = dcs.getNode(nn[1]);

    curOr[0] =  nd2.x-nd1.x;
    curOr[1] =  nd2.y-nd1.y;
    curOr[2] =  nd2.z-nd1.z;

    dcurOr[0] =  dnd2.x-dnd1.x;
    dcurOr[1] =  dnd2.y-dnd1.y;
    dcurOr[2] =  dnd2.z-dnd1.z;

    // check if beam orientation is changing

    len = dcurOr[0]*dcurOr[0] + dcurOr[1]*dcurOr[1] + dcurOr[2]*dcurOr[2];

    if (len < limit) {

      gradFrame[0][0] = gradFrame[0][1] = gradFrame[0][2] = 0.0;
      gradFrame[1][0] = gradFrame[1][1] = gradFrame[1][2] = 0.0;
      gradFrame[2][0] = gradFrame[2][1] = gradFrame[2][2] = 0.0;

      delete [] curOr;
      delete [] dcurOr;
      
      return 0;
    }

    // get memory for local arrays

    double * e1     = new double[3];
    double * e2     = new double[3];
    double * e3     = new double[3];

    double * de1    = new double[3];
    double * de2    = new double[3];
    double * de3    = new double[3];

    double * diniOr = new double[3];
    
    // scale current orientation and its derivative

    len = curOr[0]*curOr[0] + curOr[1]*curOr[1] + curOr[2]*curOr[2];
    len = sqrt(len);

    dlen = curOr[0]*dcurOr[0] + curOr[1]*dcurOr[1] + curOr[2]*dcurOr[2];
    dlen = dlen/len;

    dcurOr[0] =  (dcurOr[0]*len - curOr[0]*dlen)/(len*len);
    dcurOr[1] =  (dcurOr[1]*len - curOr[1]*dlen)/(len*len);
    dcurOr[2] =  (dcurOr[2]*len - curOr[2]*dlen)/(len*len);	   

    curOr[0] =  curOr[0]/len;
    curOr[1] =  curOr[1]/len;
    curOr[2] =  curOr[2]/len;	   

    // get rotation angle and rotation axis xxx     

    diniOr[0] =  0.0;
    diniOr[1] =  0.0;
    diniOr[2] =  0.0;

    dcrossprod(iniOr, curOr, e3, diniOr, dcurOr, de3);

    double dxcos = iniOr[0]*dcurOr[0] + iniOr[1]*dcurOr[1] + iniOr[2]*dcurOr[2];
 
    // loop for all frame axis 

    int i;
    double xloc, dxloc;
     	
    for (i=0;i<3;i++) {

      // projection onto rotation plane

      e1[0] = theFrame[i][0];  de1[0] = 0.0;
      e1[1] = theFrame[i][1];  de1[1] = 0.0;
      e1[2] = theFrame[i][2];  de1[2] = 0.0;

      len = e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2];
      len = sqrt(len);

      e1[0] = e1[0]/len;
      e1[1] = e1[1]/len;
      e1[2] = e1[2]/len;
     
      // local coordinates of rotated frame-axis

      xloc  = len;
      dxloc = len * dxcos;
  
      // check if rotation axis and frame-axis are parallel xxx    

      dcrossprod(e3, e1, e2, de3, de1, de2);

      // transformation 
      
      gradFrame[i][0] = e1[0]*dxloc + de2[0]*xloc;
      gradFrame[i][1] = e1[1]*dxloc + de2[1]*xloc;
      gradFrame[i][2] = e1[2]*dxloc + de2[2]*xloc;

   }
   
   // cleanup
   
   delete [] e1;  
   delete [] e2;  
   delete [] e3;  

   delete [] de1;  
   delete [] de2;  
   delete [] de3;  

   delete [] curOr;
   delete [] dcurOr;
   delete [] diniOr;
   
   return 1;

}   

//------------------------------------------------------------------------------

int TimoshenkoBeam::chkOptInf(CoordSet &dcs)
{
        Node &dnd1 = dcs.getNode(nn[0]);
        Node &dnd2 = dcs.getNode(nn[1]);

        double dx[2], dy[2], dz[2];

        dx[0] = dnd1.x; dy[0] = dnd1.y; dz[0] = dnd1.z;
        dx[1] = dnd2.x; dy[1] = dnd2.y; dz[1] = dnd2.z;

	double     dA  = gradprop->A;
	double     dE  = gradprop->E;
        double    dIxx = gradprop->Ixx;
	double    dIyy = gradprop->Iyy;
	double    dIzz = gradprop->Izz;
	double dalphaY = gradprop->alphaY;
	double dalphaZ = gradprop->alphaZ;
	double     dC1 = gradprop->C1;
	double     dnu = gradprop->nu;
	
	//check whether anything is changing
	
	double chg = dx[0]*dx[0] + dx[1]*dx[1]
	           + dy[0]*dy[0] + dy[1]*dy[1]
	           + dz[0]*dz[0] + dz[1]*dz[1]
		   + (dx[1]-dx[0]) * (dx[1]-dx[0])
		   + (dy[1]-dy[0]) * (dy[1]-dy[0])
		   + (dz[1]-dz[0]) * (dz[1]-dz[0])
		   +      dA*dA 
		   +      dE*dE   
                   +    dIxx*dIxx
                   +    dIyy*dIyy
                   +    dIzz*dIzz
                   + dalphaY*dalphaY
                   + dalphaZ*dalphaZ
                   +     dC1*dC1
                   +     dnu*dnu;

	return ( chg > 0 ) ? 1 : 0;
}
 
