#include <Structopt.d/Element_opt.d/Quad4.d/FourNodeQuad_opt.h>
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Driver_opt.d/GeoSource_opt.h>

// gxgetcmt  - gradient of material tensor
// gxquad4m  - gradient of element stiffness matrix
// gxq4dmas  - gradient of element mass matrix
// gxsands2  - gradient of element stresses and stains

extern "C"      
{
void _FORTRAN(gxgetcmt)(double&,double&,double&,double&,double&,double*,double*);

void _FORTRAN(gxquad4m)(double*,double*,double*,double*,double*,double*,
                        double*,double*,const int&,double*,double*,const int&,double&);

void _FORTRAN(gxq4dmas)(const int&,double&,double&,double*,double*,double*,
                        double*,double*,double*,double*,double*,double*,double*,
                        const int&,double&,double&,const int&, double&);

void _FORTRAN(gxsands2)(char*,double*,double*,double*,double*,double*,double*,
                        double*,double*,double*,double*,double*,double*,int&,
                        int&,int&,int&,int&,int&,double&,double&,double&,
                        double&,double&,double&,double&,double*,double*,double&);

void _FORTRAN(gxquad2d)(double*,double*,double*,double*,double&,double&,const int&,
                        double*,double*,const int&);
}

//------------------------------------------------------------------------------

int FourNodeQuad_opt::chkOptInf(CoordSet &d_cs)
{
  if (gradprop == 0) return 0;

  bool attributes_change =
    (gradprop->rho!=0) ||
    (gradprop->eh!=0) ||
    (gradprop->E!=0) ||
    (gradprop->W!=0) ||
    (gradprop->nu!=0);
  
  const Node &d_nd1 = d_cs.getNode(nn[0]);
  const Node &d_nd2 = d_cs.getNode(nn[1]);
  const Node &d_nd3 = d_cs.getNode(nn[2]);
  const Node &d_nd4 = d_cs.getNode(nn[3]);
  
  bool coordinates_change =
    (d_nd1.x!=0) || (d_nd1.y!=0) ||
    (d_nd2.x!=0) || (d_nd2.y!=0) ||
    (d_nd3.x!=0) || (d_nd3.y!=0) ||
    (d_nd4.x!=0) || (d_nd4.y!=0);
  
  if(coordinates_change && attributes_change) return Optvar::mixedvar;
  if(coordinates_change)	              return Optvar::coordinate;
  if(attributes_change)	                      return Optvar::attribute;  
  return 0;

}

//------------------------------------------------------------------------------

void FourNodeQuad_opt::gradstiffness(CoordSet &cs,CoordSet &d_cs,
                                 FullSquareMatrix &dd, int flg)
{

  if (prop == 0) {
     dd.zero();
     return;
  }

  const int numgauss = 2;
  const int numdof   = 8;

  static double d_x[4]={0,0,0,0}, d_y[4]={0,0,0,0};
  double d_h[4];
  double c[16],d_c[16];

  if(static_cast<GeoSource_opt*>(geoSource)->areThereMovingNodes())
    {
      // gradient orientated data
      Node &d_nd1 = d_cs.getNode(nn[0]);
      Node &d_nd2 = d_cs.getNode(nn[1]);
      Node &d_nd3 = d_cs.getNode(nn[2]);
      Node &d_nd4 = d_cs.getNode(nn[3]);
      
      d_x[0] = d_nd1.x; d_y[0] = d_nd1.y;
      d_x[1] = d_nd2.x; d_y[1] = d_nd2.y;
      d_x[2] = d_nd3.x; d_y[2] = d_nd3.y;
      d_x[3] = d_nd4.x; d_y[3] = d_nd4.y;      
    }

  d_h[0] = d_h[1] = d_h[2] = d_h[3] = gradprop->eh;

  _FORTRAN(gxgetcmt)(prop->A,prop->E,gradprop->E,prop->nu,gradprop->nu,c,d_c);

  // non-gradient orientated data
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);

  double x[4],y[4],h[4];

  x[0] = nd1.x; y[0] = nd1.y;
  x[1] = nd2.x; y[1] = nd2.y;
  x[2] = nd3.x; y[2] = nd3.y;
  x[3] = nd4.x; y[3] = nd4.y;

  h[0]  = h[1]  = h[2]  = h[3]  = prop->eh;

  double *d = new double[8*8];

  _FORTRAN(gxquad4m)(x,d_x,y,d_y,h,d_h,c,d_c,numgauss,d,dd.data(),
                     numdof,prop->A);

  delete [] d;

}

//-----------------------------------------------------------------------

double FourNodeQuad_opt::getGradMass(CoordSet& cs, CoordSet& d_cs)
{

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);
	
  double x[4],y[4];
  x[0] = nd1.x; y[0] = nd1.y;
  x[1] = nd2.x; y[1] = nd2.y;
  x[2] = nd3.x; y[2] = nd3.y;
  x[3] = nd4.x; y[3] = nd4.y;

  double h     = prop->eh;
  double rho   = prop->rho;

  double d_h   = gradprop->eh;
  double d_rho = gradprop->rho;

  double area,d_area;
  double volume,d_volume;
  
  // element area: 0.5 * (det of coord.)
  area     = 0.5 * (x[0]*y[1] - x[1]*y[0] + x[1]*y[2] - x[2]*y[1] -
                    x[0]*y[3] + x[3]*y[0] + x[2]*y[3] - x[3]*y[2]   );  

  volume   = area * h;


  Node &d_nd1 = d_cs.getNode(nn[0]);
  Node &d_nd2 = d_cs.getNode(nn[1]);
  Node &d_nd3 = d_cs.getNode(nn[2]);
  Node &d_nd4 = d_cs.getNode(nn[3]);
  
  double d_x[4]={0,0,0,0}, d_y[4]={0,0,0,0};
  d_x[0] = d_nd1.x; d_y[0] = d_nd1.y; 
  d_x[1] = d_nd2.x; d_y[1] = d_nd2.y; 
  d_x[2] = d_nd3.x; d_y[2] = d_nd3.y; 
  d_x[3] = d_nd4.x; d_y[3] = d_nd4.y; 
  
  d_area   = 0.5 * (d_x[0]*y[1] - d_x[1]*y[0] + d_x[1]*y[2] - d_x[2]*y[1] -
		    d_x[0]*y[3] + d_x[3]*y[0] + d_x[2]*y[3] - d_x[3]*y[2] + 
		    x[0]*d_y[1] - x[1]*d_y[0] + x[1]*d_y[2] - x[2]*d_y[1] -
		    x[0]*d_y[3] + x[3]*d_y[0] + x[2]*d_y[3] - x[3]*d_y[2]    );

  d_volume = d_h*area + h*d_area;
  double dtotmas = d_rho*volume + rho*d_volume;

  return dtotmas;
}

//------------------------------------------------------------------------------

void FourNodeQuad_opt::gradMassMatrix(CoordSet &cs,CoordSet &d_cs,
				      FullSquareMatrix &d_mel, double mratio)
{

  if (prop == 0) {
     d_mel.zero();
     return;
  }
	
  double x[4],y[4],h[4];

  h[0] = h[1] = h[2] = h[3] = prop->eh;

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);

  x[0] = nd1.x; y[0] = nd1.y;
  x[1] = nd2.x; y[1] = nd2.y;
  x[2] = nd3.x; y[2] = nd3.y;
  x[3] = nd4.x; y[3] = nd4.y;

  static double d_x[4]={0,0,0,0},d_y[4]={0,0,0,0};
  double d_h[4];

  if(static_cast<GeoSource_opt*>(geoSource)->areThereMovingNodes())
    {
      Node &d_nd1 = d_cs.getNode(nn[0]);
      Node &d_nd2 = d_cs.getNode(nn[1]);
      Node &d_nd3 = d_cs.getNode(nn[2]);
      Node &d_nd4 = d_cs.getNode(nn[3]);
      
      d_x[0] = d_nd1.x; d_y[0] = d_nd1.y;
      d_x[1] = d_nd2.x; d_y[1] = d_nd2.y;
      d_x[2] = d_nd3.x; d_y[2] = d_nd3.y;
      d_x[3] = d_nd4.x; d_y[3] = d_nd4.y;
    }

  d_h[0] = d_h[1] = d_h[2] = d_h[3] = gradprop->eh;

  const int numgauss = 2;

  int       grvflg   = 0,masflg = 0;
  double    mmratio = 1.0-mratio;

  double    totmas   = 0.0;
  double    d_totmas = 0.0;

  double   *gravityAcceleration = 0,*grvfor = 0,*d_grvfor = 0;
 
  _FORTRAN(gxq4dmas)(numgauss,prop->rho,gradprop->rho,d_mel.data(),h,
                     d_h,x,d_x,y,d_y,gravityAcceleration,grvfor,d_grvfor,grvflg,
                     totmas,d_totmas,masflg,mmratio);

}

//------------------------------------------------------------------------------

// FIX - nodal temps aren't implemented right now b/c not input into function

void FourNodeQuad_opt::getGradVonMises(Vector& d_stress,Vector& weight,CoordSet &cs,
                                   CoordSet &d_cs,Vector& elDisp,Vector &elGrad,
                                   int strInd,int surface, 
                                   double* ndTemp, double* dndTemp) 
{

  weight = 1.0;

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);

  double x[4],y[4],c[16];
  double d_c[16];

  x[0] = nd1.x; y[0] = nd1.y;
  x[1] = nd2.x; y[1] = nd2.y;
  x[2] = nd3.x; y[2] = nd3.y;
  x[3] = nd4.x; y[3] = nd4.y;

  static double d_x[4]={0,0,0,0},d_y[4]={0,0,0,0};
  if(static_cast<GeoSource_opt*>(geoSource)->areThereMovingNodes())
    {
      Node &d_nd1 = d_cs.getNode(nn[0]);
      Node &d_nd2 = d_cs.getNode(nn[1]);
      Node &d_nd3 = d_cs.getNode(nn[2]);
      Node &d_nd4 = d_cs.getNode(nn[3]);

      d_x[0] = d_nd1.x; d_y[0] = d_nd1.y;
      d_x[1] = d_nd2.x; d_y[1] = d_nd2.y;
      d_x[2] = d_nd3.x; d_y[2] = d_nd3.y;
      d_x[3] = d_nd4.x; d_y[3] = d_nd4.y;
    }
  double d_E  = gradprop->E;
  double d_nu = gradprop->nu;
  double d_al = gradprop->W;

  double E    = prop->E;
  double nu   = prop->nu;
  double al   = prop->W;
  double Tref = prop->Ta;

  _FORTRAN(gxgetcmt)(prop->A,E,d_E,nu,d_nu,c,d_c);

  int maxgus = 4; // maximum gauss points 
  int maxstr = 7; // maximum  
  int elm    = 1; // element number
  int maxsze = 1;

  int vmflg     = 0;
  int strainFlg = 0;

  double elStress[4][7];  
  double elStrain[4][7];  
  double d_elStress[4][7]; 
  double d_elStrain[4][7]; 

  char ESCM[7] = "EXTRAP";   // ... STRESS EXTRAPOLATION FROM GAUSS POINTS

  // sands2 to calculate stresses.

  if (strInd <= 6) {

     // Flags sands2 to calculate Von Mises stress.
     if (strInd == 6) vmflg  = 1;
   
     _FORTRAN(gxsands2)(ESCM,x,d_x,y,d_y,c,d_c,elDisp.data(),elGrad.data(),
                        reinterpret_cast<double*>(elStress),
			reinterpret_cast<double*>(d_elStress),
                        reinterpret_cast<double*>(elStrain),
			reinterpret_cast<double*>(d_elStrain),
			maxgus,maxstr,
                        elm,maxsze,vmflg,strainFlg,E,d_E,nu,d_nu,al,d_al,
			Tref,ndTemp,dndTemp,prop->A);
  			   
     d_stress[0] = d_elStress[0][strInd];
     d_stress[1] = d_elStress[1][strInd];
     d_stress[2] = d_elStress[2][strInd];
     d_stress[3] = d_elStress[3][strInd];
  } 
  else {
    // sands2 to calculate strains.

    // Flags sands2 to calculate Von Mises strain
    if(strInd == 13) strainFlg = 1;

     _FORTRAN(gxsands2)(ESCM,x,d_x,y,d_y,c,d_c,elDisp.data(),elGrad.data(),
                        reinterpret_cast<double*>(elStress),
			reinterpret_cast<double*>(d_elStress),
                        reinterpret_cast<double*>(elStrain),
			reinterpret_cast<double*>(d_elStrain),maxgus,maxstr,
                        elm,maxsze,vmflg,strainFlg,E,d_E,nu,d_nu,al,d_al,
			Tref,ndTemp,dndTemp,prop->A);
  		 
    d_stress[0] = d_elStrain[0][strInd-7];
    d_stress[1] = d_elStrain[1][strInd-7];
    d_stress[2] = d_elStrain[2][strInd-7];
    d_stress[3] = d_elStrain[3][strInd-7];
  }
}

//------------------------------------------------------------------------------

void FourNodeQuad_opt::computeGradDisp(double gp[2],double *res) 
{

  // this extrapolates gradient of displacement wrt each dof 
  // of the element to matched nodes (p) 
  // only one shape function is used for each node  b/c
  // values that shape functions are multiplied by are only non-zero for 
  // the x,y and z (X) associated with the same node
  // 4 entries (12 dofs, 3 the same for each node) must be filled up,
  // even if not that many for the element

  res[0] = (1-gp[0]) * (1-gp[1]);   // dpX/duX1(e)
  res[1] =    gp[0]  * (1-gp[1]);   // dpX/duX2(e)
  res[2] =    gp[0]  *	  gp[1];    // dpX/duX3(e)
  res[3] = (1-gp[0]) *	  gp[1];    // dpX/duX4(e)

}

//------------------------------------------------------------------------------
void FourNodeQuad_opt::getMidPoint(CoordSet &cs, double* midPoint)
{
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);
  
  midPoint[0] = ( nd1.x + nd2.x + nd3.x + nd4.x) / 4.0;
  midPoint[1] = ( nd1.y + nd2.y + nd3.y + nd4.y) / 4.0;
  midPoint[2] = ( nd1.z + nd2.z + nd3.z + nd4.z) / 4.0;
  return;
}

//------------------------------------------------------------------------------
double *
FourNodeQuad_opt::getMidPoint(CoordSet &cs)
{ 
  double * midPoint = new double[3];
  FourNodeQuad_opt::getMidPoint(cs, midPoint);
  return midPoint;
}
