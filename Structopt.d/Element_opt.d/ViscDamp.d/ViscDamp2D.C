#ifdef STRUCTOPT

#include <cstdio>
#include <cmath>
#include <algorithm>

#include <Utils.d/dofset.h>
#include <Structopt.d/Element_opt.d/ViscDamp.d/ViscDamp2D.h>
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Driver_opt.d/GeoSource_opt.h>

//----------------------------------------------------------------------
// class constants
const int    ViscDamp2D::nNodes = ViscDamp2D_nNodes;
const int    ViscDamp2D::nDim   = 2;

//----------------------------------------------------------------------
ViscDamp2D::ViscDamp2D(int* nodenums)
{
  std::copy(nodenums, nodenums+nNodes, nn);
}

//----------------------------------------------------------------------
int* ViscDamp2D::nodes(int* p)
{
  if(p == 0) 
    {
      p = new int[nNodes];
    }
  std::copy(nn, nn+nNodes, p);
  return p;
}

//----------------------------------------------------------------------
int* ViscDamp2D::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) 
    {
      p = new int[numDofs()];
    }
  dsa.number(nn[0], DofSet::Xdisp | DofSet::Ydisp, p  );
  dsa.number(nn[1], DofSet::Xdisp | DofSet::Ydisp, p+dim());
  return p;
}

//----------------------------------------------------------------------
void ViscDamp2D::markDofs( DofSetArray &dsa )
{
  dsa.mark(nn, nNodes, DofSet::Xdisp | DofSet::Ydisp);
}

//----------------------------------------------------------------------
void ViscDamp2D::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
}

//----------------------------------------------------------------------
double* ViscDamp2D::getMidPoint(CoordSet &cs)
{ 
  double* midPoint = new double[3];
  
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  
  midPoint[0] = 0.5*(nd1.x + nd2.x);
  midPoint[1] = 0.5*(nd1.y + nd2.y);
  midPoint[2] = 0.5*(nd1.z + nd2.z);
  
  return midPoint;
}

//----------------------------------------------------------------------
double ViscDamp2D::area(CoordSet& cs)
{
  Node &nd0 = cs.getNode(nn[0]);
  Node &nd1 = cs.getNode(nn[1]);
  
  const double dx    = nd1.x - nd0.x;
  const double dy    = nd1.y - nd0.y;
  const double len2  = dx*dx + dy*dy;
  return sqrt(len2)*prop->eh;
}

double ViscDamp2D::darea(CoordSet& cs, CoordSet& dcs)
{
  if(static_cast<GeoSource_opt*>(geoSource)->areThereMovingNodes())
    {
      Node &nd0    = cs.getNode( nn[0] );
      Node &nd1    = cs.getNode( nn[1] );
      
      Node &dnd0   = dcs.getNode( nn[0] );
      Node &dnd1   = dcs.getNode( nn[1] );
      
      const double dx    = nd1.x - nd0.x;
      const double dy    = nd1.y - nd0.y;
      const double len2  = dx*dx + dy*dy;
      
      const double ddx   = dnd1.x - dnd0.x;
      const double ddy   = dnd1.y - dnd0.y;
      const double dlen2 = 2.0*(dx*ddx + dy*ddy);
      
      const double len   = sqrt(len2);
      return 0.5*dlen2/len*prop->eh + sqrt(len2)*gradprop->eh;
    }
  else
    { return 0; }
}
  
//----------------------------------------------------------------------
FullSquareMatrix ViscDamp2D::dampingMatrix(CoordSet& cs, double* cel,int cmflg)
{
  FullSquareMatrix result(1, cel);
  result.setSize(numDofs());
  result.zero();

  // damping coefficients in the coordinate system, where
  // "x" axis is aligned with the bar's axis
  const double _c00 = 0.5*c00(cs);
  const double _c11 = 0.5*c11(cs);

  Node &nd0 = cs.getNode(nn[0]);
  Node &nd1 = cs.getNode(nn[1]);
  
  const double dx    = nd1.x - nd0.x;
  const double dy    = nd1.y - nd0.y;
  const double len2  = dx*dx + dy*dy;
  
  // damping matrix in the original coordinate system
  result[0][0] = (_c00*dx*dx+_c11*dy*dy)/len2;
  result[0][1] = (_c00-_c11)*dx*dy/len2;
  result[1][0] = result[0][1];
  result[1][1] = (_c11*dx*dx+_c00*dy*dy)/len2;

  result[0+dim()][0+dim()] =  result[0][0];
  result[0+dim()][1+dim()] =  result[0][1];
  result[1+dim()][0+dim()] =  result[1][0];
  result[1+dim()][1+dim()] =  result[1][1];

  return result;
}

//----------------------------------------------------------------------
Element* ViscDamp2D::clone()
{
  return new ViscDamp2D(*this);
}


//----------------------------------------------------------------------
void ViscDamp2D::gradstiffness(CoordSet&, CoordSet&, FullSquareMatrix& result, int)
{
  result.setSize(numDofs());
  result.zero();
}

//----------------------------------------------------------------------
void ViscDamp2D::gradMassMatrix(CoordSet&, CoordSet&, FullSquareMatrix& result, double mratio)
{
  result.setSize(numDofs());
  result.zero();
}

//----------------------------------------------------------------------
void ViscDamp2D::gradDampMatrix(CoordSet& cs, CoordSet& dcs, FullSquareMatrix& result, double freq)
{
  result.setSize(numDofs());
  result.zero();

  Node &nd0    = cs.getNode( nn[0] );
  Node &nd1    = cs.getNode( nn[1] );
   
  Node &dnd0   = dcs.getNode( nn[0] );
  Node &dnd1   = dcs.getNode( nn[1] );
  
  const double dx    = nd1.x - nd0.x;
  const double dy    = nd1.y - nd0.y;
  const double len2  = dx*dx + dy*dy;

  const double ddx   = dnd1.x - dnd0.x;
  const double ddy   = dnd1.y - dnd0.y;
  const double dlen2 = 2.0*(dx*ddx + dy*ddy);

  const double _c00   = 0.5*c00(cs);
  const double _c11   = 0.5*c11(cs);

  const double _dc00  = 0.5*dc00(cs, dcs);
  const double _dc11  = 0.5*dc11(cs, dcs);
  
  result[0][0] = 
    (_dc00*dx*dx+_dc11*dy*dy)/len2 
    + 2.0*(_c00*ddx*dx+_c11*ddy*dy)/len2
    - dlen2*(_c00*dx*dx+_c11*dy*dy)/(len2*len2);
  result[0][1] = 
    (_dc00-_dc11)*dx*dy/len2
    + (_c00-_c11)*ddx*dy/len2
    + (_c00-_c11)*dx*ddy/len2
    - dlen2*(_c00-_c11)*dx*dy/(len2*len2);
  result[1][0] = result[0][1];
  result[1][1] =
    (_dc11*dx*dx+_dc00*dy*dy)/len2 
    + 2.0*(_c11*ddx*dx+_c00*ddy*dy)/len2
    - dlen2*(_c11*dx*dx+_c00*dy*dy)/(len2*len2);

  result[0+dim()][0+dim()] =  result[0][0];
  result[0+dim()][1+dim()] =  result[0][1];
  result[1+dim()][0+dim()] =  result[1][0];
  result[1+dim()][1+dim()] =  result[1][1];
}

//----------------------------------------------------------------------
int ViscDamp2D::chkOptInf(CoordSet& d_cs)
{
  if(prop == 0 || gradprop == 0) { return 0; }

  Node &gnd0 = d_cs.getNode(nn[0]);
  Node &gnd1 = d_cs.getNode(nn[1]);

  double drho   = gradprop->rho;
  double dE     = gradprop->E;
  double dnu    = gradprop->nu;
  double dh     = gradprop->eh;

  double d_vd = 
    gnd0.x*gnd0.x + gnd0.y*gnd0.y + 
    gnd1.x*gnd1.x + gnd1.y*gnd1.y;

  double d_at = drho*drho + dnu*dnu + dE*dE + dh*dh;

  if ( d_vd > 0 && d_at > 0 ) { return Optvar::mixedvar; }
  if ( d_vd > 0 )	      { return Optvar::coordinate; }
  if ( d_at > 0 )	      { return Optvar::attribute; }
  return 0;  
}

extern PrioInfo examineBar2(int sub, MultiFront *mf, int *nn);
//----------------------------------------------------------------------
PrioInfo ViscDamp2D::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

