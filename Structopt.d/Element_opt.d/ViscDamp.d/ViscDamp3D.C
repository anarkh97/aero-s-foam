#ifdef STRUCTOPT

#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <cassert>

#include <Utils.d/dofset.h>
#include <Structopt.d/Element_opt.d/ViscDamp.d/ViscDamp3D.h>
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Driver_opt.d/GeoSource_opt.h>

//----------------------------------------------------------------------
// class constants
template<int nNodes> 
const int ViscDamp3D<nNodes>::nDim = 3;

//----------------------------------------------------------------------
template<int nNodes>
ViscDamp3D<nNodes>::ViscDamp3D(int* nodenums)
{
  assert(nNodes >= 3);
  std::copy(nodenums, nodenums+nNodes, nn);
}

//----------------------------------------------------------------------
template<int nNodes>
int* ViscDamp3D<nNodes>::nodes(int* p)
{
  if(p == 0) 
    {
      p = new int[nNodes];
    }
  std::copy(nn, nn+nNodes, p);
  return p;
}

//----------------------------------------------------------------------
template<int nNodes>
int* ViscDamp3D<nNodes>::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) 
    {
      p = new int[numDofs()];
    }
  for(int i = 0; i < nNodes; ++i)
    {
      dsa.number(nn[i], DofSet::XYZdisp, p + i*dim());
    }
  return p;
}

//----------------------------------------------------------------------
template<int nNodes>
void ViscDamp3D<nNodes>::markDofs( DofSetArray &dsa )
{
  dsa.mark(nn, nNodes, DofSet::XYZdisp);
}

//----------------------------------------------------------------------
template<int nNodes>
void ViscDamp3D<nNodes>::renum(int* table)
{
  for(int i=0; i < nNodes; ++i)
    {
      nn[i] = table[nn[i]];
    }
}

//----------------------------------------------------------------------
template<int nNodes>
double* ViscDamp3D<nNodes>::getMidPoint(CoordSet &cs)
{ 
  double* midPoint = new double[dim()];
  std::fill(midPoint, midPoint+dim(), 0.0);
  
  for(int i=0; i<nNodes; ++i)
    {
      Node& nd = cs.getNode(nn[i]);  
      midPoint[0] += nd.x;
      midPoint[1] += nd.y;
      midPoint[2] += nd.z;
    }
  midPoint[0] /= nNodes;
  midPoint[1] /= nNodes;
  midPoint[2] /= nNodes;
  
  return midPoint;
}

//----------------------------------------------------------------------
template<int nNodes>
FullSquareMatrix ViscDamp3D<nNodes>::dampingMatrix(CoordSet& cs, double* cel, int)
{
  FullSquareMatrix result(1, cel);
  result.setSize(numDofs());
  result.zero();

  Vector v0,  v1,  v2;
  Vector v0t, v1t, v2t;
  Vector w0t, w1t, w2t;
  orth_basis(cs, v0, v1, v2);
  transpose(v0, v1, v2, v0t, v1t, v2t);
  v0 *= c00(cs)/nNodes; v1 *= c00(cs)/nNodes; v2 *= c11(cs)/nNodes;
  transpose(v0, v1, v2, w0t, w1t, w2t);
  // damping matrix in the original coordinate system
  result[0][0] = w0t*v0t;
  result[0][1] = w0t*v1t;
  result[0][2] = w0t*v2t;
  result[1][1] = w1t*v1t;
  result[1][2] = w1t*v2t;
  result[2][2] = w2t*v2t;
  result[1][0] = result[0][1];
  result[2][0] = result[0][2];
  result[2][1] = result[1][2];
  for(int i = dim(); i < numDofs(); i += dim())
    {
      for(int j = 0; j < dim(); ++j)
	{
	  for(int k = 0; k < dim(); ++k)
	    {
	      result[i+j][i+k] = result[j][k];
	    }
	}
    }
  return result;
}

//----------------------------------------------------------------------
template<int nNodes>
Vector ViscDamp3D<nNodes>::avg_normal(CoordSet& cs)
{
  int s = 0;
  Vector sum(dim(), 0.0);
  Vector v0(dim(), 0.0);
  Vector v1(dim(), 0.0);
  for(int j=0; j<nNodes-2; ++j)
    {
      Node& nd0 = cs.getNode(nn[j]);
      Node& nd1 = cs.getNode(nn[j+1]);
      Node& nd2 = cs.getNode(nn[j+2]);
      v0[0] = nd1.x - nd0.x;
      v0[1] = nd1.y - nd0.y;
      v0[2] = nd1.z - nd0.z;      
      v1[0] = nd2.x - nd1.x;
      v1[1] = nd2.y - nd1.y;
      v1[2] = nd2.z - nd1.z;
      sum += v0.cross(v1);
      ++s;
    }
  sum *= 1.0/s;
  return sum;
}

//----------------------------------------------------------------------
template<int nNodes>
void ViscDamp3D<nNodes>::orth_basis(CoordSet& cs, Vector& v0, Vector& v1, Vector& v2)
{
  v0 = Vector(dim(), 0.0);
  v1 = Vector(dim(), 0.0);
  Node &nd0 = cs.getNode(nn[0]);
  Node &nd1 = cs.getNode(nn[1]);
  Node &nd2 = cs.getNode(nn[2]);
  v0[0] = nd1.x - nd0.x;
  v0[1] = nd1.y - nd0.y;
  v0[2] = nd1.z - nd0.z;      
  v1[0] = nd2.x - nd1.x;
  v1[1] = nd2.y - nd1.y;
  v1[2] = nd2.z - nd1.z;


  v2 = avg_normal(cs);
  double l0 = v0.magnitude();
  if( l0 == 0.0 )
    {
      fprintf(stderr, "*** ERROR: element has zero edge");
      for(int i = 0; i < nNodes; ++i)
	{
	  fprintf(stderr, " %d", nn[i]+1);
	}
      fprintf(stderr, "\n");
    }
  v0 *= 1.0/l0;

  v1.linAdd(-(v1 * v0), v0);
  double l1 = v1.magnitude();
  if( l1 == 0.0 )
    {
      fprintf(stderr, "*** ERROR: consequent edges are linearly dependent in the element");
      for(int i = 0; i < nNodes; ++i)
	{
	  fprintf(stderr, " %d", nn[i]+1);
	}
      fprintf(stderr, "\n");
    }
  v1 *= 1.0/l1;

  v2.linAdd(-(v2 * v0), v0, -(v2 * v1), v1);
  double l2 = v2.magnitude();
  if( l2 == 0.0 )
    {
      fprintf(stderr, "*** ERROR: cannot compute normal for the element");
      for(int i = 0; i < nNodes; ++i)
	{
	  fprintf(stderr, " %d", nn[i]+1);
	}
      fprintf(stderr, "\n");
    }
  v2 *= 1.0/l2;
}

//----------------------------------------------------------------------
template<int nNodes>
void ViscDamp3D<nNodes>::transpose(Vector& v0,  Vector& v1,  Vector& v2, 
				   Vector& v0t, Vector& v1t, Vector& v2t)
{
  v0t = Vector(dim(), 0.0);
  v1t = Vector(dim(), 0.0);
  v2t = Vector(dim(), 0.0);

  v0t[0] = v0[0]; v0t[1] = v1[0]; v0t[2] = v2[0];
  v1t[0] = v0[1]; v1t[1] = v1[1]; v1t[2] = v2[1];
  v2t[0] = v0[2]; v2t[1] = v1[2]; v2t[2] = v2[2];
}

//----------------------------------------------------------------------
template<int nNodes>
double ViscDamp3D<nNodes>::area(CoordSet& cs)
{
  Node& nd0 = cs.getNode(nn[0]);
  Node& nd1 = cs.getNode(nn[1]);
  Vector v0(dim(), 0.0);
  Vector v1(dim(), 0.0);

  v1[0] = nd1.x - nd0.x;
  v1[1] = nd1.y - nd0.y;
  v1[2] = nd1.z - nd0.z;
  
  double l0 = 0;
  double l1 = v1.magnitude();

  double result = 0.0;

  for(int i = 2; i<nNodes; ++i)
    {
      v0 = v1;
      l0 = l1;
      Node& nd2 = cs.getNode(nn[i]);
      v1[0] = nd2.x - nd0.x;
      v1[1] = nd2.y - nd0.y;
      v1[2] = nd2.z - nd0.z;
      l1 = v1.magnitude();
      double v0v1  = v0*v1;
      double l0l1  = l0*l1;
      double area3 = 0.5*sqrt(l0l1*l0l1-v0v1*v0v1);
      result += area3;
    }
  return result;
}


//----------------------------------------------------------------------
template<int nNodes>
double ViscDamp3D<nNodes>::darea(CoordSet& cs, CoordSet& dcs)
{
  if(static_cast<GeoSource_opt*>(geoSource)->areThereMovingNodes())
    {
      int i = 1;
      Node& nd0 = cs.getNode(nn[0]);
      Node& nd1 = cs.getNode(nn[1]);
      Vector v0(dim(), 0.0);
      Vector v1(dim(), 0.0);
      
      v1[0] = nd1.x - nd0.x;
      v1[1] = nd1.y - nd0.y;
      v1[2] = nd1.z - nd0.z;
      
      Node& dnd0 = dcs.getNode(nn[0]);
      Node& dnd1 = dcs.getNode(nn[1]);
      Vector dv0(dim(), 0.0);
      Vector dv1(dim(), 0.0);
      
      dv1[0] = dnd1.x - dnd0.x;
      dv1[1] = dnd1.y - dnd0.y;
      dv1[2] = dnd1.z - dnd0.z;
      
      double l0_2 = 0;
      double l1_2 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
      double dl0_2 = 0;
      double dl1_2 = 2.0*(v1[0]*dv1[0] + v1[1]*dv1[1] + v1[2]*dv1[2]);
      
      double result = 0.0;
      for(int i = 2; i<nNodes; ++i)
	{
	  v0 = v1; dv0 = dv1;
	  l0_2 = l1_2; dl0_2 = dl1_2;
	  Node& nd2 = cs.getNode(nn[i]);
	  Node& dnd2 = dcs.getNode(nn[i]);
	  v1[0] = nd2.x - nd0.x;
	  v1[1] = nd2.y - nd0.y;
	  v1[2] = nd2.z - nd0.z;
	  
	  dv1[0] = dnd2.x - dnd0.x;
	  dv1[1] = dnd2.y - dnd0.y;
	  dv1[2] = dnd2.z - dnd0.z;
	  
	  l1_2 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
	  dl1_2 = 2.0*(v1[0]*dv1[0] + v1[1]*dv1[1] + v1[2]*dv1[2]);
	  
	  double v0v1   = v0*v1;      
	  double d_v0v1 = dv0*v1 + v0*dv1;
	  double area3  = 0.5*sqrt(l0_2*l1_2-v0v1*v0v1);
	  double darea3 = 0.25*(dl0_2*l1_2+l0_2*dl1_2 - 2.0*v0v1*d_v0v1)/area3;
	  result += darea3;
	}
      return result;
    }
  else
    {
      return 0;
    }
}

//----------------------------------------------------------------------
template<int nNodes>
Element* ViscDamp3D<nNodes>::clone()
{
  return new ViscDamp3D<nNodes>(*this);
}


//----------------------------------------------------------------------
template<int nNodes>
void ViscDamp3D<nNodes>::gradstiffness(CoordSet&, CoordSet&, FullSquareMatrix& result, int)
{
  result.setSize(numDofs());
  result.zero();
}

//----------------------------------------------------------------------
template<int nNodes>
void ViscDamp3D<nNodes>::gradMassMatrix(CoordSet&, CoordSet&, FullSquareMatrix& result, double mratio)
{
  result.setSize(numDofs());
  result.zero();
}

//----------------------------------------------------------------------
template<int nNodes>
void ViscDamp3D<nNodes>::gradDampMatrix(CoordSet& cs, CoordSet& dcs, FullSquareMatrix& result, double freq)
{
  result.setSize(numDofs());
  result.zero();

  Vector v0,  v1,  v2;
  Vector v0t, v1t, v2t;
  Vector w0t, w1t, w2t;
  orth_basis(cs, v0, v1, v2);
  transpose(v0, v1, v2, v0t, v1t, v2t);
  v0 *= dc00(cs, dcs)/nNodes; v1 *= dc00(cs, dcs)/nNodes; v2 *= dc11(cs, dcs)/nNodes;
  transpose(v0, v1, v2, w0t, w1t, w2t);
  
  // damping matrix in the original coordinate system
  result[0][0] = w0t*v0t;
  result[0][1] = w0t*v1t;
  result[0][2] = w0t*v2t;
  result[1][1] = w1t*v1t;
  result[1][2] = w1t*v2t;
  result[2][2] = w2t*v2t;
  result[1][0] = result[0][1];
  result[2][0] = result[0][2];
  result[2][1] = result[1][2];

  for(int i = dim(); i < numDofs(); i += dim())
    {
      for(int j = 0; j < dim(); ++j)
	{
	  for(int k = 0; k < dim(); ++k)
	    {
	      result[i+j][i+k] = result[j][k];
	    }
	}
    }
}


//----------------------------------------------------------------------
template<int nNodes>
Vector ViscDamp3D<nNodes>::davg_normal(CoordSet& cs, CoordSet& dcs)
{
  int s = 0;
  Vector sum(dim(), 0.0);
  Vector  v0(dim(), 0.0);
  Vector  v1(dim(), 0.0);
  Vector dv0(dim(), 0.0);
  Vector dv1(dim(), 0.0);

  for(int i=0; i<nNodes-2; ++i)
    {
      Node&  nd0 =  cs.getNode(nn[i]);
      Node& dnd0 = dcs.getNode(nn[i]);
      for(int j=i+1; j<nNodes-1; ++j)
	{
	  Node&  nd1 =  cs.getNode(nn[j]);
	  Node& dnd1 = dcs.getNode(nn[j]);

	  v0[0]  = nd1.x  - nd0.x;
	  v0[1]  = nd1.y  - nd0.y;
	  v0[2]  = nd1.z  - nd0.z;      

	  dv0[0] = dnd1.x - dnd0.x;
	  dv0[1] = dnd1.y - dnd0.y;
	  dv0[2] = dnd1.z - dnd0.z;      
	  for(int k=j+1; k<nNodes; ++k)
	    {
	      Node&  nd2 =  cs.getNode(nn[k]);
	      Node& dnd2 = dcs.getNode(nn[k]);
	      v1[0]  = nd2.x  - nd1.x;
	      v1[1]  = nd2.y  - nd1.y;
	      v1[2]  = nd2.z  - nd1.z;

	      dv1[0] = dnd2.x - dnd1.x;
	      dv1[1] = dnd2.y - dnd1.y;
	      dv1[2] = dnd2.z - dnd1.z;

	      sum += dv0.cross(v1);
	      sum += v0.cross(dv1);
	      ++s;
	    }
	}
    }
  sum *= 1.0/s;
  return sum;
}


//----------------------------------------------------------------------
template<int nNodes>
int ViscDamp3D<nNodes>::chkOptInf(CoordSet& d_cs)
{
  if(prop == 0 || gradprop == 0) { return 0; }

  double drho   = gradprop->rho;
  double dE     = gradprop->E;
  double dnu    = gradprop->nu;

  double d_vd = 0.0;
  for(int i = 0; i < nNodes; ++i)
    {
      Node& gnd = d_cs.getNode(nn[i]);
      d_vd += gnd.x*gnd.x + gnd.y*gnd.y + gnd.z*gnd.z;
    }

  double d_at = drho*drho + dnu*dnu + dE*dE;

  if( d_vd > 0 )
    {
      fprintf( stderr, "ERROR: derivatives w.r.t. changes in nodal positions are not implemented!!!\n" );
    }

  if ( d_at > 0 ) { return Optvar::attribute; }
  return 0;  
}

extern PrioInfo examineQuad4(int sub, MultiFront *mf, int *nn);
//----------------------------------------------------------------------
template<> inline
PrioInfo ViscDamp3D<4>::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

#endif
