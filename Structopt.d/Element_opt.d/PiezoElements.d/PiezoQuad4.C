#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/PiezoElements.d/PiezoQuad4.h>
#include <Structopt.d/Optvar.h>
#include <cassert>
#include <algorithm>
#include <cstdlib>

//------------------------------------------------------------------------------
int* PiezoQuad4::dofs(DofSetArray& dsa, int *p)
{
  if(p == 0) 
    {
      p = new int[numDofs()];
      assert(p!=0);
    }
  for(int i = 0; i < numNodes(); ++i)
    { dsa.number(nn[i], DofSet::Xdisp|DofSet::Ydisp|DofSet::Temp, p + i*dim()); }
  return p;
}

//------------------------------------------------------------------------------
FullSquareMatrix PiezoQuad4::stiffness(CoordSet& cs, double *d, int)
{
  FullSquareMatrix result(numDofs(), d);
  result.zero();
  
  // [kappa] is a 2nd rank dielectric tensor
  double kappa[2][2];
  // [e] is a 3rd rank tensor of piezoelectricity
  double e[3][2];
  fillPiezoTensors(e, kappa);

  double X[4], Y[4], h[4];
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);
  
  X[0] = nd1.x; Y[0] = nd1.y;
  X[1] = nd2.x; Y[1] = nd2.y;
  X[2] = nd3.x; Y[2] = nd3.y;
  X[3] = nd4.x; Y[3] = nd4.y;

  h[0] = h[1] = h[2] = h[3] = prop->eh;

  // order of gaussian quadrature
  int p = 2;

  // now we are ready to evaluate the stiffness mx
  piezo_quad4v(X, Y, h, p, kappa, e, d);

  return result;
}

extern "C" {
  void _FORTRAN(qgauss)(int*, int*, int*, int*, double*, double*, double*);
  void _FORTRAN(q4shpe)(double*, double*, const double*, const double*, 
			double*, double*, double*, double*);
}

//------------------------------------------------------------------------------
void PiezoQuad4::piezo_quad4v(const double* x, const double* y, const double* h,
			      int p, 
			      const double kappa[2][2],
			      const double e[3][2],
			      double* d)
{
  // for numbering of [kappa], [e] see paper by C.M. Landis, IJNME 2002, 55:613-628, p. 615

  const int cNumNodes = 4;
  const int cDim      = 3;
  const int cNumDofs  = numDofs();
  const int dofX=0, dofY=1, dofT=2;
  const int ls[cNumNodes][cDim] = {{0,1,2},{3,4,5},
				   {6,7,8},{9,10,11}};
  double xi, eta, det, w, weight;
  double q[cNumNodes], qx[cNumNodes], qy[cNumNodes];
  
  std::fill(d, d+numDofs()*numDofs(), 0);
  for(int k=1; k<=p; ++k)
    {
      for(int l=1; l<=p; ++l)
	{
	  _FORTRAN(qgauss)(&p,&k,&p,&l,&xi,&eta,&weight);
	  _FORTRAN(q4shpe)(&xi,&eta,x,y,q,qx,qy,&det);
	  if(det<=0)
	    { 
	      fprintf(stderr, "Negative Jacobian determinant in PiezoQuad4\n");
	      exit(-1);
	    }
	  w = weight*det*(h[0]*q[0]+h[1]*q[1]+h[2]*q[2]+h[3]*q[3]);
	  // calculate the piezo-piezo part
	  for(int j=0; j<cNumNodes; ++j)
	    {
	      int jT = ls[j][dofT];
	      double kappa_x = (kappa[0][0]*qx[j]+kappa[0][1]*qy[j])*w;
	      double kappa_y = (kappa[1][0]*qx[j]+kappa[1][1]*qy[j])*w;
	      for(int i=j; i<cNumNodes; ++i)
		{
		  int iT = ls[i][dofT];
		  d[iT*cNumDofs+jT] -= kappa_x*qx[i] + kappa_y*qy[i];
		  d[jT*cNumDofs+iT] = d[iT*cNumDofs+jT];
		}
	    }
	  
	  // calculate the piezo-structure part
	  for(int j=0; j<cNumNodes; ++j)
	    {
	      int jx = ls[j][dofX];
	      int jy = ls[j][dofY];
	      int jT = ls[j][dofT];
	      
	      double e_0x = (e[0][0]*qx[j]+e[2][0]*qy[j])*w;
	      double e_0y = (e[2][0]*qx[j]+e[1][0]*qy[j])*w;
	      
	      double e_1x = (e[0][1]*qx[j]+e[2][1]*qy[j])*w;
	      double e_1y = (e[2][1]*qx[j]+e[1][1]*qy[j])*w;
	      
	      for(int i=0; i<cNumNodes; ++i)
		{
		  int ix = ls[i][dofX];
		  int iy = ls[i][dofY];
		  int iT = ls[i][dofT];
		  
		  d[jx*cNumDofs+iT] += e_0x*qx[i]+e_1x*qy[i];
		  d[jy*cNumDofs+iT] += e_0y*qx[i]+e_1y*qy[i];
		  
		  d[iT*cNumDofs+jx] = d[jx*cNumDofs+iT];
		  d[iT*cNumDofs+jy] = d[jy*cNumDofs+iT];
		}
	    }
	}
    }
  
  return;
}

//------------------------------------------------------------------------------
void PiezoQuad4::piezo_dquad4v(const double* x,  const double*  y, const double*  h,
			       const double* dx, const double* dy, const double* dh,
			       int p, 
			       const double kappa[2][2], const double dkappa[2][2],
			       const double e[3][2],     const double de[3][2],
			       double* d)
{ 
  const int cNumNodes = 4;
  const int cDim      = 3;
  const int cNumDofs  = numDofs();
  const int dofX=0, dofY=1, dofT=2;
  const int ls[cNumNodes][cDim] = {{0,1,2},{3,4,5},
				   {6,7,8},{9,10,11}};
  double xi, eta, det, w, dw, weight;
  bool cvar = false;
  double q[cNumNodes], qx[cNumNodes], qy[cNumNodes];

  for(int i=0; i<cNumNodes && !cvar; ++i) 
    { cvar |= (dx[i] != 0) | (dy[i] != 0); }
  
  std::fill(d, d+numDofs()*numDofs(), 0);

  if(!cvar)
    {
      for(int k=1; k<=p; ++k)
	{
	  for(int l=1; l<=p; ++l)
	    {
	      _FORTRAN(qgauss)(&p,&k,&p,&l,&xi,&eta,&weight);
	      _FORTRAN(q4shpe)(&xi,&eta,x,y,q,qx,qy,&det);
	      if(det<=0)
		{ 
		  fprintf(stderr, "Negative Jacobian determinant in PiezoQuad4\n");
		  exit(-1);
		}
	      w  = weight*det*(h[0]*q[0]+h[1]*q[1]+h[2]*q[2]+h[3]*q[3]);
	      dw = weight*det*(dh[0]*q[0]+dh[1]*q[1]*dh[2]*q[2]*dh[3]*q[3]);
	      // calculate the piezo-piezo part
	      for(int j=0; j<cNumNodes; ++j)
		{
		  int jT = ls[j][dofT];
		  double dkappa_x = 
		    (dkappa[0][0]*qx[j]+dkappa[0][1]*qy[j])*w +
		    (kappa[0][0]*qx[j]+kappa[0][1]*qy[j])*dw;
		  double dkappa_y = 
		    (dkappa[1][0]*qx[j]+dkappa[1][1]*qy[j])*w +
		    (kappa[1][0]*qx[j]+kappa[1][1]*qy[j])*dw;
		  for(int i=j; i<cNumNodes; ++i)
		    {
		      int iT = ls[i][dofT];
		      d[iT*cNumDofs+jT] -= dkappa_x*qx[i] + dkappa_y*qy[i];
		      d[jT*cNumDofs+iT] = d[iT*cNumDofs+jT];
		    }
		}
	      // calculate the piezo-structure part
		  for(int j=0; j<cNumNodes; ++j)
		    {
		      int jx = ls[j][dofX];
		      int jy = ls[j][dofY];
		      int jT = ls[j][dofT];
		  
		      double de_0x = 
			(de[0][0]*qx[j]+de[2][0]*qy[j])*w +
			(e[0][0]*qx[j]+e[2][0]*qy[j])*dw;			
		      double de_0y = 
			(de[2][0]*qx[j]+de[1][0]*qy[j])*w +
			(e[2][0]*qx[j]+e[1][0]*qy[j])*dw;
		      
		      double de_1x = 
			(de[0][1]*qx[j]+de[2][1]*qy[j])*w +
			(e[0][1]*qx[j]+e[2][1]*qy[j])*dw;
		      double de_1y = 
			(de[2][1]*qx[j]+de[1][1]*qy[j])*w +
			(e[2][1]*qx[j]+e[1][1]*qy[j])*dw;
		      
		      for(int i=0; i<cNumNodes; ++i)
			{
			  int ix = ls[i][dofX];
			  int iy = ls[i][dofY];
			  int iT = ls[i][dofT];
			  
			  d[jx*cNumDofs+iT] += de_0x*qx[i]+de_1x*qy[i];
			  d[jy*cNumDofs+iT] += de_0y*qx[i]+de_1y*qy[i];
			  
			  d[iT*cNumDofs+jx] = d[jx*cNumDofs+iT];
			  d[iT*cNumDofs+jy] = d[jy*cNumDofs+iT];
			}
		    }		  
	    }
	}
    }
  else
    {
      // not implemented yet
      assert(0);
    }
  return;
}

//------------------------------------------------------------------------------
int PiezoQuad4::chkOptInf(CoordSet& dcs)
{
  if(cCoefs==0 || getDCCoefs()==0 || prop==0 || gradprop==0) { return 0; }

  double dX[4], dY[4];
  // gradient orientated data
  Node &d_nd1 = dcs.getNode(nn[0]);
  Node &d_nd2 = dcs.getNode(nn[1]);
  Node &d_nd3 = dcs.getNode(nn[2]);
  Node &d_nd4 = dcs.getNode(nn[3]);
  
  dX[0] = d_nd1.x; dY[0] = d_nd1.y;
  dX[1] = d_nd2.x; dY[1] = d_nd2.y;
  dX[2] = d_nd3.x; dY[2] = d_nd3.y;
  dX[3] = d_nd4.x; dY[3] = d_nd4.y;      

  double dvd = 0;
  for(int i=0; i<numNodes(); ++i) { dvd += dX[i]*dX[i] + dY[i]*dY[i]; }
  double dat = 0;
  for(int j=0; j<6; ++j) { for(int i=0; i<=j; ++i) { dat += getDCCoefs()[i   + j*6]*getDCCoefs()[i   + j*6]; } }
  dat += (gradprop->eh!=0);
  
  if ( dvd > 0 && dat >0 ) return Optvar::mixedvar;
  if ( dvd > 0 )           return Optvar::coordinate;
  if ( dat > 0 )           return Optvar::attribute;

  return 0;
}

//------------------------------------------------------------------------------
void PiezoQuad4::gradstiffness(CoordSet& cs,CoordSet& dcs,FullSquareMatrix& dd, int flg)
{
  // [kappa] is a 2nd rank dielectric tensor
  double kappa[2][2], dkappa[2][2];
  // [e] is a 3rd rank tensor of piezoelectricity
  double e[3][2], de[3][2];
  fillPiezoTensors(e, kappa);
  fillDPiezoTensors(de, dkappa);

  double  X[4],  Y[4],  h[4];
  double dX[4], dY[4], dh[4];

  std::fill(h, h+numNodes(), prop->eh);
  std::fill(dh, dh+numNodes(), gradprop->eh);

  // gradient orientated data
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);
  
  X[0] = nd1.x; Y[0] = nd1.y;
  X[1] = nd2.x; Y[1] = nd2.y;
  X[2] = nd3.x; Y[2] = nd3.y;
  X[3] = nd4.x; Y[3] = nd4.y;

  // gradient orientated data
  Node &d_nd1 = dcs.getNode(nn[0]);
  Node &d_nd2 = dcs.getNode(nn[1]);
  Node &d_nd3 = dcs.getNode(nn[2]);
  Node &d_nd4 = dcs.getNode(nn[3]);
  
  dX[0] = d_nd1.x; dY[0] = d_nd1.y;
  dX[1] = d_nd2.x; dY[1] = d_nd2.y;
  dX[2] = d_nd3.x; dY[2] = d_nd3.y;
  dX[3] = d_nd4.x; dY[3] = d_nd4.y;      

  // order of gaussian quadrature
  int p = 2;

  // now we are ready to evaluate the stiffness mx
  piezo_dquad4v(X, Y, h, dX, dY, dh, p, kappa, dkappa, e, de, dd.data());
  return;
}

//------------------------------------------------------------------------------
/*
  in the input file: use COEF keyword
  Since the matrix is assumed to be symmetric there, we have "only" 21 coefficients.
  Therefore, please input coefficients as follows:
  1 1 1. -> e11
  1 2 2. -> e12
  1 3 3. -> e13
  1 4 4. -> e21
  1 5 5. -> e22
  1 6 6. -> e23
  2 2 7. -> e31
  2 3 8. -> e32
  2 4 9. -> e33
  2 5 10. -> e41
  2 6 11. -> e42
  3 3 12. -> e43
  3 4 13. -> e51
  3 5 14. -> e52
  3 6 15. -> e53
  4 4 16. -> e61
  4 5 17. -> e62
  4 6 18. -> e63
  5 5 19. -> k11
  5 6 20. -> k22
  6 6 21. -> k33
  First 18 coefficients go into the piezo-electric matrix, 6x3.
  Last  3  coefficients end up on the diagonal of the dielectric matrix.
 */
void PiezoQuad4::fillPiezoTensors(double e[3][2], double kappa[2][2])
{  
  assert(cCoefs != 0);

  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  kappa[i][j] = 0.0;
	}
    }

  e[0][0] = cCoefs[0   + 0*6];
  e[0][1] = cCoefs[0   + 1*6];

  e[1][0] = cCoefs[0   + 2*6];
  e[1][1] = cCoefs[0   + 3*6];

  e[2][0] = cCoefs[0   + 4*6];
  e[2][1] = cCoefs[0   + 5*6];

  kappa[0][0] = cCoefs[1   + 1*6];
  kappa[1][1] = cCoefs[1   + 2*6];

  if(cFrame != 0)
    {
      //FIXIT: tensors are not rotated!!!!!
    }
  return;
}

void PiezoQuad4::fillDPiezoTensors(double de[3][2], double dkappa[2][2])
{
  assert(getDCCoefs() != 0);

  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  dkappa[i][j] = 0.0;
	}
    }

  de[0][0] = getDCCoefs()[0   + 0*6];
  de[0][1] = getDCCoefs()[0   + 1*6];

  de[1][0] = getDCCoefs()[0   + 2*6];
  de[1][1] = getDCCoefs()[0   + 3*6];

  de[2][0] = getDCCoefs()[0   + 4*6];
  de[2][1] = getDCCoefs()[0   + 5*6];

  dkappa[0][0] = getDCCoefs()[1   + 1*6];
  dkappa[1][1] = getDCCoefs()[1   + 2*6];

  if(cFrame != 0)
    {
      //FIXIT: tensors are not rotated!!!!!
    }
  return;
}

#endif
