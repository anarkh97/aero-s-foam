#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/PiezoElements.d/PiezoBrick8.h>
#include <Structopt.d/Optvar.h>
#include <cassert>
#include <algorithm>
#include <cstdlib>

//------------------------------------------------------------------------------
int* PiezoBrick8::dofs(DofSetArray& dsa, int *p)
{
  if(p == 0) 
    {
      p = new int[numDofs()];
      assert(p!=0);
    }
  for(int i = 0; i < numNodes(); ++i)
    { dsa.number(nn[i], DofSet::XYZdisp|DofSet::Temp, p + i*dim()); }
  return p;
}

//------------------------------------------------------------------------------
FullSquareMatrix PiezoBrick8::stiffness(CoordSet& cs, double *d, int)
{
  FullSquareMatrix result(numDofs(), d);
  result.zero();
  
  // [kappa] is a 2nd rank dielectric tensor
  double kappa[3][3];
  // [e] is a 3rd rank tensor of piezoelectricity
  double e[6][3];
  fillPiezoTensors(e, kappa);

  double X[8], Y[8], Z[8];
  cs.getCoordinates(nn, numNodes(), X, Y, Z);

  // order of gaussian quadrature
  int p = 2;

  // now we are ready to evaluate the stiffness mx
  piezo_br8v(X, Y, Z, p, kappa, e, d);

  return result;
}

extern "C" {
  void _FORTRAN(hxgaus)(int*, int*, int*, int*, int*, int*, 
			double*, double*, double*, double*);
  void _FORTRAN(h8shpe)(double*, double*, double*, const double*, const double*, const double*, 
			double*, double*, double*, double*, double*);
}

//------------------------------------------------------------------------------
void PiezoBrick8::piezo_br8v(const double* x, const double* y, const double* z,
			     int p, 
			     const double kappa[3][3],
			     const double e[6][3],
			     double* d)
{
  // for numbering of [kappa], [e] see paper by C.M. Landis, IJNME 2002, 55:613-628, p. 615

  const int cNumNodes = 8;
  const int cDim      = 4;
  const int cNumDofs  = numDofs();
  const int dofX=0, dofY=1, dofZ=2, dofT=3;
  const int ls[cNumNodes][cDim] = {{0,1,2,3},{4,5,6,7},
				   {8,9,10,11},{12,13,14,15},
				   {16,17,18,19},{20,21,22,23},
				   {24,25,26,27},{28,29,30,31}};
  double xi, eta, emu, det, w, weight;
  double q[cNumNodes], qx[cNumNodes], qy[cNumNodes], qz[cNumNodes];
  
  std::fill(d, d+numDofs()*numDofs(), 0);
  for(int k=1; k<=p; ++k)
    {
      for(int l=1; l<=p; ++l)
	{
	  for(int m=1; m<=p; ++m)
	    {
	      _FORTRAN(hxgaus)(&p,&k,&p,&l,&p,&m,&xi,&eta,&emu,&weight);
	      _FORTRAN(h8shpe)(&xi,&eta,&emu,x,y,z,q,qx,qy,qz,&det);
	      if(det==0)
		{ 
		  fprintf(stderr, "Zero Jacobian determinant in PiezoBrick8\n");
		  exit(-1);
		}
	      w = weight*det;
	      // calculate the piezo-piezo part
	      for(int j=0; j<cNumNodes; ++j)
		{
		  int jT = ls[j][dofT];
		  double kappa_x = (kappa[0][0]*qx[j]+kappa[0][1]*qy[j]+kappa[0][2]*qz[j])*w;
		  double kappa_y = (kappa[1][0]*qx[j]+kappa[1][1]*qy[j]+kappa[1][2]*qz[j])*w;
		  double kappa_z = (kappa[2][0]*qx[j]+kappa[2][1]*qy[j]+kappa[2][2]*qz[j])*w;
		  for(int i=j; i<cNumNodes; ++i)
		    {
		      int iT = ls[i][dofT];
		      d[iT*cNumDofs+jT] -= kappa_x*qx[i] + kappa_y*qy[i] + kappa_z*qz[i];
		      d[jT*cNumDofs+iT] = d[iT*cNumDofs+jT];
		    }
		}

	      // calculate the piezo-structure part
	      for(int j=0; j<cNumNodes; ++j)
		{
		  int jx = ls[j][dofX];
		  int jy = ls[j][dofY];
		  int jz = ls[j][dofZ];
		  int jT = ls[j][dofT];
		  
		  double e_0x = (e[0][0]*qx[j]+e[3][0]*qy[j]+e[5][0]*qz[j])*w;
		  double e_0y = (e[3][0]*qx[j]+e[1][0]*qy[j]+e[4][0]*qz[j])*w;
		  double e_0z = (e[5][0]*qx[j]+e[4][0]*qy[j]+e[2][0]*qz[j])*w;

		  double e_1x = (e[0][1]*qx[j]+e[3][1]*qy[j]+e[5][1]*qz[j])*w;
		  double e_1y = (e[3][1]*qx[j]+e[1][1]*qy[j]+e[4][1]*qz[j])*w;
		  double e_1z = (e[5][1]*qx[j]+e[4][1]*qy[j]+e[2][1]*qz[j])*w;

		  double e_2x = (e[0][2]*qx[j]+e[3][2]*qy[j]+e[5][2]*qz[j])*w;
		  double e_2y = (e[3][2]*qx[j]+e[1][2]*qy[j]+e[4][2]*qz[j])*w;
		  double e_2z = (e[5][2]*qx[j]+e[4][2]*qy[j]+e[2][2]*qz[j])*w;

		  for(int i=0; i<cNumNodes; ++i)
		    {
		      int ix = ls[i][dofX];
		      int iy = ls[i][dofY];
		      int iz = ls[i][dofZ];
		      int iT = ls[i][dofT];

		      d[jx*cNumDofs+iT] += e_0x*qx[i]+e_1x*qy[i]+e_2x*qz[i];
		      d[jy*cNumDofs+iT] += e_0y*qx[i]+e_1y*qy[i]+e_2y*qz[i];
		      d[jz*cNumDofs+iT] += e_0z*qx[i]+e_1z*qy[i]+e_2z*qz[i];

		      d[iT*cNumDofs+jx] = d[jx*cNumDofs+iT];
		      d[iT*cNumDofs+jy] = d[jy*cNumDofs+iT];
		      d[iT*cNumDofs+jz] = d[jz*cNumDofs+iT];
		    }
		}
	    }
	}
    }
  
  return;
}

//------------------------------------------------------------------------------
void PiezoBrick8::piezo_dbr8v(const double* x,  const double*  y, const double*  z,
			      const double* dx, const double* dy, const double* dz,
			      int p, 
			      const double kappa[3][3], const double dkappa[3][3],
			      const double e[6][3],     const double de[6][3],
			      double* d)
{ 
  const int cNumNodes = 8;
  const int cDim      = 4;
  const int cNumDofs  = numDofs();
  const int dofX=0, dofY=1, dofZ=2, dofT=3;
  const int ls[cNumNodes][cDim] = {{0,1,2,3},{4,5,6,7},
				   {8,9,10,11},{12,13,14,15},
				   {16,17,18,19},{20,21,22,23},
				   {24,25,26,27},{28,29,30,31}};
  double xi, eta, emu, det, w, weight;
  bool cvar = false;
  double q[cNumNodes], qx[cNumNodes], qy[cNumNodes], qz[cNumNodes];

  for(int i=0; i<cNumNodes && !cvar; ++i) 
    { cvar |= (dx[i] != 0) | (dy[i] != 0) | (dz[i] != 0); }
  
  std::fill(d, d+numDofs()*numDofs(), 0);

  if(!cvar)
    {
      for(int k=1; k<=p; ++k)
	{
	  for(int l=1; l<=p; ++l)
	    {
	      for(int m=1; m<=p; ++m)
		{
		  _FORTRAN(hxgaus)(&p,&k,&p,&l,&p,&m,&xi,&eta,&emu,&weight);
		  _FORTRAN(h8shpe)(&xi,&eta,&emu,x,y,z,q,qx,qy,qz,&det);
		  if(det==0)
		    { 
		      fprintf(stderr, "Zero Jacobian determinant in PiezoBrick8\n");
		      exit(-1);
		    }
		  w = weight*det;
		  // calculate the piezo-piezo part
		  for(int j=0; j<cNumNodes; ++j)
		    {
		      int jT = ls[j][dofT];
		      double dkappa_x = (dkappa[0][0]*qx[j]+dkappa[0][1]*qy[j]+dkappa[0][2]*qz[j])*w;
		      double dkappa_y = (dkappa[1][0]*qx[j]+dkappa[1][1]*qy[j]+dkappa[1][2]*qz[j])*w;
		      double dkappa_z = (dkappa[2][0]*qx[j]+dkappa[2][1]*qy[j]+dkappa[2][2]*qz[j])*w;
		      for(int i=j; i<cNumNodes; ++i)
			{
			  int iT = ls[i][dofT];
			  d[iT*cNumDofs+jT] -= dkappa_x*qx[i] + dkappa_y*qy[i] + dkappa_z*qz[i];
			  d[jT*cNumDofs+iT] = d[iT*cNumDofs+jT];
			}
		    }
		  // calculate the piezo-structure part
		  for(int j=0; j<cNumNodes; ++j)
		    {
		      int jx = ls[j][dofX];
		      int jy = ls[j][dofY];
		      int jz = ls[j][dofZ];
		      int jT = ls[j][dofT];
		  
		      double de_0x = (de[0][0]*qx[j]+de[3][0]*qy[j]+de[5][0]*qz[j])*w;
		      double de_0y = (de[3][0]*qx[j]+de[1][0]*qy[j]+de[4][0]*qz[j])*w;
		      double de_0z = (de[5][0]*qx[j]+de[4][0]*qy[j]+de[2][0]*qz[j])*w;
		      
		      double de_1x = (de[0][1]*qx[j]+de[3][1]*qy[j]+de[5][1]*qz[j])*w;
		      double de_1y = (de[3][1]*qx[j]+de[1][1]*qy[j]+de[4][1]*qz[j])*w;
		      double de_1z = (de[5][1]*qx[j]+de[4][1]*qy[j]+de[2][1]*qz[j])*w;
		      
		      double de_2x = (de[0][2]*qx[j]+de[3][2]*qy[j]+de[5][2]*qz[j])*w;
		      double de_2y = (de[3][2]*qx[j]+de[1][2]*qy[j]+de[4][2]*qz[j])*w;
		      double de_2z = (de[5][2]*qx[j]+de[4][2]*qy[j]+de[2][2]*qz[j])*w;
		      
		      for(int i=0; i<cNumNodes; ++i)
			{
			  int ix = ls[i][dofX];
			  int iy = ls[i][dofY];
			  int iz = ls[i][dofZ];
			  int iT = ls[i][dofT];
			  
			  d[jx*cNumDofs+iT] += de_0x*qx[i]+de_1x*qy[i]+de_2x*qz[i];
			  d[jy*cNumDofs+iT] += de_0y*qx[i]+de_1y*qy[i]+de_2y*qz[i];
			  d[jz*cNumDofs+iT] += de_0z*qx[i]+de_1z*qy[i]+de_2z*qz[i];
			  
			  d[iT*cNumDofs+jx] = d[jx*cNumDofs+iT];
			  d[iT*cNumDofs+jy] = d[jy*cNumDofs+iT];
			  d[iT*cNumDofs+jz] = d[jz*cNumDofs+iT];
			}
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
int PiezoBrick8::chkOptInf(CoordSet& dcs)
{
  if(cCoefs == 0 || getDCCoefs() == 0) { return 0; }

  double dX[8], dY[8], dZ[8];
  dcs.getCoordinates(nn, numNodes(), dX, dY, dZ);
  double dvd = 0;
  for(int i=0; i<numNodes(); ++i) { dvd += dX[i]*dX[i] + dY[i]*dY[i] + dZ[i]*dZ[i]; }
  double dat = 0;
  for(int j=0; j<6; ++j) { for(int i=0; i<=j; ++i) { dat += getDCCoefs()[i   + j*6]*getDCCoefs()[i   + j*6]; } }
  
  if ( dvd > 0 && dat >0 ) return Optvar::mixedvar;
  if ( dvd > 0 )           return Optvar::coordinate;
  if ( dat > 0 )           return Optvar::attribute;

  return 0;
}

//------------------------------------------------------------------------------
void PiezoBrick8::gradstiffness(CoordSet& cs,CoordSet& dcs,FullSquareMatrix& dd, int flg)
{
  // [kappa] is a 2nd rank dielectric tensor
  double kappa[3][3], dkappa[3][3];
  // [e] is a 3rd rank tensor of piezoelectricity
  double e[6][3], de[6][3];
  fillPiezoTensors(e, kappa);
  fillDPiezoTensors(de, dkappa);

  double  X[8],  Y[8],  Z[8];
  double dX[8], dY[8], dZ[8];

   cs.getCoordinates(nn, numNodes(),  X,  Y,  Z);
  dcs.getCoordinates(nn, numNodes(), dX, dY, dZ);

  // order of gaussian quadrature
  int p = 2;

  // now we are ready to evaluate the stiffness mx
  piezo_dbr8v(X, Y, Z, dX, dY, dZ, p, kappa, dkappa, e, de, dd.data());
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
void PiezoBrick8::fillPiezoTensors(double e[6][3], double kappa[3][3])
{  
  assert(cCoefs != 0);

  for(int i=0; i<3; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  kappa[i][j] = 0.0;
	}
    }

  e[0][0] = cCoefs[0   + 0*6];
  e[0][1] = cCoefs[0   + 1*6];
  e[0][2] = cCoefs[0   + 2*6];

  e[1][0] = cCoefs[0   + 3*6];
  e[1][1] = cCoefs[0   + 4*6];
  e[1][2] = cCoefs[0   + 5*6];

  e[2][0] = cCoefs[1   + 1*6];
  e[2][1] = cCoefs[1   + 2*6];
  e[2][2] = cCoefs[1   + 3*6];

  e[3][0] = cCoefs[1   + 4*6];
  e[3][1] = cCoefs[1   + 5*6];
  e[3][2] = cCoefs[2   + 2*6];

  e[4][0] = cCoefs[2   + 3*6];
  e[4][1] = cCoefs[2   + 4*6];
  e[4][2] = cCoefs[2   + 5*6];

  e[5][0] = cCoefs[3   + 3*6];
  e[5][1] = cCoefs[3   + 4*6];
  e[5][2] = cCoefs[3   + 5*6];

  kappa[0][0] = cCoefs[4   + 4*6];
  kappa[1][1] = cCoefs[4   + 5*6];
  kappa[2][2] = cCoefs[5   + 5*6];

  if(cFrame != 0)
    {
      double eOut[6][3];
      double kappaOut[3][3];
      // rotate tensors
      rotateConstitutiveMatrix33(kappa, cFrame, kappaOut);
      rotateConstitutiveMatrix63(e,     cFrame,     eOut);
      for(int i=0; i<3; ++i) { for(int j=0; j<3; ++j) { kappa[i][j] = kappaOut[i][j]; } }
      for(int i=0; i<6; ++i) { for(int j=0; j<3; ++j) { e[i][j] = eOut[i][j]; } }
    }
  return;
}

void PiezoBrick8::fillDPiezoTensors(double de[6][3], double dkappa[3][3])
{
  assert(getDCCoefs() != 0);

  for(int i=0; i<3; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  dkappa[i][j] = 0.0;
	}
    }

  de[0][0] = getDCCoefs()[0   + 0*6];
  de[0][1] = getDCCoefs()[0   + 1*6];
  de[0][2] = getDCCoefs()[0   + 2*6];

  de[1][0] = getDCCoefs()[0   + 3*6];
  de[1][1] = getDCCoefs()[0   + 4*6];
  de[1][2] = getDCCoefs()[0   + 5*6];

  de[2][0] = getDCCoefs()[1   + 1*6];
  de[2][1] = getDCCoefs()[1   + 2*6];
  de[2][2] = getDCCoefs()[1   + 3*6];

  de[3][0] = getDCCoefs()[1   + 4*6];
  de[3][1] = getDCCoefs()[1   + 5*6];
  de[3][2] = getDCCoefs()[2   + 2*6];

  de[4][0] = getDCCoefs()[2   + 3*6];
  de[4][1] = getDCCoefs()[2   + 4*6];
  de[4][2] = getDCCoefs()[2   + 5*6];

  de[5][0] = getDCCoefs()[3   + 3*6];
  de[5][1] = getDCCoefs()[3   + 4*6];
  de[5][2] = getDCCoefs()[3   + 5*6];

  dkappa[0][0] = getDCCoefs()[4   + 4*6];
  dkappa[1][1] = getDCCoefs()[4   + 5*6];
  dkappa[2][2] = getDCCoefs()[5   + 5*6];

  if(cFrame != 0)
    {
      double eOut[6][3];
      double kappaOut[3][3];
      // rotate tensors
      rotateConstitutiveMatrix33(dkappa, cFrame, kappaOut);
      rotateConstitutiveMatrix63(de,     cFrame,     eOut);
      for(int i=0; i<3; ++i) { for(int j=0; j<3; ++j) { dkappa[i][j] = kappaOut[i][j]; } }
      for(int i=0; i<6; ++i) { for(int j=0; j<3; ++j) { de[i][j] = eOut[i][j]; } }
    }
  return;
}

#endif
