#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/PMLWave.d/PMLWaveElement.h>
#include <Structopt.d/Element_opt.d/PMLWave.d/PMLWaveElement3D.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Structopt.d/Optvar.h>

#include <algorithm>
#include <cmath>

//----------------------------------------------------------------------
Element* PMLWaveElement3D::clone()
{
  return new PMLWaveElement3D(*this);
}

//----------------------------------------------------------------------
extern "C"
{
  void  _FORTRAN(zbrkcmt)(double&, double&, DComplex&, DComplex&, DComplex&, DComplex*);  
  void  _FORTRAN(zbrik8v)(double*, double*, double*, DComplex*, const int&,
			  DComplex*, int &);
}

#define PMLWAVE2D_NDIM 3
#define PMLWAVE2D_NNODES 8
//----------------------------------------------------------------------
FullSquareMatrixC PMLWaveElement3D::complexStiffness(CoordSet &cs, DComplex *d, int flg)
{
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);
  Node &nd5 = cs.getNode(nn[4]);
  Node &nd6 = cs.getNode(nn[5]);
  Node &nd7 = cs.getNode(nn[6]);
  Node &nd8 = cs.getNode(nn[7]);
  
  double x[PMLWAVE3D_NNODES], y[PMLWAVE3D_NNODES], z[PMLWAVE3D_NNODES];
  DComplex c[PMLWAVE3D_NDIM*PMLWAVE3D_NDIM*PMLWAVE3D_NDIM*PMLWAVE3D_NDIM];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z; 
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
  x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
  x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
  x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
  x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;
	
  DComplex epsx = eps_x(prop);
  DComplex epsy = eps_y(prop);
  DComplex epsz = eps_z(prop);

  _FORTRAN(zbrkcmt)(prop->E, prop->nu, epsx, epsy, epsz, c);

  const int numgauss = 2;
  const int numdof   = numDofs();
  int status;
  
  _FORTRAN(zbrik8v)(x, y, z, c, numgauss, d, status);

  FullSquareMatrixC ret(numdof, d);
  return ret;
}

//----------------------------------------------------------------------
FullSquareMatrixC PMLWaveElement3D::complexMassMatrix(CoordSet &cs, DComplex *mel, double mratio)
{
  const DComplex epsx = eps_x(prop);
  const DComplex epsy = eps_y(prop);
  const DComplex epsz = eps_z(prop);
  const DComplex epsxyz = epsx*epsy*epsz;
  const int numdof   = numDofs();

  FullSquareMatrix tmp(numdof, reinterpret_cast<double*>(mel));
  EightNodeBrick::massMatrix(cs, tmp.data());
  Element::lumpMatrix(tmp, mratio);
  
  for(int i=numdof*numdof-1; i>=0; --i)
    { mel[i] = epsxyz*reinterpret_cast<double*>(mel)[i]; }
  return FullSquareMatrixC(numdof, mel);
}

#endif
