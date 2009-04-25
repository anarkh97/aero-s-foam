#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/PMLWave.d/PMLWaveElement.h>
#include <Structopt.d/Element_opt.d/PMLWave.d/PMLWaveElement2D.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Structopt.d/Optvar.h>

#include <algorithm>
#include <cmath>

//----------------------------------------------------------------------
Element* PMLWaveElement2D::clone()
{
  return new PMLWaveElement2D(*this);
}

//----------------------------------------------------------------------
extern "C"
{
  void _FORTRAN(zgetcmt)(double &, double&, double&, DComplex&, DComplex&, DComplex*);
  void _FORTRAN(zquad4m)(double*, double*, double*, DComplex*, const int&,
			 DComplex*, const int&);
}

#define PMLWAVE2D_NDIM 2
#define PMLWAVE2D_NNODES 4
//----------------------------------------------------------------------
FullSquareMatrixC PMLWaveElement2D::complexStiffness(CoordSet &cs, DComplex *d, int flg)
{
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);
  Node &nd4 = cs.getNode(nn[3]);
  
  DComplex c[PMLWAVE2D_NDIM*PMLWAVE2D_NDIM*PMLWAVE2D_NDIM*PMLWAVE2D_NDIM];
  double x[PMLWAVE2D_NNODES], y[PMLWAVE2D_NNODES], h[PMLWAVE2D_NNODES];

  x[0] = nd1.x; y[0] = nd1.y; 
  x[1] = nd2.x; y[1] = nd2.y;
  x[2] = nd3.x; y[2] = nd3.y; 
  x[3] = nd4.x; y[3] = nd4.y;

  h[0] = h[1] = h[2] = h[3] = prop->eh;

  DComplex epsx = eps_x(prop);
  DComplex epsy = eps_y(prop);
  _FORTRAN(zgetcmt)(prop->A, prop->E, prop->nu, epsx, epsy, c);

  const int numgauss = 2;
  const int numdof   = numDofs();

  std::fill(d, d+numdof*numdof, 0.0);
  _FORTRAN(zquad4m)(x, y, h, c, numgauss, d, numdof);
  return FullSquareMatrixC (numdof, d);
}

//----------------------------------------------------------------------
FullSquareMatrixC PMLWaveElement2D::complexMassMatrix(CoordSet &cs, DComplex *mel, double mratio)
{
  const DComplex epsx = eps_x(prop);
  const DComplex epsy = eps_y(prop);
  const DComplex epsxy = epsx*epsy;
  const int numdof   = numDofs();

  Element::massMatrix(cs, reinterpret_cast<double*>(mel), mratio);
  for(int i=numdof*numdof-1; i>=0; --i)
    { mel[i] = epsxy*reinterpret_cast<double*>(mel)[i]; }
  return FullSquareMatrixC(numdof, mel);
}

#endif
