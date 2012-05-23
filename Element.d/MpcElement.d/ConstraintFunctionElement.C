#include <Utils.d/dofset.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Element.d/MpcElement.d/ConstraintFunction.h>

template<template <typename S> class ConstraintFunctionTemplate>
ConstraintFunctionElement<ConstraintFunctionTemplate>::ConstraintFunctionElement(int _nNodes, DofSet nodalDofs, int* _nn, int _type)
 : MpcElement(_nNodes, nodalDofs, _nn)
{
  type = _type;
}

template<template <typename S> class ConstraintFunctionTemplate>
ConstraintFunctionElement<ConstraintFunctionTemplate>::ConstraintFunctionElement(int _nNodes, DofSet *nodalDofs, int* _nn, int _type)
 : MpcElement(_nNodes, nodalDofs, _nn)
{
  type = _type;
}

template<template <typename S> class ConstraintFunctionTemplate>
void
ConstraintFunctionElement<ConstraintFunctionTemplate>::buildFrame(CoordSet& c0)
{
  // instantiate the constraint function object
  Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst);
  ConstraintFunctionTemplate<double> f(sconst,iconst);

  // prepare the constraint function inputs
  const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> u = Eigen::Matrix<double,N,1>::Zero();
  double t = 0;

  // evaluate the constraint function and store -ve value in LMPCons::rhs
  rhs.r_value = -f(u,t);

  // instantiate the constraint jacobian object
  SpatialJacobian<double,ConstraintFunctionTemplate> dfdu(sconst,iconst,t);

  // evaluate the constraint jacobian (partial derivatives w.r.t. the spatial variables)
  // and store coefficients in LMPCons::terms array
  Eigen::Matrix<double,N,1> J;
  dfdu(u, &J);
  for(int i = 0; i < nterms; ++i) terms[i].coef.r_value = J[i];

  // TODO: if the function is quadratic then we should compute and store the hessian for future reference
}

template<template <typename S> class ConstraintFunctionTemplate>
void 
ConstraintFunctionElement<ConstraintFunctionTemplate>::update(GeomState& c1, CoordSet& c0, double t) 
{
  // instantiate the constraint function object
  Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst);
  ConstraintFunctionTemplate<double> f(sconst,iconst);

  // prepare the constraint function inputs
  const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> u;
  int k = 0;
  for (int i = 0; i < nterms; i++) {
    switch(terms[i].dofnum) {
      case 0 :
        u[k] = c1[terms[i].nnum].x - c0[terms[i].nnum]->x;
        break;
      case 1 : 
        u[k] = c1[terms[i].nnum].y - c0[terms[i].nnum]->y;
        break;
      case 2 :
        u[k] = c1[terms[i].nnum].z - c0[terms[i].nnum]->z;
        break;
      case 3 : case 4 : case 5 : {
        double theta[3];
        mat_to_vec(c1[terms[i].nnum].R, theta);
        u[k] = theta[terms[i].dofnum-3];
      } break;
    }
    k++;
  }

  // evaluate the constraint function and store -ve value in LMPCons::rhs
  rhs.r_value = -f(u,t);

  // instantiate the constraint jacobian object
  SpatialJacobian<double,ConstraintFunctionTemplate> dfdu(sconst,iconst,t);

  // evaluate the constraint jacobian (partial derivatives w.r.t. the spatial variables)
  // and store coefficients in LMPCons::terms array
  Eigen::Matrix<double,N,1> J;
  dfdu(u, &J);
  for(int i = 0; i < nterms; ++i) terms[i].coef.r_value = J[i];
}

template<template <typename S> class ConstraintFunctionTemplate>
void 
ConstraintFunctionElement<ConstraintFunctionTemplate>::getHessian(GeomState &c1, CoordSet& c0, FullSquareMatrix& B, double t) 
{
#ifdef USE_SACADO
  // prepare the function inputs
  const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> u;
  int k = 0;
  for (int i = 0; i < nterms; i++) {
    switch(terms[i].dofnum) {
      case 0 :
        u[k] = c1[terms[i].nnum].x - c0[terms[i].nnum]->x;
        break;
      case 1 :
        u[k] = c1[terms[i].nnum].y - c0[terms[i].nnum]->y;
        break;
      case 2 :
        u[k] = c1[terms[i].nnum].z - c0[terms[i].nnum]->z;
        break;
      case 3 : case 4 : case 5 : {
        double theta[3];
        mat_to_vec(c1[terms[i].nnum].R, theta);
        u[k] = theta[terms[i].dofnum-3];
      } break;
    }
    k++;
  }

  // instantiate the constraint jacobian object
  Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst);
  SpatialJacobian<double,ConstraintFunctionTemplate> dfdu(sconst,iconst,t);

  // instantiate the constraint hessian object
  SacadoReverseJacobian<SpatialJacobian<double,ConstraintFunctionTemplate> > d2fdu2(dfdu);

  // evaluate the constraint hessian and store the coefficients in B
  Eigen::Matrix<double,N,N> H;
  d2fdu2(u, H);
  for(int i = 0; i < nterms; ++i)
    for(int j = 0; j < nterms; ++j)
      B[i][j] = H(i,j);
#else
  B.zero();
#endif
}

template<template <typename S> class ConstraintFunctionTemplate>
void
ConstraintFunctionElement<ConstraintFunctionTemplate>::computePressureForce(CoordSet& c0, Vector& elPressureForce,
                                                                            GeomState *c1, int cflg, double t)
{
  // instantiate the constraint function object
  Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst);
  ConstraintFunctionTemplate<double> f(sconst,iconst);

  // prepare the constraint function inputs
  const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> u = Eigen::Matrix<double,N,1>::Zero();

  // evaluate the function
  rhs.r_value = -f(u,t);

  MpcElement::computePressureForce(c0, elPressureForce, c1, cflg, t);
}

template<template <typename S> class ConstraintFunctionTemplate>
double
ConstraintFunctionElement<ConstraintFunctionTemplate>::getVelocityConstraintRhs(GeomState *c1, CoordSet& c0,
                                                                                double t)
{
  // prepare the constraint function inputs
  const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> u;
  int k = 0;
  for (int i = 0; i < nterms; i++) {
    switch(terms[i].dofnum) {
      case 0 :
        u[k] = (*c1)[terms[i].nnum].x - c0[terms[i].nnum]->x;
        break;
      case 1 :
        u[k] = (*c1)[terms[i].nnum].y - c0[terms[i].nnum]->y;
        break;
      case 2 :
        u[k] = (*c1)[terms[i].nnum].z - c0[terms[i].nnum]->z;
        break;
      case 3 : case 4 : case 5 : {
        double theta[3];
        mat_to_vec((*c1)[terms[i].nnum].R, theta);
        u[k] = theta[terms[i].dofnum-3];
      } break;
    }
    k++;
  }

  // instantiate the first time derivative of the constraint function object
  Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst);
  TemporalJacobian<double,ConstraintFunctionTemplate> dfdt(sconst,iconst,u);

  // evaluate the first time derivative of the constraint function
  Eigen::Matrix<double,1,1> x, v;
  x[0] = t;
  dfdt(x, &v);
  return -v[0];
}

template<template <typename S> class ConstraintFunctionTemplate>
double
ConstraintFunctionElement<ConstraintFunctionTemplate>::getAccelerationConstraintRhs(GeomState *c1, CoordSet& c0,
                                                                                    double t)
{
#ifdef USE_SACADO
  // prepare the constraint function inputs
  const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> u;
  int k = 0;
  for (int i = 0; i < nterms; i++) {
    switch(terms[i].dofnum) {
      case 0 :
        u[k] = (*c1)[terms[i].nnum].x - c0[terms[i].nnum]->x;
        break;
      case 1 :
        u[k] = (*c1)[terms[i].nnum].y - c0[terms[i].nnum]->y;
        break;
      case 2 :
        u[k] = (*c1)[terms[i].nnum].z - c0[terms[i].nnum]->z;                 
        break;
      case 3 : case 4 : case 5 : {
        double theta[3];
        mat_to_vec((*c1)[terms[i].nnum].R, theta);
        u[k] = theta[terms[i].dofnum-3];
      } break;
    }
    k++;
  }

  // instantiate the time derivative of the constraint function object
  Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst);
  TemporalJacobian<double,ConstraintFunctionTemplate> dfdt(sconst,iconst,u);

  // instantiate the second time derivative of the constraint function object
  SacadoReverseJacobian<TemporalJacobian<double,ConstraintFunctionTemplate> > d2fdt2(dfdt);

  // evaluate the second time derivative of the constraint function
  Eigen::Matrix<double,1,1> x, y;
  x[0] = t;
  d2fdt2(x,y);

  // instantiate the time derivative of the constraint jacobian object
  TemporalViewOfSpatialJacobian<double,ConstraintFunctionTemplate> J(sconst,iconst,u);
  SacadoReverseJacobian<TemporalViewOfSpatialJacobian<double,ConstraintFunctionTemplate> > dJdt(J);

  // evaluate the time derivative of the constraint jacobian
  Eigen::Matrix<double,N,1> z;
  dJdt(x,z);

  // instantiate the spatial derivative of the constraint jacobian velocity product
  Eigen::Matrix<double,N,1> v;
  for(int i = 0; i < nterms; ++i) v[i] = (*c1)[terms[i].nnum].v[terms[i].dofnum];
  SpatialJacobianVelocityProduct<double,ConstraintFunctionTemplate> Jv(sconst,iconst,t,v); // assuming v is constant for now
  SacadoReverseJacobian<SpatialJacobianVelocityProduct<double,ConstraintFunctionTemplate> > dJvdu(Jv);

  // evaluate the spatial derivative of the constraint jacobian velocity product
  Eigen::Matrix<double,1,N> w;
  dJvdu(u,w);

  return -w.dot(v) - 2*z.dot(v) -y[0];
#else
  return 0;
#endif
}

