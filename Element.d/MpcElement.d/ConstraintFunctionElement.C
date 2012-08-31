#ifdef USE_EIGEN3
#include <Utils.d/dofset.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Element.d/MpcElement.d/ConstraintFunction.d/ConstraintFunction.h>
#include <iostream>
#include <unsupported/Eigen/NumericalDiff>

template<template <typename S> class ConstraintFunctionTemplate>
ConstraintFunctionElement<ConstraintFunctionTemplate>
::ConstraintFunctionElement(int _nNodes, DofSet nodalDofs, int* _nn, int _type, int _rotdescr)
 : MpcElement(_nNodes, nodalDofs, _nn)
{
  type = _type;
  rotdescr = _rotdescr;
}

template<template <typename S> class ConstraintFunctionTemplate>
ConstraintFunctionElement<ConstraintFunctionTemplate>
::ConstraintFunctionElement(int _nNodes, DofSet *nodalDofs, int* _nn, int _type, int _rotdescr)
 : MpcElement(_nNodes, nodalDofs, _nn)
{
  type = _type;
  rotdescr = _rotdescr;
}

template<template <typename S> class ConstraintFunctionTemplate>
void
ConstraintFunctionElement<ConstraintFunctionTemplate>
::getInputs(Eigen::Matrix<double,ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q,
            CoordSet& c0, GeomState *curState, GeomState *refState)
{
  // prepare the constraint function inputs
  int k = 0;
  if(curState == NULL) { // in this case the function will be evaluated in the undeformed configuration
    for (int i = 0; i < nterms; i++) {
      switch(terms[i].dofnum) {
        case 0 :
          q[k] = 0;
          break;
        case 1 :
          q[k] = 0;
          break;
        case 2 :
          q[k] = 0;
          break;
        case 3 : case 4 : case 5 : {
          q[k] = 0;
        } break;
      }
      k++;
    }
  }
  else if(rotdescr == 0) { // total lagrangian description of rotations
    for (int i = 0; i < nterms; i++) {
      switch(terms[i].dofnum) {
        case 0 :
          q[k] = (*curState)[terms[i].nnum].x - c0[terms[i].nnum]->x;
          break;
        case 1 :
          q[k] = (*curState)[terms[i].nnum].y - c0[terms[i].nnum]->y;
          break;
        case 2 :
          q[k] = (*curState)[terms[i].nnum].z - c0[terms[i].nnum]->z;
          break;
        case 3 : case 4 : case 5 : {
          double theta[3];
          mat_to_vec((*curState)[terms[i].nnum].R, theta);
          q[k] = theta[terms[i].dofnum-3];
        } break;
      }
      k++;
    }
  }
  else if(rotdescr == 1) { // updated lagrangian description of spatial rotations
    for (int i = 0; i < nterms; i++) {
      switch(terms[i].dofnum) {
        case 0 :
          q[k] = (*curState)[terms[i].nnum].x - c0[terms[i].nnum]->x;
          break;
        case 1 :
          q[k] = (*curState)[terms[i].nnum].y - c0[terms[i].nnum]->y;
          break;
        case 2 :
          q[k] = (*curState)[terms[i].nnum].z - c0[terms[i].nnum]->z;
          break;
        case 3 : case 4 : case 5 : {
          // spatial incremental rotation matrix and vector
          double dR[3][3], dtheta[3];
          // curState = dR*refState  --> dR = curState*refState^T
          mat_mult_mat((*curState)[terms[i].nnum].R, (*refState)[terms[i].nnum].R, dR, 2);
          mat_to_vec(dR, dtheta);
          q[k] = dtheta[terms[i].dofnum-3];
        } break;
      }
      k++;
    }
  }
  else { // eulerian description of rotations
    for (int i = 0; i < nterms; i++) {
      switch(terms[i].dofnum) {
        case 0 :
          q[k] = (*curState)[terms[i].nnum].x - c0[terms[i].nnum]->x;
          break;
        case 1 :
          q[k] = (*curState)[terms[i].nnum].y - c0[terms[i].nnum]->y;
          break;
        case 2 :
          q[k] = (*curState)[terms[i].nnum].z - c0[terms[i].nnum]->z;
          break;
        case 3 : case 4 : case 5 : {
          q[k] = 0; // spin
        } break;
      }
      k++;
    }
  }
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
  Eigen::Matrix<double,N,1> q; // = Eigen::Matrix<double,N,1>::Zero();
  getInputs(q, c0, NULL, NULL);
  double t = 0;

  // evaluate the constraint function and store -ve value in LMPCons::rhs
  original_rhs.r_value = rhs.r_value = -f(q,t);

  // instantiate the constraint jacobian object
  ConstraintJacobian<double,ConstraintFunctionTemplate> dfdq(sconst,iconst,t);

  // evaluate the constraint jacobian (partial derivatives w.r.t. the spatial variables)
  // and store coefficients in LMPCons::terms array
  Eigen::Matrix<double,N,1> J;
  dfdq(q, &J);
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
  getConstants(c0, sconst, iconst, &c1);
  ConstraintFunctionTemplate<double> f(sconst,iconst);

  // prepare the constraint function inputs
  const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, &c1, NULL);

  // evaluate the constraint function and store -ve value in LMPCons::rhs
  rhs.r_value = -f(q,t);

  // instantiate the constraint jacobian object
  ConstraintJacobian<double,ConstraintFunctionTemplate> dfdq(sconst,iconst,t);

  // evaluate the constraint jacobian (partial derivatives w.r.t. the spatial variables)
  // and store coefficients in LMPCons::terms array
  Eigen::Matrix<double,N,1> J;
  dfdq(q, &J);
  for(int i = 0; i < nterms; ++i) terms[i].coef.r_value = J[i];
}

template<template <typename S> class ConstraintFunctionTemplate>
void 
ConstraintFunctionElement<ConstraintFunctionTemplate>::getHessian(GeomState &c1, CoordSet& c0, FullSquareMatrix& B, double t) 
{
  // instantiate the constraint function object
  Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst, &c1);
  ConstraintFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, &c1, NULL);

  // instantiate the constraint jacobian object
  ConstraintJacobian<double,ConstraintFunctionTemplate> dfdq(sconst,iconst,t);

  Eigen::Matrix<double,N,N> H;
  int flag = prop->constraint_hess;
  switch (flag) {
   
    default : case 0 :
      H.setZero();
      break;
    case 1: {
      // instantiate the constraint hessian object
      SacadoReverseJacobian<ConstraintJacobian<double,ConstraintFunctionTemplate>,true> d2fdq2(dfdq);
      // evaluate the constraint hessian
      d2fdq2(q, H);
    } break;
#ifdef USE_SACADO
    case 2: {
      // instantiate the constraint hessian object
      SacadoReverseJacobian<ConstraintJacobian<double,ConstraintFunctionTemplate>,false> d2fdq2(dfdq);
      // evaluate the constraint hessian
      d2fdq2(q, H);
    } break;
#endif
    case 4 : {
      // instantiate the forward difference object
      Eigen::NumericalDiff<ConstraintJacobian<double,ConstraintFunctionTemplate>,Eigen::Forward> fd(dfdq, prop->constraint_hess_eps);
      // evaluate the forward difference approximation to the hessian
      fd.df(q, H);
    } break;
    case 5 : {
      // instantiate the central difference object
      Eigen::NumericalDiff<ConstraintJacobian<double,ConstraintFunctionTemplate>,Eigen::Central> cd(dfdq, prop->constraint_hess_eps);
      // evaluate the central difference approximation to the hessian
      cd.df(q, H);
    } break;
  }
  for(int i = 0; i < nterms; ++i)
    for(int j = 0; j < nterms; ++j)
      B[i][j] = H(i,j);
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
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, NULL, NULL);

  // evaluate the function
  rhs.r_value = -f(q,t);

  MpcElement::computePressureForce(c0, elPressureForce, c1, cflg, t);
}

template<template <typename S> class ConstraintFunctionTemplate>
double
ConstraintFunctionElement<ConstraintFunctionTemplate>::getVelocityConstraintRhs(GeomState *c1, CoordSet& c0,
                                                                                double t)
{
  // instantiate the constraint function object
  Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst, c1);
  ConstraintFunctionTemplate<double> f(sconst,iconst);

  // prepare the constraint function inputs
  const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, c1, NULL);

  // instantiate the first time derivative of the constraint function object
  PartialTimeDerivative<double,ConstraintFunctionTemplate> dfdt(sconst,iconst,q);

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
  // instantiate the constraint function object
  Eigen::Array<double, ConstraintFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, ConstraintFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst, c1);
  ConstraintFunctionTemplate<double> f(sconst,iconst);

  // prepare the constraint function inputs
  const int N = ConstraintFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, c1, NULL);

  // instantiate the partial time derivative of the constraint function object
  PartialTimeDerivative<double,ConstraintFunctionTemplate> dfdt(sconst,iconst,q);

  // instantiate the second partial time derivative of the constraint function object
  SacadoReverseJacobian<PartialTimeDerivative<double,ConstraintFunctionTemplate> > d2fdt2(dfdt);

  // evaluate the second partial time derivative of the constraint function
  Eigen::Matrix<double,1,1> x, y;
  x[0] = t;
  d2fdt2(x,y);

  // instantiate the partial time derivative of the constraint jacobian object
  TemporalViewOfConstraintJacobian<double,ConstraintFunctionTemplate> J(sconst,iconst,q);
  SacadoReverseJacobian<TemporalViewOfConstraintJacobian<double,ConstraintFunctionTemplate> > dJdt(J);

  // evaluate the partial time derivative of the constraint jacobian
  Eigen::Matrix<double,N,1> z;
  dJdt(x,z);

  // instantiate the jacobian of the constraint jacobian velocity product
  Eigen::Matrix<double,N,1> v;
  for(int i = 0; i < nterms; ++i) v[i] = (*c1)[terms[i].nnum].v[terms[i].dofnum];
  ConstraintJacobianVelocityProduct<double,ConstraintFunctionTemplate> Jv(sconst,iconst,t,v); // assuming v is constant for now
  SacadoReverseJacobian<ConstraintJacobianVelocityProduct<double,ConstraintFunctionTemplate> > dJvdq(Jv);

  // evaluate the jacobian of the constraint jacobian velocity product
  Eigen::Matrix<double,1,N> w;
  dJvdq(q,w);

  return -w.dot(v) - 2*z.dot(v) -y[0];
#else
  return 0;
#endif
}
#endif
