#ifdef USE_EIGEN3
#include <Utils.d/dofset.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <iostream>
#include <Element.d/Function.d/VectorValuedFunction.h>
#include <Element.d/Function.d/SpaceJacobian.h>
#include <Element.d/Function.d/utilities.hpp>
#include <unsupported/Eigen/NumericalDiff>

template<template <typename S> class VectorValuedFunctionTemplate>
ForceFunctionElement<VectorValuedFunctionTemplate>
::ForceFunctionElement(int _nNodes, DofSet nodalDofs, int* _nn)
 : BoundaryElement(_nNodes, nodalDofs, _nn)
{
  SO3param = 2;
}

template<template <typename S> class VectorValuedFunctionTemplate>
ForceFunctionElement<VectorValuedFunctionTemplate>
::ForceFunctionElement(int _nNodes, DofSet *nodalDofs, int* _nn)
 : BoundaryElement(_nNodes, nodalDofs, _nn)
{
  SO3param = 2;
}

template<template <typename S> class VectorValuedFunctionTemplate>
ForceFunctionElement<VectorValuedFunctionTemplate>
::ForceFunctionElement(int _nNodes, DofSet *nodalInputDofs, DofSet *nodalOutputDofs, int* _nn)
 : BoundaryElement(_nNodes, nodalInputDofs, nodalOutputDofs, _nn)
{
  SO3param = 2;
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
ForceFunctionElement<VectorValuedFunctionTemplate>
::getInputs(Eigen::Matrix<double,VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q,
            CoordSet& c0, GeomState *curState, GeomState *refState)
{
  // prepare the constraint function inputs
  if(curState == NULL) { // in this case the function will be evaluated in the undeformed configuration
    q.setZero();
    //q.setConstant(std::numeric_limits<double>::epsilon()); // XXXX
    //q.setConstant(1e-8);
  }
  else {
    switch(SO3param) {
      case 0 : // total lagrangian description of rotations
        for (int k = 0; k < inputs.size(); k++) {
          int i = inputs[k];
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
              Eigen::Vector3d theta;
              mat_to_vec((*curState)[terms[i].nnum].R, theta.data());
              q[k] = theta[terms[i].dofnum-3];
            } break;
          }
        }
        break;
      case 1 : // updated lagrangian description of spatial rotations
        for (int k = 0; k < inputs.size(); k++) {
          int i = inputs[k];
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
        }
        break;
      case 3 :
        // TODO
        break;
      case 2: case 4: case 5: default : // eulerian description of rotations 
        for (int k = 0; k < inputs.size(); k++) {
          int i = inputs[k];
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
              q[k] = 0; //10*std::numeric_limits<double>::epsilon(); //XXXX 0; // spin
            } break;
          }
        }
        break;
    }
  }
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
ForceFunctionElement<VectorValuedFunctionTemplate>::computePressureForce(CoordSet& c0, Vector& F,
                                                                         GeomState *c1, int, double t)
{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst, (GeomState*)NULL, t);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, c1, NULL);

  // evaluate the function and store values terms
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  Eigen::Matrix<double,M,1> fval = f(q,t);

  // fill F
  for(int i=0; i<numDofs(); ++i) F[i] = 0;

  for(int i=0; i<M; ++i) {
    int k = outputs[i];
    F[k] = fval[i];
  }
}

template<template <typename S> class VectorValuedFunctionTemplate>
FullSquareMatrix
ForceFunctionElement<VectorValuedFunctionTemplate>::stiffness(CoordSet& c0, double* karray, int)
{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q; // = Eigen::Matrix<double,N,1>::Zero();
  getInputs(q, c0, NULL, NULL);
  double t = 0;

  // evaluate the jacobian (partial derivatives w.r.t. the spatial variables)
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  FullM J(M,N);
  GeomState c1(c0);
  getJacobian(NULL, c1, c0, J, t);

  // fill K
  FullSquareMatrix K(numDofs(), karray);
  K.zero();
  for(int i=0; i<M; ++i) {
    int k = outputs[i];
    for(int j=0; j<N; ++j) {
      int l = inputs[j];
      K[k][l] = J[i][j];
    }
  }
  return K;
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
ForceFunctionElement<VectorValuedFunctionTemplate>::getInternalForce(GeomState *refState, GeomState& c1, CoordSet& c0, FullSquareMatrix&,
                                                                     double* F, double, double t)
{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  GeomState *cx;
  if(SO3param == 2 || SO3param == 4 || SO3param == 5) cx = &c1; // current state for eulerian
  else if(SO3param == 1 || SO3param == 3) cx = refState; // reference state for update lagrangian
  else cx = NULL; // not used for total lagrangian
  getConstants(c0, sconst, iconst, cx, t);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q; // = Eigen::Matrix<double,N,1>::Zero();
  getInputs(q, c0, &c1, refState);

  // evaluate the function and store values terms
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  Eigen::Matrix<double,M,1> fval = f(q,t);

  // fill F
  for(int i=0; i<numDofs(); ++i) F[i] = 0;

  for(int i=0; i<M; ++i) {
    int k = outputs[i];
    F[k] = fval[i];
  }
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
ForceFunctionElement<VectorValuedFunctionTemplate>::getStiffAndForce(GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan,
                                                                     double* F, double dt, double t)

{
  getStiffAndForce((GeomState*)NULL, c1, c0, Ktan, F, dt, t);
}

template<template <typename S> class VectorValuedFunctionTemplate>
void 
ForceFunctionElement<VectorValuedFunctionTemplate>::getStiffAndForce(GeomState *refState, GeomState& c1, CoordSet& c0, FullSquareMatrix& Ktan,
                                                                     double* F, double, double t)

{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  GeomState *cx;
  if(SO3param == 2 || SO3param == 4 || SO3param == 5) cx = &c1; // current state for eulerian
  else if(SO3param == 1 || SO3param == 3) cx = refState; // reference state for update lagrangian
  else cx = NULL; // not used for total lagrangian
  getConstants(c0, sconst, iconst, cx, t);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q; // = Eigen::Matrix<double,N,1>::Zero();
  getInputs(q, c0, &c1, refState);

  // evaluate the function and store values terms
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  Eigen::Matrix<double,M,1> fval = f(q,t);
  //std::cerr << "fval = \n" << fval.transpose() << std::endl;

  // evaluate the jacobian (partial derivatives w.r.t. the spatial variables)
  FullM J(M,N);
  getJacobian(refState, c1, c0, J, t);

  // fill Ktan and F
  Ktan.zero();
  for(int i=0; i<numDofs(); ++i) F[i] = 0;

  for(int i=0; i<M; ++i) {
    int k = outputs[i];
    F[k] = fval[i];
    for(int j=0; j<N; ++j) {
      int l = inputs[j];
      Ktan[k][l] = J[i][j];
    }
  }
}

template<template <typename S> class VectorValuedFunctionTemplate>
void 
ForceFunctionElement<VectorValuedFunctionTemplate>::getJacobian(GeomState *refState, GeomState &c1, CoordSet& c0, FullM& B, double t) 
{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  GeomState *cx;
  if(SO3param == 2 || SO3param == 4 || SO3param == 5) cx = &c1; // current state for eulerian
  else if(SO3param == 1 || SO3param == 3) cx = refState; // reference state for update lagrangian
  else cx = NULL; // not used for total lagrangian
  getConstants(c0, sconst, iconst, cx, t);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, &c1, refState);

  // instantiate the jacobian object
  Simo::SpaceJacobian<double,VectorValuedFunctionTemplate> dfdq(sconst,iconst);

  // evaluate the jacobian
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  Eigen::Matrix<double,M,N> J;
  J = dfdq(q, t);
  //std::cerr << "here are the eigenvalues of J\n" << J.template block<M,M>(0,0).eigenvalues().transpose() << std::endl;
  //std::cerr << "J = \n" << J << std::endl;

/* TODO
  Eigen::Matrix<double,N,N> H;
  int flag = prop->constraint_hess;
  switch (flag) {
   
    default : case 0 :
      H.setZero();
      break;
    case 1: {
      // instantiate the constraint hessian object
      SacadoReverseJacobian<ConstraintJacobian<double,VectorValuedFunctionTemplate>,true> d2fdq2(dfdq);
      // evaluate the constraint hessian
      d2fdq2(q, H);
    } break;
#ifdef USE_SACADO
    case 2: {
      // instantiate the constraint hessian object
      SacadoReverseJacobian<ConstraintJacobian<double,VectorValuedFunctionTemplate>,false> d2fdq2(dfdq);
      // evaluate the constraint hessian
      d2fdq2(q, H);
    } break;
#endif
    case 4 : {
      // instantiate the forward difference object
      Eigen::NumericalDiff<ConstraintJacobian<double,VectorValuedFunctionTemplate>,Eigen::Forward> fd(dfdq, prop->constraint_hess_eps);
      // evaluate the forward difference approximation to the hessian
      fd.df(q, H);
    } break;
    case 5 : {
      // instantiate the central difference object
      Eigen::NumericalDiff<ConstraintJacobian<double,VectorValuedFunctionTemplate>,Eigen::Central> cd(dfdq, prop->constraint_hess_eps);
      // evaluate the central difference approximation to the hessian
      cd.df(q, H);
    } break;
  }
*/

  for(int i = 0; i < M; ++i)
    for(int j = 0; j < N; ++j) {
      B[i][j] = J(i,j);
    }

}

#endif
