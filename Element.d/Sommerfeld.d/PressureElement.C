#ifdef USE_EIGEN3
#include <Utils.d/dofset.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Element.d/Function.d/Function.h>
#include <iostream>
#include <unsupported/Eigen/NumericalDiff>

template<template <typename S> class VectorValuedFunctionTemplate>
PressureElement<VectorValuedFunctionTemplate>
::PressureElement(int _nNodes, DofSet nodalDofs, int* _nn, PressureBCond* _pbc)
 : nNodes(_nNodes), pbc(_pbc)
{
  // this constructor is for a force involving the same DofSet on each node, used for both inputs and outputs
  nn = new int[nNodes];
  for(int i = 0; i < nNodes; ++i)
    nn[i] = _nn[i];
  addTerms(nodalDofs);
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
PressureElement<VectorValuedFunctionTemplate>
::addTerms(DofSet nodalDofs)
{
  terms.clear();
  nterms = 0;
  for(int i = 0; i < nNodes; ++i) {
    int j = 0;
    for(int k = 0; ; ++k)
      if(nodalDofs.contains(1 << k)) {
        BCond t;
        t.setData(nn[i], k, 0.0);
        terms.push_back(t);
        nterms++;
        j++;
        if(j == nodalDofs.count()) break;
      }
  }
}

template<template <typename S> class VectorValuedFunctionTemplate>
PressureElement<VectorValuedFunctionTemplate>
::~PressureElement()
{
  delete [] nn;
}

template<template <typename S> class VectorValuedFunctionTemplate>
int
PressureElement<VectorValuedFunctionTemplate>
::numNodes()
{
  return nNodes;
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
PressureElement<VectorValuedFunctionTemplate>
::renum(int* table)
{
  for(int i = 0; i < numNodes(); ++i)
    if(nn[i] > -1)
      nn[i] = table[nn[i]];
  for(int i = 0; i < nterms; ++i)
    terms[i].nnum = table[terms[i].nnum];
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
PressureElement<VectorValuedFunctionTemplate>
::renum(EleRenumMap& table)
{
  for(int i = 0; i < numNodes(); ++i)
    if(nn[i] > -1)
      nn[i] = table[nn[i]];
  for(int i = 0; i < nterms; ++i)
    terms[i].nnum = table[terms[i].nnum];
}

template<template <typename S> class VectorValuedFunctionTemplate>
int*
PressureElement<VectorValuedFunctionTemplate>
::nodes(int* p)
{
  if(p == 0) p = new int[numNodes()];
  for(int i = 0; i < numNodes(); ++i) p[i] = nn[i];
  return p;
}

template<template <typename S> class VectorValuedFunctionTemplate>
int
PressureElement<VectorValuedFunctionTemplate>
::numDofs()
{
  return nterms;
}

template<template <typename S> class VectorValuedFunctionTemplate>
int *
PressureElement<VectorValuedFunctionTemplate>
::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[numDofs()];
  for(int i = 0; i < nterms; i++)
    dsa.number(terms[i].nnum, 1 << terms[i].dofnum, p+i);
  return p;
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
PressureElement<VectorValuedFunctionTemplate>
::markDofs(DofSetArray &dsa)
{
  for(int i = 0; i < nterms; i++)
    dsa.mark(terms[i].nnum, 1 << terms[i].dofnum);
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
PressureElement<VectorValuedFunctionTemplate>
::getInputs(Eigen::Matrix<double,VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates,1> &q,
            CoordSet& c0, GeomState *curState)
{
  // prepare the constraint function inputs
  if(curState == NULL) {
    // in this case the function will be evaluated in the undeformed configuration
    q.setZero();
  }
  else {
    for (int i = 0; i < terms.size(); i++) {
      switch(terms[i].dofnum) {
        case 0 :
          q[i] = (*curState)[terms[i].nnum].x - c0[terms[i].nnum]->x;
          break;
        case 1 :
          q[i] = (*curState)[terms[i].nnum].y - c0[terms[i].nnum]->y;
          break;
        case 2 :
          q[i] = (*curState)[terms[i].nnum].z - c0[terms[i].nnum]->z;
          break;
        case 3 : case 4 : case 5 : {
          q[i] = 0;
        } break;
      }
    }
  }
}

template<template <typename S> class VectorValuedFunctionTemplate>
void
PressureElement<VectorValuedFunctionTemplate>::neumVector(CoordSet& c0, Vector& F, int, GeomState *c1, double t)
{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, c1);

  // evaluate the function and store values terms
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  Eigen::Matrix<double,M,1> fval = f(q,t);

  // fill F
  for(int i=0; i<M; ++i) {
    F[i] = fval[i];
  }
}

template<template <typename S> class VectorValuedFunctionTemplate>
void 
PressureElement<VectorValuedFunctionTemplate>::neumVectorJacobian(CoordSet& c0, FullSquareMatrix& Ktan, int, GeomState* c1, double t)
{
  // instantiate the function object
  Eigen::Array<double, VectorValuedFunctionTemplate<double>::NumberOfScalarConstants, 1> sconst;
  Eigen::Array<int, VectorValuedFunctionTemplate<double>::NumberOfIntegerConstants, 1> iconst;
  getConstants(c0, sconst, iconst);
  VectorValuedFunctionTemplate<double> f(sconst,iconst);

  // prepare the function inputs
  const int N = VectorValuedFunctionTemplate<double>::NumberOfGeneralizedCoordinates;
  Eigen::Matrix<double,N,1> q;
  getInputs(q, c0, c1);

  // instantiate the jacobian object
  Simo::Jacobian<double,VectorValuedFunctionTemplate> dfdq(sconst,iconst);

  // evaluate the jacobian
  const int M = VectorValuedFunctionTemplate<double>::NumberOfValues;
  Eigen::Map<Eigen::Matrix<double,M,N,Eigen::RowMajor> > J(Ktan.data());
  J = dfdq(q, t);
}

template<template <typename S> class VectorValuedFunctionTemplate>
FullSquareMatrix
PressureElement<VectorValuedFunctionTemplate>
::sommerMatrix(CoordSet &cs, double *d)
{
  FullSquareMatrix sommerM(numDofs(),d);
  sommerM.zero();

  return sommerM;
}

template<template <typename S> class VectorValuedFunctionTemplate>
int
PressureElement<VectorValuedFunctionTemplate>
::findAndSetEle(CoordSet& cs, Elemset &eset, Connectivity *nodeToElem, int *eleTouch, int *eleCount, int myNum, int it)
{
  // overriding SommerElement::findAndSetEle because the normal should not be reversed
  this->iEle = findEle(nodeToElem, eleTouch, eleCount, myNum, &eset, it);
  if(iEle == -1) {
    std::cerr << "PressureElement::findAndSetEle could not find the corresponding element.\n";
    return 0;
  }
  el = eset[iEle];
  return -1;
}

#endif
