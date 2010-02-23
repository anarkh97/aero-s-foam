#include <Corotational.d/SuperCorotator.h>

SuperCorotator::SuperCorotator(SuperElement *_superElem)
{
  superElem = _superElem;
  nSubElems = superElem->getNumSubElems();
  subElemCorotators = new Corotator * [nSubElems];
  origK = 0;
  sub_vld = 0;
  sub_vlr = 0;
}

SuperCorotator::~SuperCorotator() 
{ 
  delete [] subElemCorotators; 
  if(origK) delete origK;
  if(sub_vld) {
    int i;
    for(i=0; i<nSubElems; ++i)
      if(sub_vld[i]) delete [] sub_vld[i];
    delete [] sub_vld;
  }
  if(sub_vlr) {
    int i;
    for(i=0; i<nSubElems; ++i)
      if(sub_vlr[i]) delete [] sub_vlr[i];
    delete [] sub_vlr;
  }
}


void 
SuperCorotator::getStiffAndForce(GeomState &geomState, CoordSet &cs,
                                 FullSquareMatrix &elK, double *f) 
{
  int i, j;
  elK.zero();
  for(i=0; i<elK.dim(); ++i) f[i] = 0.0;

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    FullSquareMatrix subK(ndofs);
    subK.zero();  // this is necessary because getStiffAndForce may not be implemented for all subelems
    double *subf = new double[ndofs];
    for(j=0; j<ndofs; ++j) subf[j] = 0.0;
    subElemCorotators[i]->getStiffAndForce(geomState, cs, subK, subf);
    int *subElemDofs = superElem->getSubElemDofs(i);
    elK.add(subK, subElemDofs);
    for(j=0; j<ndofs; ++j) f[subElemDofs[j]] += subf[j];
    delete [] subf;
  }
  //cerr << "here in SuperCorotator::getStiffAndForce\n";
  //cerr << "elK = \n"; elK.print();
  //cerr << "f = \n"; for(int i=0; i<elK.dim(); ++i) cerr << f[i] << " "; cerr << endl;
}

void 
SuperCorotator::formGeometricStiffness(GeomState &geomState, CoordSet &cs, 
                                       FullSquareMatrix &elK, double *f)
{
  int i, j;
  elK.zero();
  for(i=0; i<superElem->numDofs(); ++i) f[i] = 0.0;

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    FullSquareMatrix subK(ndofs);
    subK.zero();  // this is necessary because formGeometricStiffness may not be implemented for all subelems
    double *subf = new double[ndofs];
    for(j=0; j<ndofs; ++j) subf[j] = 0.0;
    subElemCorotators[i]->formGeometricStiffness(geomState, cs, subK, subf);
    int *subElemDofs = superElem->getSubElemDofs(i);
    elK.add(subK, subElemDofs);
    for(j=0; j<ndofs; ++j) f[subElemDofs[j]] += subf[j];
    delete [] subf;
  }
}

double *
SuperCorotator::getOriginalStiffness() 
{
  if(!origK) origK = new FullM(superElem->numDofs());
  origK->zero();
  int i;
  for(i=0; i<nSubElems; ++i) {
    double *subKdata = subElemCorotators[i]->getOriginalStiffness();
    if(subKdata) {
      int ndofs = superElem->getSubElemNumDofs(i);
      FullM subK(subKdata, ndofs, ndofs);
      origK->add(subK, superElem->getSubElemDofs(i));
    }
  }
  return origK->data();
}

void 
SuperCorotator::extractDeformations(GeomState &geomState, CoordSet &cs, double *vld, int &nlflag) 
{
  int i, j;
  if(!sub_vld) { 
    sub_vld = new double * [nSubElems];
    for(i=0; i<nSubElems; ++i) sub_vld[i] = 0;
  }
  for(i=0; i<superElem->numDofs(); ++i) vld[i] = 0.0; // this shouldn't be used, should always use sub_vld for non-linear
                                                      // because the same dof can have 2 different vld displacements in different sub elements
  nlflag = 999;

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    if(!sub_vld[i]) sub_vld[i] = new double[ndofs];
    for(j=0; j<ndofs; ++j) sub_vld[i][j] = 0.0;
    int subnlflag;
    subElemCorotators[i]->extractDeformations(geomState, cs, sub_vld[i], subnlflag);
    if(subnlflag < nlflag) nlflag = subnlflag;  // always take the lowest
  }
}

void 
SuperCorotator::getNLVonMises(Vector &stress, Vector &weight, GeomState &geomState,
                              CoordSet &cs, int strInd)
{
  int i;
  stress.zero();
  weight.zero();

  for(i=0; i<nSubElems; ++i) {
    int nnodes = superElem->getSubElemNumNodes(i);
    Vector subStress(nnodes);
    subStress.zero();
    Vector subWeight(nnodes);
    subWeight.zero();
    subElemCorotators[i]->getNLVonMises(subStress, subWeight, geomState, cs, strInd);
    int *subElemNodes = superElem->getSubElemNodes(i);
    stress.add(subStress, subElemNodes);
    weight.add(subWeight, subElemNodes);
  }
}

void 
SuperCorotator::getNLAllStress(FullM &stress, Vector &weight, GeomState &geomState, 
                               CoordSet &cs, int strInd) 
{
  int i;
  stress.zero();
  weight.zero();

  for(i=0; i<nSubElems; ++i) {
    int nnodes = superElem->getSubElemNumNodes(i);
    FullM subStress(nnodes, 9);
    subStress.zero();
    Vector subWeight(nnodes);
    subWeight.zero();
    subElemCorotators[i]->getNLAllStress(subStress, subWeight, geomState, cs, strInd);
    int *subElemNodes = superElem->getSubElemNodes(i);
    stress.addrows(subStress, subElemNodes);
    weight.add(subWeight, subElemNodes);
  }
}

double 
SuperCorotator::getElementEnergy(GeomState &geomState, CoordSet &cs) 
{
  int i;
  double ret = 0.0;
  for(i=0; i<nSubElems; ++i)
    ret += subElemCorotators[i]->getElementEnergy(geomState, cs);
  return ret;
}

void 
SuperCorotator::extractRigidBodyMotion(GeomState &geomState, CoordSet &cs, double *vlr) 
{
  int i,j;
  if(!sub_vlr) {
    sub_vlr = new double * [nSubElems];
    for(i=0; i<nSubElems; ++i) sub_vlr[i] = 0;
  }
  for(i=0; i<superElem->numDofs(); ++i) vlr[i] = 0.0; // this shouldn't be used, should always use sub_vlr for non-linear
                                                      // because the same dof can have 2 different vlr displacements in different sub elements
  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    if(!sub_vlr[i]) sub_vlr[i] = new double[ndofs];
    for(j=0; j<ndofs; ++j) sub_vlr[i][j] = 0.0;
    subElemCorotators[i]->extractRigidBodyMotion(geomState, cs, sub_vlr[i]);
  }
}

