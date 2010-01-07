#include <Element.d/Joint.d/TranslationalJoint.h>
#include <Element.d/Joint.d/ParallelAxisConstraintType1.h>
#include <Element.d/Joint.d/DotConstraintType1.h>
#include <Element.d/Joint.d/ParallelAxisConstraintType2.h>

TranslationalJoint::TranslationalJoint(int* _nn)
{
  nSubElems = 3, nnodes = 2; int nDofsPerNode = 6;
  subElems = new Element * [nSubElems];
  subElems[0] = new ParallelAxisConstraintType1(_nn);
  subElems[1] = new DotConstraintType1(_nn, 2, 1);
  subElems[2] = new ParallelAxisConstraintType2(_nn);
  for(int i = 0; i < nSubElems; ++i) subElems[i]->setGlNum(-1);

  subElemNodes = new int * [nSubElems];
  int l = 0;
  for(int i = 0; i < nSubElems; ++i) {
    subElemNodes[i] = new int[subElems[i]->numNodes()];
    for(int j = 0; j < nnodes; ++j) subElemNodes[i][j] = j;
    for(int j = nnodes; j < nnodes + subElems[i]->numInternalNodes(); ++j) subElemNodes[i][j] = nnodes + l++;
  }
  nn = new int[nnodes+l];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  
  subElemDofs = new int * [nSubElems];
  l = 0;
  for(int i = 0; i < nSubElems; ++i) {
    subElemDofs[i] = new int[subElems[i]->numDofs()];
    int m = (subElems[i]->numDofs()-subElems[i]->numInternalNodes())/nnodes;
    for(int j = 0; j < nnodes; ++j) for(int k = 0; k < m; ++k) subElemDofs[i][j*m+k] = j*nDofsPerNode+k;
    for(int j = nnodes*m; j < subElems[i]->numDofs(); ++j) subElemDofs[i][j] = nnodes*nDofsPerNode + l++;
  }

  ndofs = nnodes*nDofsPerNode + l;
  nnodes = nnodes + l; 
}

int 
TranslationalJoint::getTopNumber() 
{ 
  return 106; 
}
