#include <Element.d/Brick.d/EightNodeBrickRigid.h>
#include <Element.d/Truss.d/TwoNodeTrussRigid.h>

// Rigid eight node brick superelement comprising 18 rigid two-node truss elements connected together

EightNodeBrickRigid::EightNodeBrickRigid(int *nodenums)
{
  nnodes = 8;
  ndofs = 24;

  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = nodenums[i];

  nSubElems = 18;

  subElems = new Element * [nSubElems];
  subElemNodes = new int * [nSubElems];
  subElemDofs = new int * [nSubElems];

  int indices[18][2] = {{0,1},
                        {1,2},
                        {2,3},
                        {4,5},
                        {5,6},
                        {6,7},
                        {0,2},
                        {4,6},
                        {0,3},
                        {4,7},
                        {0,4},
                        {1,5},
                        {2,6},
                        {3,7},
                        {0,5},
                        {1,6},
                        {2,7},
                        {0,7}};

  for(int i = 0; i < nSubElems; ++i) {
    subElemNodes[i] = new int[2];
    subElemNodes[i][0] = indices[i][0]; subElemNodes[i][1] = indices[i][1]; 
    subElemDofs[i] = new int[6];
    for(int j = 0; j < 3; ++j) {
      subElemDofs[i][j] = 3*subElemNodes[i][0]+j;
      subElemDofs[i][j+3] = 3*subElemNodes[i][1]+j;
    }
    subElems[i] = new TwoNodeTrussRigid(subElemNodes[i]);
    subElems[i]->setGlNum(-1);
  }
}

Element *
EightNodeBrickRigid::clone()
{
  return new EightNodeBrickRigid(*this);
}

