#include <Element.d/Solid.d/RigidMpcSolid.h>
#include <Element.d/Truss.d/TwoNodeTrussRigidMpc.h>

// Rigid solid superelement comprising a number of rigid trusses connected together

RigidMpcSolid::RigidMpcSolid(int numnodes, int *nodenums)
{
  int i,j;
  nnodes = numnodes;
  ndofs = numnodes*3;

  nn = new int[nnodes];
  for(i=0; i<nnodes; ++i) nn[i] = nodenums[i];

  if(nnodes == 1) {
    cerr << " *** WARNING: RigidMpcSolid should have 2 or more nodes \n";
    nSubElems = 0;
    subElems = 0;
    subElemNodes = 0;
    subElemDofs = 0; 
  }
  else {
    if(nnodes == 2) nSubElems = 1;
    else nSubElems = nnodes;

    subElems = new Element * [nSubElems];
    subElemNodes = new int * [nSubElems];
    subElemDofs = new int * [nSubElems];

    for(i=0; i<(nnodes-1); ++i) {
      subElemNodes[i] = new int[2];
      subElemNodes[i][0] = i; subElemNodes[i][1] = i+1; 
      subElemDofs[i] = new int[6];
      for(j=0; j<3; ++j) {
        subElemDofs[i][j] = 3*i+j;
        subElemDofs[i][j+3] = 3*(i+1)+j;
      }
      int glSubElemNodes[2];
      glSubElemNodes[0] = nn[subElemNodes[i][0]];
      glSubElemNodes[1] = nn[subElemNodes[i][1]];
      subElems[i] = new TwoNodeTrussRigidMpc(glSubElemNodes);
      subElems[i]->setGlNum(-1);
    }
    if(nnodes > 2) { // close the loop
      i = nnodes-1;
      subElemNodes[i] = new int[2];
      subElemNodes[i][0] = i; subElemNodes[i][1] = 0;
      subElemDofs[i] = new int[6];
      for(j=0; j<3; ++j) {
        subElemDofs[i][j] = 3*i+j;
        subElemDofs[i][j+3] = j;
      }
      int glSubElemNodes[2];
      glSubElemNodes[0] = nn[subElemNodes[i][0]];
      glSubElemNodes[1] = nn[subElemNodes[i][1]];
      subElems[i] = new TwoNodeTrussRigidMpc(glSubElemNodes);
      subElems[i]->setGlNum(-1);
    }
  }
}

Element *
RigidMpcSolid::clone()
{
  return new RigidMpcSolid(*this);
}

