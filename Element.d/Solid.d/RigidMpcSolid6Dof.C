#include <Element.d/Solid.d/RigidMpcSolid6Dof.h>
#include <Element.d/Beam.d/RigidMpcBeam.h>

// Rigid solid superelement comprising a number of rigid beams connected together

RigidMpcSolid6Dof::RigidMpcSolid6Dof(int numnodes, int *nodenums)
{
  int i,j;
  nnodes = numnodes;
  ndofs = numnodes*6;

  nn = new int[nnodes];
  for(i=0; i<nnodes; ++i) nn[i] = nodenums[i];

  if(nnodes == 1) {
    cerr << " *** WARNING: RigidMpcSolid6Dof should have 2 or more nodes \n";
    nSubElems = 0;
    subElems = 0;
    subElemNodes = 0;
    subElemDofs = 0;
  }
  else {
    nSubElems = nnodes - 1;
    subElems = new Element * [nSubElems];
    subElemNodes = new int * [nSubElems];
    subElemDofs = new int * [nSubElems];

    for(i=0; i<nSubElems; ++i) {
      subElemNodes[i] = new int[2];
      subElemNodes[i][0] = i; subElemNodes[i][1] = i+1; 
      subElemDofs[i] = new int[12];
      for(j=0; j<6; ++j) {
        subElemDofs[i][j] = 6*i+j;
        subElemDofs[i][j+6] = 6*(i+1)+j;
      }
      int glSubElemNodes[2];
      glSubElemNodes[0] = nn[subElemNodes[i][0]];
      glSubElemNodes[1] = nn[subElemNodes[i][1]];
      subElems[i] = new RigidMpcBeam(glSubElemNodes);
      subElems[i]->setGlNum(-1);
    }
  }
}

Element *
RigidMpcSolid6Dof::clone()
{
  return new RigidMpcSolid6Dof(*this);
}

