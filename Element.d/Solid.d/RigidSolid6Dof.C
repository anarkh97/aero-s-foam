#include <Element.d/Solid.d/RigidSolid6Dof.h>
#include <Element.d/Beam.d/RigidBeam.h>

// Rigid solid superelement comprising a number of rigid beams connected together

RigidSolid6Dof::RigidSolid6Dof(int numnodes, int *nodenums)
{
  nnodes = numnodes;
  ndofs = numnodes*6;

  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = nodenums[i];

  if(nnodes < 2) {
    cerr << " *** ERROR: RigidSolid6Dof should have 2 or more nodes. Exiting... \n";
    exit(-1);
  }
  else {
    nSubElems = nnodes - 1;
    subElems = new Element * [nSubElems];
    subElemNodes = new int * [nSubElems];
    subElemDofs = new int * [nSubElems];

    for(int i = 0; i < nSubElems; ++i) {
      subElemNodes[i] = new int[2];
      subElemNodes[i][0] = i; subElemNodes[i][1] = i+1; 
      subElemDofs[i] = new int[12];
      for(int j = 0; j < 6; ++j) {
        subElemDofs[i][j] = 6*i+j;
        subElemDofs[i][j+6] = 6*(i+1)+j;
      }
      subElems[i] = new RigidBeam(subElemNodes[i]);
      subElems[i]->setGlNum(-1);
    }
  }
}

Element *
RigidSolid6Dof::clone()
{
  return new RigidSolid6Dof(*this);
}

