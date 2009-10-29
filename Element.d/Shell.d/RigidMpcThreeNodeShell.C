#include <Element.d/Shell.d/RigidMpcThreeNodeShell.h>
#include <Element.d/Beam.d/RigidMpcBeam.h>

// Rigid three node shell superelement comprising two rigid beams connected together

RigidMpcThreeNodeShell::RigidMpcThreeNodeShell(int *nodenums)
{
  nnodes = 3;
  ndofs = nnodes*6;

  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = nodenums[i];

  nSubElems = 2;
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
    subElems[i] = new RigidMpcBeam(subElemNodes[i]);
    subElems[i]->setGlNum(-1);
  }
}

int
RigidMpcThreeNodeShell::getTopNumber()
{
  return 108;
}

