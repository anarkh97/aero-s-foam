#include <Element.d/Rigid.d/RigidThreeNodeShell.h>
#include <Element.d/Rigid.d/RigidBeam.h>

// Rigid three node shell superelement comprising two rigid beams connected together

RigidThreeNodeShell::RigidThreeNodeShell(int *_nn)
{
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  for(int i = 0; i < nSubElems; ++i) {
    int indices[2] = { i, i+1 };
    subElems[i] = new RigidBeam(indices);
  }
  initialize(3, _nn);
}
