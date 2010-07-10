#include <Element.d/Rigid.d/RigidThreeNodeShell.h>
#include <Element.d/Rigid.d/RigidBeam.h>

RigidThreeNodeShell::RigidThreeNodeShell(int *_nn)
 : SuperElement(true)
{
  nnodes = 3;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  for(int i = 0; i < nSubElems; ++i) {
    int indices[2] = { i+1, 0 };
    subElems[i] = new RigidBeam(indices);
  }
}
