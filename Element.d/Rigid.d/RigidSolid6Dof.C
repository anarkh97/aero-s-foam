#include <Element.d/Rigid.d/RigidSolid6Dof.h>
#include <Element.d/Rigid.d/RigidBeam.h>

// Rigid solid superelement comprising a number of rigid beams connected together

RigidSolid6Dof::RigidSolid6Dof(int _nnodes, int *_nn)
{
  nSubElems = _nnodes - 1;
  subElems = new Element * [nSubElems];
  for(int i = 0; i < nSubElems; ++i) {
    int indices[2] = { i+1, 0 }; // all sub-elements have the same master node
    subElems[i] = new RigidBeam(indices);
  }
  initialize(_nnodes, _nn);
}

int
RigidSolid6Dof::getTopNumber()
{
  switch(nnodes-numInternalNodes()) {
    case 4 : return 123; // 4-node tetra
    case 6 : return 124; // 6-node penta
    case 8 : return 117; // 8-node hexa
    case 10 : return 125; // 10-node tetra
    case 15 : return 197; // 15-node penta
    case 20 : return 172; // 20-node hexa
    case 26 : return 192; // 26-node penta
    case 32 : return 191; // 32-node hexa
    default : return 101;
  }
}

int
RigidSolid6Dof::numTopNodes()
{
  switch(nnodes-numInternalNodes()) {
    case 4 : return 4; // 4-node tetra
    case 6 : return 6; // 6-node penta
    case 8 : return 8; // 8-node hexa
    case 10 : return 10; // 10-node tetra
    case 15 : return 15; // 15-node penta
    case 20 : return 20; // 20-node hexa
    case 26 : return 26; // 26-node penta
    case 32 : return 32; // 32-node hexa
    default : return 2;
  }
}

