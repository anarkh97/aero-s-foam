#include <Element.d/Solid.d/RigidMpcSolid6Dof.h>
#include <Element.d/Beam.d/RigidMpcBeam.h>

// Rigid solid superelement comprising a number of rigid beams connected together

RigidMpcSolid6Dof::RigidMpcSolid6Dof(int numnodes, int *nodenums)
{
  nnodes = numnodes;
  ndofs = numnodes*6;

  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = nodenums[i];

  if(nnodes == 1) {
    cerr << " *** ERROR: RigidMpcSolid6Dof should have 2 or more nodes. Exiting...\n";
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
      subElems[i] = new RigidMpcBeam(subElemNodes[i]);
      subElems[i]->setGlNum(-1);
    }
  }
}

int
RigidMpcSolid6Dof::getTopNumber()
{
  switch(nnodes) {
    case 4 : return 123; // 4-node tetra
    case 6 : return 124; // 6-node penta
    case 8 : return 117; // 8-node hexa
    case 10 : return 125; // 10-node tetra
    case 15 : return 197; // 15-node hexa
    case 20 : return 172; // 20-node brick
    case 26 : return 192; // 26-node hexa 
    case 32 : return 191; // 32-node wedge
    default : return 101;
  }
}

int
RigidMpcSolid6Dof::numTopNodes()
{
  switch(nnodes) {
    case 4 : return 4; // 4-node tetra
    case 6 : return 6; // 6-node penta
    case 8 : return 8; // 8-node hexa
    case 10 : return 10; // 10-node tetra
    case 15 : return 15; // 15-node hexa
    case 20 : return 20; // 20-node brick
    case 26 : return 26; // 26-node hexa 
    case 32 : return 32; // 32-node wedge
    default : return 2;
  }
}

