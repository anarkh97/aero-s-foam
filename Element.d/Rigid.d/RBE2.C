#include <Element.d/Rigid.d/RBE2.h>
#include <Element.d/Rigid.d/GenRigidBeam.h>

// General rigid solid superelement comprising a number of general rigid beams connected together
// same as Nastran RBE2 element
// implemented by lagrange multipliers, used by direct solver
// 1st term in "nodenums" array is the master node
// 2nd term in "nodenums" array is the dof list eg 123456 for rigid all trans + all rot 
//                                                 123 for rigid all trans
// remaining terms in "nodenums" array are the slave nodes

RBE2::RBE2(int numnodes, int *nodenums)
{
  int i,j;
  nnodes = numnodes-1;
  sprintf(cdofs,"%d\0",nodenums[1]+1);
  numcdofs = 0;
  for(int i=0; i<8; ++i) {
    if(cdofs[i] != '\0') numcdofs++;
    else break;
  }
  // cerr << "RBE2 constructor: master node = " << nodenums[0] 
  //      << ", numcdofs = " << numcdofs << ", cdofs = " << cdofs << endl;;

  ndofs = 6*nnodes + numcdofs*(nnodes-1);

  nn = new int[2*nnodes-1];
  nn[0] = nodenums[0];
  for(i=1; i<nnodes; ++i) nn[i] = nodenums[i+1];

  // cerr << "nn = "; for(int i=0; i<nnodes; ++i) cerr << nn[i] << " "; cerr << endl;

  if(nnodes == 1) {
    cerr << " *** WARNING: RBE2 should have 2 or more nodes \n";
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
      subElemNodes[i] = new int[3];
      subElemNodes[i][0] = 0; // this is deliberate, all sub-elements have same master node
      subElemNodes[i][1] = i+1; 
      subElemNodes[i][2] = nnodes+i;
      subElemDofs[i] = new int[12 + numcdofs];
      for(j=0; j<6; ++j) {
        subElemDofs[i][j] = 6*subElemNodes[i][0]+j;
        subElemDofs[i][j+6] = 6*subElemNodes[i][1]+j;
      }
      for(j=0; j<numcdofs; ++j) {
        subElemDofs[i][12+j] = 6*nnodes+numcdofs*i+j;
      }
      int glSubElemNodes[2];
      glSubElemNodes[0] = nn[subElemNodes[i][0]];
      glSubElemNodes[1] = nn[subElemNodes[i][1]];
      subElems[i] = new GenRigidBeam(glSubElemNodes, numcdofs, cdofs);
      subElems[i]->setGlNum(-1);
    }
  }
  nnodes += nSubElems; // these are the internal nodes
}

Element *
RBE2::clone()
{
  return new RBE2(*this);
}

