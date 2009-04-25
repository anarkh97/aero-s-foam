#include <Element.d/Rigid.d/RBE2Mpc.h>
#include <Element.d/Rigid.d/GenRigidMpcBeam.h>

// General rigid solid superelement comprising a number of general rigid beams connected together
// same as Nastran RBE2 element
// implemented as linear multiple point equations & used by FETI solvers
// 1st term in "nodenums" array is the master node
// 2nd term in "nodenums" array is the dof list eg 123456 for rigid all trans + all rot 
//                                                 123 for rigid all trans
// remaining terms in "nodenums" array are the slave nodes

RBE2Mpc::RBE2Mpc(int numnodes, int *nodenums)
{
  int i,j;
  nnodes = numnodes-1;
  sprintf(cdofs,"%d\0",nodenums[1]+1);
  numcdofs = 0;
  for(int i=0; i<8; ++i) {
    if(cdofs[i] != '\0') numcdofs++;
    else break;
  }
  // cerr << "RBE2Mpc constructor: master node = " << nodenums[0]
  //      << ", numcdofs = " << numcdofs << ", cdofs = " << cdofs << endl;;

  ndofs = (nnodes-1)*2*numcdofs;

  nn = new int[nnodes];
  nn[0] = nodenums[0];
  for(i=1; i<nnodes; ++i) nn[i] = nodenums[i+1];

  if(nnodes == 1) {
    cerr << " *** WARNING: RBE2Mpc should have 2 or more nodes \n";
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
      subElemNodes[i][0] = 0; // this is deliberate, all sub-elements have same master node
      subElemNodes[i][1] = i+1; 
      subElemDofs[i] = new int[2*numcdofs];
      for(j=0; j<numcdofs; ++j) {
        subElemDofs[i][j] = numcdofs*subElemNodes[i][0]+j;
        subElemDofs[i][j+numcdofs] = numcdofs*subElemNodes[i][1]+j;
      }
      int glSubElemNodes[2];
      glSubElemNodes[0] = nn[subElemNodes[i][0]];
      glSubElemNodes[1] = nn[subElemNodes[i][1]];
      subElems[i] = new GenRigidMpcBeam(glSubElemNodes, numcdofs, cdofs);
      subElems[i]->setGlNum(-1);
    }
  }
}

Element *
RBE2Mpc::clone()
{
  return new RBE2Mpc(*this);
}

