#include <stdio.h>
#include <Utils.d/MyComplex.h>
#include <HelmAxi.d/MPCData.h>


MPCData *globalMPCs = 0;


AxiMPC zeroMPC = MPC(-1, 0.0, 0.0, 1, DComplex(1.0,0.0), 0.0, 0.0, 0.0);


MPCData::MPCData() : 
          term(zeroMPC)
{
  numMPCSet = 0;
  totalNumMPC = 0;
}


void MPCData::addMPC(AxiMPC &mpc) {
  term[numMPCSet++] = mpc;
}


