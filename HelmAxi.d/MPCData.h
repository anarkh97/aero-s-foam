#ifndef _MPCDATA_H_
#define _MPCDATA_H_


#include <Utils.d/resize_array.h>


struct AxiMPC {
  int nnum; 
  double angle1;
  double angle2;
  int slices;
// Complex value of type c*e^(i*k*(x*dx+y*dy+z*dz))
  double cr, ci;
  double dx,dy,dz;
};


inline AxiMPC
MPC(int node, double angle1, double angle2, int slices, DComplex _c, 
    double _dx, double _dy, double _dz) {
  AxiMPC r;
  r.nnum = node;
  r.angle1 = angle1;
  r.angle2 = angle2;
  r.slices = slices;
  r.cr = _c.real(); 
  r.ci = _c.imag();
  r.dx = _dx; 
  r.dy = _dy; 
  r.dz = _dz;
  return r;
}


class MPCData {

public : 

  int numMPCSet;
  int totalNumMPC;
  ResizeArray<AxiMPC> term;   

  MPCData();
  void addMPC(AxiMPC &mpc);
 
};


extern MPCData *globalMPCs;


#endif
