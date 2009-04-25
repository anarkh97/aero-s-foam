#ifndef _MODES_H_
#define _MODES_H_

// Structure to contain Mode Data
// i.e. number of modes, number of nodes, 
// array of frequencies, and matrix of mode shapes.

struct ModeData {
  int numModes;
  int numNodes;
  double *frequencies;
  double (**modes)[6];
  int *nodes;

  /*int multY(int numY, BCond *Y, BCond *&iDis, int ndof = 6) {
    if(iDis == 0) iDis = new BCond[numNodes*ndof];
    for(int j = 0; j < numNodes; ++j) 
      for(int k = 0; k < ndof; ++k) 
        iDis[ndof*j + k].setData(nodes[j], k, 0.0);

    for(int i = 0; i<numY; ++i) 
      for(int j = 0; j < numNodes; ++j) 
        for(int k = 0; k < ndof; ++k) 
          iDis[ndof*j + k].val += Y[i].val*modes[Y[i].nnum][j][k];

    return numNodes*ndof;
  }*/

  template<class VecType>
  void addMultY(int numY, BCond *Y, VecType& x, DofSetArray* dsa,  int ndof=6) { //HB: x <- x + X.y 
    // The following implementation assumes that all the modes have been given for the SAME nodes
    // to reduce the number of calls to DofSetArray::locate
    for(int j = 0; j < numNodes; ++j) // loop over nodes
      for(int k = 0; k < ndof; ++k) { // loop over dofs
         int dof = dsa->locate(nodes[j], 1 << k);
         if(dof >= 0)
           for(int i = 0; i<numY; ++i) // loop over selected modes
             x[dof] += Y[i].val*modes[Y[i].nnum][j][k];
      }
  }  
};

#endif
