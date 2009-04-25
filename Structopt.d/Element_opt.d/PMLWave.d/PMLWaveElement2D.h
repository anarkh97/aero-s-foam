#ifndef _PMLWAVE2D_H_
#define _PMLWAVE2D_H_
#ifdef STRUCTOPT 

#include <Structopt.d/Element_opt.d/Quad4.d/FourNodeQuad_opt.h>

class PMLWaveElement2D: public FourNodeQuad_opt
{
 public:
  PMLWaveElement2D(int* n) : FourNodeQuad_opt(n) {}
  Element* clone();
    
  bool isComplex() { return true; }
  FullSquareMatrixC complexStiffness(CoordSet&, DComplex*, int flg=1);
  FullSquareMatrixC complexMassMatrix(CoordSet&, DComplex*, double mratio);

  int chkOptInf(CoordSet&) { return 0; }
};

#endif
#endif
