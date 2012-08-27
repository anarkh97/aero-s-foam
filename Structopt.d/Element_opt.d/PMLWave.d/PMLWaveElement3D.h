#ifndef _PMLWAVE3D_H_
#define _PMLWAVE3D_H_
#ifdef STRUCTOPT 

#include <Structopt.d/Element_opt.d/Brick.d/EightNodeBrick_opt.h>

#define PMLWAVE3D_NNODES 8
#define PMLWAVE3D_NDIM   3
class PMLWaveElement3D: public EightNodeBrick_opt
{
 private:
  const static int nDim;

 public:
  PMLWaveElement3D(int* n) : EightNodeBrick_opt(n) {}
  Element* clone();

  bool isComplex() { return true; }
  FullSquareMatrix  stiffness (CoordSet& cs, double* s, int flg=1) { return Element::stiffness (cs, s, flg); }
  FullSquareMatrix  massMatrix(CoordSet& cs, double* m, int flg=1) { return Element::massMatrix(cs, m, flg); }
  FullSquareMatrixC complexStiffness(CoordSet&, DComplex*, int flg=1);
  FullSquareMatrixC complexMassMatrix(CoordSet&, DComplex*, double mratio);

  int chkOptInf(CoordSet&) { return 0; }
};

#endif
#endif
