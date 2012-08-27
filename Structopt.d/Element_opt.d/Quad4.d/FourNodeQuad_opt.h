#ifndef _FOURNODEQUAD_OPT_H_
#define _FOURNODEQUAD_OPT_H_

#ifdef STRUCTOPT

#include <Element.d/Quad4.d/FourNodeQuad.h>
#include <Structopt.d/Element_opt.d/Element_opt.h>

class FourNodeQuad_opt : public Element_opt, public FourNodeQuad
{
private:
  double alpha, beta;
public:
  FourNodeQuad_opt(int* nn) : FourNodeQuad(nn) {}
  int chkOptInf(CoordSet&);  

  double* getMidPoint(CoordSet &cs);
  void getMidPoint(CoordSet &cs, double* result);

  double getGradMass(CoordSet& cs, CoordSet& dcs);  
  void getGradVonMises(Vector &dstress, Vector &weight, 
		       CoordSet &cs, CoordSet &dcs, 
		       Vector &elDisp, Vector &elGrad,
		       int strInd, int surface, 
		       double* ndTemp=0, double* dndTemp=0);

  void gradstiffness(CoordSet& cs, CoordSet& gradcs, 
		     FullSquareMatrix & gradkarray, int flg=1);  

  void gradMassMatrix(CoordSet &cs,CoordSet &dcs,FullSquareMatrix &dmel, double mratio=1);  

  void computeGradDisp(double[2],double*);  

  int getTopNumber() { return 2; }
};
#endif
#endif
