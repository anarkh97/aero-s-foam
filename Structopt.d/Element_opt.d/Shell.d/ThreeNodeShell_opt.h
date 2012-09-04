#ifndef _THREENODESHELL_OPT_H_
#define _THREENODESHELL_OPT_H_

#ifdef STRUCTOPT

#include <Element.d/Shell.d/ThreeNodeShell.h>
#include <Structopt.d/Element_opt.d/Element_opt.h>

class ThreeNodeShell_opt : public Element_opt, public ThreeNodeShell 
{
private:
  double alpha, beta;
public:
  ThreeNodeShell_opt(int* nn, double w=3.0, double l_alpha = 1.5, double l_beta = 0.32) : 
    ThreeNodeShell(nn, w), alpha(l_alpha), beta(l_beta) {}
  int chkOptInf(CoordSet&);  

  double* getMidPoint(CoordSet &cs);
  void getMidPoint(CoordSet &cs, double* midPoint);

  double getGradMass(CoordSet& cs, CoordSet& dcs);  
  double getGradDCmass(CoordSet &,CoordSet &,Vector &,Vector &,double&);  

  void computeDfDsPressure(CoordSet&, CoordSet&, Vector&, 
			   GeomState *gs=0, int cflg = 0);  

  void getGradVonMises(Vector &dstress, Vector &weight, 
		       CoordSet &cs, CoordSet &dcs, 
		       Vector &elDisp, Vector &elGrad,
		       int strInd, int surface, 
		       double* ndTemp=0, double* dndTemp=0);

  void getGradVonMisesInt(CoordSet &,CoordSet &, Vector &, Vector &, 
			  double &, double &, int, double &, double &, 
			  double &, double &, 
			  double* ndTemp=0, double* dndTemp=0);

  void gradstiffness(CoordSet& cs, CoordSet& gradcs, 
		     FullSquareMatrix & gradkarray, int flg=1);  

  void gradMassMatrix(CoordSet &cs,CoordSet &dcs,FullSquareMatrix &dmel, double mratio=1);  

  void computeGradDisp(double[2],double*);  

  int getTopNumber() { return 4; }
};
#endif
#endif
