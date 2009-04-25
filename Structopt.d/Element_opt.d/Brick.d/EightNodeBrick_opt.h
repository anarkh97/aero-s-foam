#ifndef _EIGHTNODEBRICK_OPT_HPP_
#define _EIGHTNODEBRICK_OPT_HPP_
#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/Element_opt.h>
#include <Element.d/Brick.d/EightNodeBrick.h>

class EightNodeBrick_opt: public EightNodeBrick, public Element_opt
{
public:
  EightNodeBrick_opt(int* nn) : EightNodeBrick(nn) {}
  int chkOptInf(CoordSet&);
  double* getMidPoint(CoordSet &cs);
  void getMidPoint(CoordSet &cs, double* result);
  void gradstiffness(CoordSet &, CoordSet &, FullSquareMatrix &, int =1); 
  void gradMassMatrix(CoordSet &, CoordSet &, FullSquareMatrix &, double mratio=1);
        
  void getGradVonMises (Vector &, Vector &, CoordSet &, CoordSet &, 
			Vector &, Vector &, int , int =0, 
			double* ndTemp=0, double*  dndTemp=0 );
  void getGradVonMisesInt(CoordSet &, CoordSet &, Vector &, Vector &,
			  double &,  double &, int, double &, double &, 
			  double &, double &, 
			  double* ndTemp=0, double*  dndTemp=0);
  double getGradMass(CoordSet &, CoordSet &);
  int getTopNumber() { return 3; }
};

#endif
#endif

