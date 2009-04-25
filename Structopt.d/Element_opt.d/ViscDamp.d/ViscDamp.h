#ifndef _VISCDAMP_H_
#define _VISCDAMP_H_
#ifdef STRUCTOPT

#include <Element.d/Element.h>
#include <Structopt.d/Element_opt.d/Element_opt.h>


class ViscDamp: virtual public Element, public Element_opt
{
 private:
  const static double a;
  const static double b;

  FullSquareMatrix stiffness(CoordSet&, double* m, int cmflag=1)  
  { FullSquareMatrix ret(numDofs(),m); ret.zero(); return ret; }
  FullSquareMatrix massMatrix(CoordSet&, double* m, int cmflag=1) 
  { FullSquareMatrix ret(numDofs(),m); ret.zero(); return ret; }
  bool hasDamping() { return true; }

 protected:
  double V_p();
  double V_s();
  virtual double area(CoordSet& cs) = 0;
  double c00(CoordSet&);
  double c11(CoordSet&);

  double dV_p();
  double dV_s();
  virtual double darea(CoordSet&, CoordSet&) = 0;
  double dc00(CoordSet&, CoordSet&);
  double dc11(CoordSet&, CoordSet&);
};

#endif
#endif
