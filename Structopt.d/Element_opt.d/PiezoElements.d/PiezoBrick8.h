#ifndef __PiezoBrick8.h__
#define __PiezoBrick8.h__
#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/Brick.d/EightNodeBrick_opt.h>
#include <Utils.d/dofset.h>

/* Four node piezo-patch; Temp-dof represents scalar potential in
 * our case. Has no structural properties, but only piezo & piezo-structural
 * couplings.
 */
class PiezoBrick8: public EightNodeBrick_opt
{
public:
  PiezoBrick8(int* n) : EightNodeBrick_opt(n), d_cCoefs(0) /*, d_cFrame(0)*/ {}
  Element* clone() { return new PiezoBrick8(*this); }

  FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
  void gradstiffness(CoordSet& cs,CoordSet& dcs,FullSquareMatrix& dd, int flg = 1);
  int chkOptInf(CoordSet& dcs);
	
  double getMass(CoordSet&) { return 0.0; }  
  void   markDofs(DofSetArray& dsa) {dsa.mark(nn, numNodes(), DofSet::XYZdisp|DofSet::Temp);}
  int*   dofs(DofSetArray&, int *p=0);
  int    numDofs() { return numNodes()*dim(); }
  int    dim()     { return 4; }

  void gradMassMatrix(CoordSet&cs, CoordSet&dcs, FullSquareMatrix&m, double mratio=1) 
    { return Element_opt::gradMassMatrix(cs, dcs, m, mratio); }
  double getGradMass(CoordSet& cs, CoordSet& dcs) { return 0; }

protected:
  FullSquareMatrix massMatrix(CoordSet& cs, double *d, int cmf=1) { return Element::massMatrix(cs, d, cmf); }
  void fillPiezoTensors(double e[6][3], double kappa[3][3]);
  void fillDPiezoTensors(double e[6][3], double kappa[3][3]);
  void piezo_br8v(const double* x, const double* y, const double* z,
		  int p, 
		  const double kappa[3][3],
		  const double e[6][3],
		  double* d);

  void piezo_dbr8v(const double* x,  const double*  y, const double*  z,
		   const double* dx, const double* dy, const double* dz,
		   int p, 
		   const double kappa[3][3], const double dkappa[3][3],
		   const double e[6][3],     const double de[6][3],
		   double* d);

  //virtual const double* getDCFrame() const { return d_cFrame; }
  //virtual void setDCFrame(double* d) { d_cFrame = d; }
  virtual const double* getDCCoefs() const { return d_cCoefs; }
  virtual void setDCCoefs(double* d) { d_cCoefs = d; }

private:
  // data:
  double* d_cCoefs;
  //double* d_cFrame;
};

#endif
#endif
