#ifndef __PiezoQuad4.h__
#define __PiezoQuad4.h__
#ifdef STRUCTOPT

#include <Structopt.d/Element_opt.d/Quad4.d/FourNodeQuad_opt.h>
#include <Utils.d/dofset.h>

/* Four node piezo-patch; Temp-dof represents scalar potential in
 * our case. Has no structural properties, but only piezo & piezo-structural
 * couplings.
 */
class PiezoQuad4: public FourNodeQuad_opt
{
public:
  PiezoQuad4(int* n) : FourNodeQuad_opt(n), cCoefs(0), cFrame(0), d_cCoefs(0) /*, d_cFrame(0)*/ {}
  Element* clone() { return new PiezoQuad4(*this); }

  FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
  void gradstiffness(CoordSet& cs,CoordSet& dcs,FullSquareMatrix& dd, int flg = 1);
  int chkOptInf(CoordSet& dcs);
	
  double getMass(CoordSet&) { return 0.0; }  
  void   markDofs(DofSetArray& dsa) {dsa.mark(nn, numNodes(), DofSet::Xdisp|DofSet::Ydisp|DofSet::Temp);}
  int*   dofs(DofSetArray&, int *p=0);
  int    numDofs() { return numNodes()*dim(); }
  int    dim()     { return 3; }

  void gradMassMatrix(CoordSet&cs, CoordSet&dcs, FullSquareMatrix&m, double mratio=1) 
    { return Element_opt::gradMassMatrix(cs, dcs, m, mratio); }
  double getGradMass(CoordSet& cs, CoordSet& dcs) { return 0; }

  void setCompositeData(int _type, int nlays, double *lData,
			double *coefs, double *frame) { cCoefs = coefs; cFrame = frame; }

protected:
  FullSquareMatrix massMatrix(CoordSet& cs, double *d, int cmf=1) { return Element::massMatrix(cs, d, cmf); }
  void fillPiezoTensors(double e[3][2], double kappa[2][2]);
  void fillDPiezoTensors(double e[3][2], double kappa[2][2]);
  void piezo_quad4v(const double* x, const double* y, const double* z,
		    int p, 
		    const double kappa[2][2],
		    const double e[3][2],
		    double* d);

  void piezo_dquad4v(const double* x,  const double*  y, const double*  z,
		     const double* dx, const double* dy, const double* dz,
		     int p, 
		     const double kappa[2][2], const double dkappa[2][2],
		     const double e[3][2],     const double de[3][2],
		     double* d);

  //virtual const double* getDCFrame() const { return d_cFrame; }
  //virtual void setDCFrame(double* d) { d_cFrame = d; }
  virtual const double* getDCCoefs() const { return d_cCoefs; }
  virtual void setDCCoefs(double* d) { d_cCoefs = d; }

private:
  // data:
  double* cCoefs;
  double* cFrame;
  double* d_cCoefs;
  //double* d_cFrame;
};

#endif
#endif
