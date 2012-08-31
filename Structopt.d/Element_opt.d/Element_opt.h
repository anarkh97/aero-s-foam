#ifndef	_ELEMENT_OPT_HPP_
#define	_ELEMENT_OPT_HPP_

#ifdef STRUCTOPT

#include <Element.d/Element.h>

class Element_opt : virtual public Element
{
protected:  
  StructProp *gradprop;	// gradient of structural properties 
  
public:
  StructProp * getGradProp() { return gradprop; }

  virtual double getGradDensity() { if (gradprop) return gradprop->rho;
    else return 0.0; }
  
  void setGradProp(StructProp *gradp) { gradprop = gradp; }

  virtual int chkOptInf(CoordSet &) { return 0; }

  virtual int  chkOptInfTH() { return 0; }
  virtual int  chkOptInfEL() { return 0; }

  virtual double getGradMass(CoordSet&, CoordSet&)
  { fprintf(stderr,"WARNING: analytical deriv. of element mass not implemented\n");
    return 0.0; }
 	
  virtual double getGradDCmass(CoordSet &,CoordSet &,Vector &,Vector &, double&)
  { fprintf(stderr,"WARNING: analytical deriv. of element mass in def. config. not implemented\n");
    return 0.0; }

  virtual void computeDfDsPressure(CoordSet& cs, CoordSet& dcs, Vector& eleDfDs, GeomState *geomState=0, int cflg = 0);

  virtual void computeDfDsThermal(CoordSet&, CoordSet& dcs, Vector &,
				  Vector &, Vector &, int , 
				  GeomState *gs=0,int flag=0); 

  virtual void getGradVonMises(Vector &, Vector &, CoordSet &, CoordSet &, 
			       Vector &, Vector &, int, int, double* ndTemp=0, double*  dndTemp=0 );

  virtual void getGradVonMises(ComplexVector &, Vector &, CoordSet &, CoordSet &, 
			       ComplexVector &, ComplexVector &, int, int, double* ndTemp=0, double*  dndTemp=0 );

  virtual void getGradVonMisesInt(CoordSet &,CoordSet &,Vector &,Vector &,double &,
				  double &,int, double &,double &,double &,double &, 
				  double* ndTemp=0, double*  dndTemp=0);       

  virtual void gradstiffness(CoordSet &, CoordSet &, FullSquareMatrix &, int =1);
  virtual void gradMassMatrix(CoordSet &, CoordSet &, FullSquareMatrix &,  double mratio=1.0);
  virtual void gradDampMatrix(CoordSet &, CoordSet &, FullSquareMatrix &,  double freq=0.0);
  //virtual void gradstiffness(CoordSet &, CoordSet &, FullSquareMatrixC&, int =1);
  //virtual void gradMassMatrix(CoordSet &, CoordSet &, FullSquareMatrixC &, double mratio=1.0);
  //virtual void gradDampMatrix(CoordSet &, CoordSet &, FullSquareMatrixC &, double freq=0.0);

  virtual void setCompositeGrad( double * ) {}

  virtual void setdFrame(EFrame *) {} 
  virtual double * getdFrame()     { return 0; } 

  virtual void zerodFrame() {} 
  virtual void getGradElectricField(CoordSet&,CoordSet&,Vector&,Vector&,
				    Vector&,Vector&,Vector&,Vector&,
				    Vector&);

  virtual int  doesElemComputeEfield() { return 0; }
  virtual int  doesElemComputeEforce() { return 0; }
  virtual int  doesElemComputeEstiff() { return 0; }
  virtual int  doesElemComputeEflux()  { return 0; }

  virtual void getGradElectricForce(CoordSet&,CoordSet&,Vector&,Vector&,
				    Vector&,Vector&,Vector&,Vector&,
				    Vector&,Vector&,Vector&);

  virtual	void computeGradDisp(double[2],double*)  {}
        
  virtual void computeGradSourceFlux(Vector&,CoordSet&,Vector&,
				     Vector&, Vector&); 
  virtual void computeDConvectiveLoadDs(Vector& load,CoordSet&,CoordSet&,
					Vector&,Vector&,int flag =1)
  { load.zero(); }
  virtual void computeDRadiationLoadDs(Vector &df , CoordSet& , CoordSet& , 
				       Vector& , Vector&){ df.zero();}

  virtual const double* getDCCoefs() const { return 0; }
  virtual void setDCCoefs(double*) { assert(0); }
  //virtual const double* getDCFrame() const { return 0; }
  //virtual void setDCFrame(double*) {}

public:
  static void rotateConstitutiveMatrix33(double Cin[3][3], double *T33, double Cout[3][3]);
  static void rotateConstitutiveMatrix63(double Cin[6][3], double *T33, double Cout[6][3]);

};

class ElementFactory_opt : public ElementFactory
{
public:
  ElementFactory_opt() : ElementFactory() {}
  virtual Element* elemadd(int num, int type, int nnodes, int *nodes,
			   BlockAlloc& ba);
};


#endif
#endif

