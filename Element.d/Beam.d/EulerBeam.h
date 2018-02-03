#ifndef _EULERBEAM_H_
#define _EULERBEAM_H_

#include <Element.d/Element.h>
#include <Element.d/Beam.d/BeamElementTemplate.hpp>

class GeomState;

class EulerBeam : public Element,
                  public BeamElementTemplate<double>
{
	EFrame *elemframe;
	int nn[3];
	double *offset;
	double c0[3][3];
	PressureBCond *pbc;
public:

	explicit EulerBeam(int*);

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	void setFrame(EFrame *ef) override { elemframe = ef; }
	const EFrame *getFrame() const override { return elemframe; }
	void buildFrame(CoordSet&) override;

	FullSquareMatrix stiffness(CoordSet&, double *kel,int flg) override;
	FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg) override;
	double getMass(CoordSet& cs) override;
	void   getGravityForce(CoordSet&, double *g, Vector &f, int gravflg,
	                       GeomState *gs) override;
	void   getIntrnForce(Vector& elForce, CoordSet& cs,
	                     double* elDisp,int forceIndex, double *ndTemps) override;

	void   markDofs(DofSetArray &) override;
	int*   dofs(DofSetArray &, int *p) override;
	int    numDofs() const override;

	int    numNodes() const override;
	int*   nodes(int *) const override;

	Corotator *getCorotator(CoordSet &, double*, int, int) override;
	int getTopNumber() override;

	void setPressure(PressureBCond *_pbc) override { pbc = _pbc; }
	PressureBCond* getPressure() override { return pbc; }
	void computePressureForce(CoordSet&, Vector& elPressureForce,
	                          GeomState *gs = 0, int cflg = 0, double t = 0) override;

	void getThermalForce(CoordSet &cs, Vector &ndTemps,
	                     Vector &ThermalForce, int glflag, GeomState *gs) override;
	void computeDisp(CoordSet &cs, State &state, const InterpPoint &,
	                 double *res, GeomState *gs) override;
	void getFlLoad(CoordSet &cs,  const InterpPoint &, double *flF,
	               double *resF, GeomState *gs) override;
	void setOffset(double *o) override { offset = o; }

	void getVonMises(Vector& stress, Vector& weight,CoordSet &cs, Vector& elDisp,
	                 int strInd,int surface, double *ndTemps,
	                 double ylayer, double zlayer, int avgnum) override;
#ifdef USE_EIGEN3
	void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *,
	                                        CoordSet &cs, Vector &elDisp, int strInd, int surface,
	                                        double *, int avgnum, double ylayer, double zlayer) override;
#endif
	// Routines for the decomposer
	PrioInfo examine(int sub, MultiFront *) override;
	bool hasRot() override { return true; }

	int getMassType() override { return 2; } // both consistent and lumped

private:

	void    getLength(CoordSet&, double &length);
	void    updTransMatrix(CoordSet&, GeomState *gs, double t[3][3], double &len, double weight = 0.5);
	void    offsetAxis(FullSquareMatrix& mat);

};
#endif

