#ifndef _FOURNODEQUAD_H_
#define _FOURNODEQUAD_H_

#include <Element.d/Element.h>
#include <Element.d/Quad4.d/QuadElementTemplate.hpp>

class FourNodeQuad: virtual public Element,
                    public QuadElementTemplate<double>
{
protected:
	int nn[4];
public:
	explicit FourNodeQuad(int*);

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(CoordSet&, double *d, int flg) override;
	FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg) override;
	double           getMass(CoordSet&) override;
	double getMassSensitivityWRTthickness(CoordSet&);
	double weight(CoordSet&, double *) override;
	double weightDerivativeWRTthickness(CoordSet&, double *, int=1);

	void             getGravityForce(CoordSet&,double *gravity, Vector& f, int gravflg,
	                                 GeomState *gs) override;
	void getGravityForceSensitivityWRTthickness(CoordSet&,double *gravity, int senMethod, Vector& f, int gravflg,
	                                            GeomState *gs);
	void             getVonMises (Vector &stress, Vector &weight,
	                              CoordSet &cs, Vector &elDisp,
	                              int strInd, int surface,
	                              double *ndTemps,
	                              double ylayer, double zlayer, int avgnum) override;
	void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double>*, CoordSet &cs,
	                                        Vector &elDisp, int strInd, int surface, int senMethod, double *ndTemps,
	                                        int avgnum, double ylayer, double zlayer);

	void             getAllStress(FullM &stress, Vector &weight,
	                              CoordSet &cs, Vector &elDisp,
	                              int strInd, int surface,
	                              double *ndTemps) override;

	void             markDofs(DofSetArray &) override;
	int*             dofs(DofSetArray &, int *p) override;
	int              numDofs() const override;

	int              numNodes() const override;
	int*             nodes(int *) const override;

	int getTopNumber() override;


	void computeDisp(CoordSet&cs, State &state, const InterpPoint &,
	                 double*res, GeomState *gs) override;
	void getFlLoad(CoordSet &, const InterpPoint &,  double *flF,
	               double *resF, GeomState *gs) override;
	void getThermalForce(CoordSet &cs, Vector &ndTemps,
	                     Vector &ThermalForce, int glflag, GeomState *gs) override;
	int dim() const override { return 2; }

	void spy();

	// DEC
	PrioInfo examine(int sub, MultiFront *) override;

	Corotator * getCorotator(CoordSet &, double*, int, int) override { return nullptr; }
};
#endif

