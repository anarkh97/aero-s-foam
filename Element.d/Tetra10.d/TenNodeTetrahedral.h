#ifndef _TENNODETETRAHEDRAL_H_
#define _TENNODETETRAHEDRAL_H_

#include <Element.d/Element.h>

#include <complex>
using std::complex;

class TenNodeTetrahedral: public Element
{
	int nn[10];
	double *cCoefs;
	double *cFrame;
	NLMaterial *mat;

public:
	TenNodeTetrahedral(int*);
	~TenNodeTetrahedral();

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(CoordSet&, double *kel, int flg) override;
	FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg) override;
	void aRubberStiffnessDerivs(CoordSet&, complex<double> *d, int n, double omega) override;
	double getMass(CoordSet& cs) override;

	void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg, GeomState *gs) override;
	void getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &force, int glflag, GeomState *gs) override;

	void getVonMises(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
	                 int surface, double *ndTemps, double ylayer, double zlayer, int avgnum) override;

	void getAllStress(FullM &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
	                  int surface, double *ndTemps) override;

	void markDofs(DofSetArray &) override;
	int* dofs(DofSetArray &, int *p) override;
	 int numDofs() const override;

	int numNodes() const override;
	int* nodes(int *) const override;

	int getTopNumber() override;

	PrioInfo examine(int sub, MultiFront *) override;
	int nDecFaces() override { return 4; }
	int getDecFace(int iFace, int *fn) override;

	int getFace(int iFace, int *fn) override;

	void setCompositeData(int _type, int nlays, double *lData, double *coefs, double *frame) override
	{ cCoefs = coefs; cFrame = frame; }

	double* setCompositeData2(int, int, double*, double*, CoordSet&, double) override
	{ fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
				"              for TenNodeTetrahedral el.\n"); return (double *) 0;
	}
	void getCFrame(CoordSet &cs, double cFrame[3][3]) const override;

	void getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
	                      int surface=0, double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0);

	void getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
	                       Vector &elDisp, int strInd, int surface=0, double *ndTemps=0);

	void setMaterial(NLMaterial *) override;
	int numStates() override;
	void initStates(double *st) override;
	Corotator *getCorotator(CoordSet &cs, double *kel, int=2, int=2) override;
};

#endif
