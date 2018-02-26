#ifndef    _TIMOSHENKOBEAM_H_
#define    _TIMOSHENKOBEAM_H_

#include <Element.d/Element.h>
#include <Element.d/Beam.d/BeamElementTemplate.hpp>

class GeomState;

class TimoshenkoBeam : public Element,
                       public BeamElementTemplate<double> {
	EFrame *elemframe;
	double oeframe[3][3];
	bool myElemFrame;
	int nn[3];
	double *iniOr;
	PressureBCond *pbc;

public:
	explicit TimoshenkoBeam(int *);

	~TimoshenkoBeam() override;

	TimoshenkoBeam *clone() override;

	void renum(int *) override;

	void renum(EleRenumMap &) override;

	void setFrame(EFrame *ef) override {
		elemframe = ef;
		myElemFrame = false;
	}

	const EFrame *getFrame() const override { return elemframe; }

	void buildFrame(CoordSet &) override;

	FullSquareMatrix stiffness(const CoordSet &cs, double *kel, int flg) const;

	FullSquareMatrix massMatrix(const CoordSet &cs, double *mel, int cmflg) const;

	void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs);

	double getMass(const CoordSet &cs) const;

	void getMassNodalCoordinateSensitivity(CoordSet &cs, Vector &dMassdx);

	double weight(CoordSet &cs, double *gravityAcceleration);

	void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet &cs, double *gravityAcceleration);

	void getGravityForce(CoordSet &, double *gravity, Vector &, int gravflg, GeomState *gs);

	void getGravityForceNodalCoordinateSensitivity(CoordSet &cs, double *gravityAcceleration,
	                                               GenFullM<double> &dGfdx, int gravflg, GeomState *geomState);

	void getIntrnForce(Vector &elForce, CoordSet &cs,
	                   double *elDisp, int forceIndex, double *ndTemps);

	void markDofs(DofSetArray &) const override;

	int *dofs(DofSetArray &, int *p) const override;

	int numDofs() const override;

	int numNodes() const override;

	int *nodes(int *) const override;

	Corotator *getCorotator(CoordSet &, double *, int, int);

	int getTopNumber() override;

	void setPressure(PressureBCond *_pbc) { pbc = _pbc; }

	PressureBCond *getPressure() { return pbc; }

	void computePressureForce(CoordSet &, Vector &elPressureForce,
	                          GeomState *gs, int cflg, double t);

	void getThermalForce(CoordSet &, Vector &, Vector &, int, GeomState *geomState);

	void getVonMises(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp,
	                 int strInd, int surface, double *ndTemps,
	                 double ylayer, double zlayer, int avgnum);

#ifdef USE_EIGEN3

	void
	getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *, CoordSet &cs,
	                                   Vector &elDisp, int strInd, int surface, double *ndTemps, int avgnum,
	                                   double ylayer, double zlayer) override;

	void getVonMisesDisplacementSensitivity(GenFullM<DComplex> &dStdDisp, ComplexVector &weight, CoordSet &cs,
	                                        ComplexVector &elDisp,
	                                        int strInd, int surface, double *ndTemps, int avgnum,
	                                        double ylayer,
	                                        double zlayer);

	void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs, Vector &elDisp,
	                                           int strInd, int surface,
	                                           double *ndTemps, int avgnum, double ylayer,
	                                           double zlayer);

#endif

	int getMassType() const override { return 0; } // lumped only
private:
	TimoshenkoBeam(const TimoshenkoBeam &);

	TimoshenkoBeam &operator=(const TimoshenkoBeam &);

	void updTransMatrix(CoordSet &, GeomState *gs, double t[3][3], double &len);

	// Routines for the decomposer
	PrioInfo examine(int sub, MultiFront *) override;

	bool hasRot() const override { return true; }
};

#endif
