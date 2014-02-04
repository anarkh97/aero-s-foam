#ifndef _FELIPPASHELL_H_
#define _FELIPPASHELL_H_

#ifdef USE_EIGEN3
#include <Element.d/Element.h>
#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>
#include <Corotational.d/Shell3Corotator.h>

class FelippaShell : public Element, 
                     public ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle>, 
                     public Shell3Corotator
{
        int      nn[3];
        int      type;
        double  *cFrame;
        PressureBCond *pbc;

public:
	FelippaShell(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
#ifdef USE_EIGEN3
        FullSquareMatrix getStiffnessThicknessSensitivity(CoordSet& cs, Vector &elDisp, double *d, int flg=1, int senMethod=0);
#endif
        FullSquareMatrix massMatrix(CoordSet& cs,double *mel,int cmflg=1);

        void getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
	                     GeomState *gs);

        void getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
                         Vector &elDisp, int strInd, int surface = 0,
                         double *ndTemps = 0, double ylayer = 0, double zlayer = 0,
                         int avgnum = 0);
#ifdef USE_EIGEN3
        void getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs, Vector &elDisp, 
                                             int strInd, int surface, int senMethod = 1, double * = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);

        void getVonMisesThicknessSensitivity(ComplexVector &dStdThick, ComplexVector &weight, CoordSet &cs, ComplexVector &elDisp, 
                                             int strInd, int surface, int senMethod = 1, double * = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);

        void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, CoordSet &cs, Vector &elDisp, 
                                                int strInd, int surface, int senMethod = 1, double * = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);

        void getVonMisesDisplacementSensitivity(GenFullM<DComplex> &dStdDisp, ComplexVector &weight, CoordSet &cs, ComplexVector &elDisp, 
                                                int strInd, int surface, int senMethod = 1, double * = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);

        void getStiffnessThicknessSensitivity(GenFullM<double> &dStiffdThick, CoordSet &cs, Vector &elDisp, 
                                               int senMethod, double *, double ylayer = 0, double zlayer = 0);

#endif

        void getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
                         Vector &elDisp, int strInd, int surface = 0,
                         double *ndTemps = 0);

        void setProp(StructProp *p, bool myProp = false);
        void setCompositeData(int _type, int nlays, double *lData,
                              double *coefs, double *frame);
        double * setCompositeData2(int _type, int nlays, double *lData,
                                   double *coefs, CoordSet &cs, double theta);
        void setMaterial(NLMaterial *);
        int numStates();

	void markDofs(DofSetArray &);
        int* dofs(DofSetArray &, int *p=0);
        int  numDofs();
	int  getTopNumber();

        int  numNodes();
        int* nodes(int * = 0);
	double getMass(CoordSet &);
        double weight(CoordSet&, double *, int);
        double weightDerivativeWRTthickness(CoordSet&, double *, int, int senMethod = 1);

        void computeDisp(CoordSet&, State &, const InterpPoint &,
                         double*, GeomState *gs=0);
        void getFlLoad(CoordSet &, const InterpPoint &, double *flF, 
                       double *resF, GeomState *gs=0);

        void setPressure(PressureBCond *_pbc) { pbc = _pbc; }
        PressureBCond* getPressure() { return pbc; }
        void computePressureForce(CoordSet&, Vector& elPressureForce,
                                  GeomState *gs = 0, int cflg = 0, double t = 0);

        // Nonlinear
        Corotator* getCorotator(CoordSet&, double*, int, int);
        void getStiffAndForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
                              FullSquareMatrix &elK, double *f, double dt, double t);
        void getInternalForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
                               FullSquareMatrix &elK, double *f, double dt, double t);
        void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0);
        void initStates(double *);
        double getDissipatedEnergy(GeomState &, CoordSet &);

        // Routines for the decomposer
        PrioInfo examine(int sub, MultiFront *);

	bool hasRot() { return true; }

        int getMassType() { return 0; } // lumped only

};
#endif
#endif

