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

public:
	FelippaShell(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet& cs,double *mel,int cmflg=1);

        void getGravityForce(CoordSet&,double *gravity, Vector&, int gravflg,
	                     GeomState *gs);

        void getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
                         Vector &elDisp, int strInd, int surface = 0,
                         double *ndTemps = 0, double ylayer = 0, double zlayer = 0,
                         int avgnum = 0);

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

        void computeDisp(CoordSet&, State &, const InterpPoint &,
                         double*, GeomState *gs=0);
        void getFlLoad(CoordSet &, const InterpPoint &, double *flF, 
                       double *resF, GeomState *gs=0);
        void computePressureForce(CoordSet&, Vector& elPressureForce,
                                  GeomState *gs = 0, int cflg = 0, double t = 0);

        // Nonlinear
        Corotator* getCorotator(CoordSet&, double*, int, int);
        void getStiffAndForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
                              FullSquareMatrix &elK, double *f, double dt, double t);
        void getInternalForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
                               FullSquareMatrix &elK, double *f, double dt, double t);
        void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0);

        // Routines for the decomposer
        PrioInfo examine(int sub, MultiFront *);

	bool hasRot() { return true; }
	//double weight() { return 3; }
        //double trueWeight() { return 3; }

        int getMassType() { return 0; } // lumped only

};
#endif
#endif

