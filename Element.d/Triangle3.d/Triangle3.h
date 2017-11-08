#ifndef _TRIANGLE3_H_
#define _TRIANGLE3_H_

#include <Element.d/Element.h>
#include <Element.d/Triangle3.d/Triangle3ElementTemplate.hpp>

class Triangle3: public Element,
                 public Triangle3ElementTemplate<double>
{

	int nn[3];
public:
	Triangle3(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg = 1);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);
        double           getMass(CoordSet&);
        double getMassThicknessSensitivity(CoordSet&);

        void             getGravityForce(CoordSet&,double *gravity, Vector& f, int gravflg,
	                                       GeomState *gs);
        void getGravityForceThicknessSensitivity(CoordSet&,double *gravity, Vector& f, int gravflg,
                                                 GeomState *gs);

        virtual void     getVonMises (Vector &stress, Vector &weight, 
                                      CoordSet &cs, Vector &elDisp, 
                                      int strInd, int surface=0,
                                      double *ndTemps=0,
                                      double ylayer=0.0, double zlayer=0.0, int avgnum=0);

        void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp,
                                                CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                double *ndTemps, int avgnum, double ylayer, double zlayer);

        virtual void     getAllStress(FullM &stress, Vector &weight, 
                                      CoordSet &cs, Vector &elDisp, 
                                      int strInd, int surface=0,
                                      double *ndTemps=0);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

	int getTopNumber();
	
	// Routines for the decomposer
        PrioInfo examine(int sub, MultiFront *);

};
#endif

