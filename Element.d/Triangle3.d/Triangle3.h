#ifndef _TRIANGLE3_H_
#define _TRIANGLE3_H_

#include <Element.d/Element.h>

class Triangle3: public Element {

	int nn[3];
public:
	Triangle3(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg = 1);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);
        double           getMass(CoordSet&);
        double weight(CoordSet&, double *, int);
        double weightDerivativeWRTthickness(CoordSet&, double *, int, int=1);

        void             getGravityForce(CoordSet&,double *gravity, Vector& f, int gravflg,
	                                 GeomState *gs);

        virtual void     getVonMises (Vector &stress, Vector &weight, 
                                      CoordSet &cs, Vector &elDisp, 
                                      int strInd, int surface=0,
                                      double *ndTemps=0,
				      double ylayer=0.0, double zlayer=0.0, int avgnum=0);

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

	//bool hasRot() { return true; }
	//double weight() { return 3; }
        //double trueWeight() { return 3; }
};
#endif

