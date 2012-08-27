#ifndef _FOURNODEQUAD_H_
#define _FOURNODEQUAD_H_

#include <Element.d/Element.h>

class FourNodeQuad: virtual public Element {
 protected:
	int nn[4];
public:
	FourNodeQuad(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double           getMass(CoordSet&);

        void             getGravityForce(CoordSet&,double *gravity, Vector& f, int gravflg,
	                                 GeomState *gs);
        void             getVonMises (Vector &stress, Vector &weight,
                                      CoordSet &cs, Vector &elDisp, 
                                      int strInd, int surface=0,
                                      double *ndTemps=0,
				      double ylayer=0.0, double zlayer=0.0, int avgnum=0);

        void             getAllStress(FullM &stress, Vector &weight,
                                      CoordSet &cs, Vector &elDisp,
                                      int strInd, int surface=0,
                                      double *ndTemps=0);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

	int getTopNumber();


        void computeDisp(CoordSet&cs, State &state, const InterpPoint &,
                         double*res, GeomState *gs);
        void getFlLoad(CoordSet &, const InterpPoint &,  double *flF,
                       double *resF, GeomState *gs=0);
        void getThermalForce(CoordSet &cs, Vector &ndTemps, 
                             Vector &ThermalForce, int glflag, GeomState *gs=0);
        int dim() { return 2; }

        void spy();
	
	// DEC
	PrioInfo examine(int sub, MultiFront *);
};
#endif

