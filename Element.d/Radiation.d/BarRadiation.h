#ifndef _BARRADIATION_H_
#define _BARRADIATION_H_

#include <Element.d/Element.h>

class BarRadiation: public virtual Element {

        int nn[2];
public:

	BarRadiation(int*);

        Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&,double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double           getMass(CoordSet&);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        Corotator* getCorotator(CoordSet &, double*, int, int);
        
        //void getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &ThermalForce, int glflag,
        //                   GeomState *gs) { ThermalForce.zero(); }

        int              numNodes();
        int*             nodes(int * = 0);
	PrioInfo examine(int sub, MultiFront *);
	int 		getTopNumber();
};
#endif
