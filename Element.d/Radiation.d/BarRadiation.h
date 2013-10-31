#ifndef _BARRADIATION_H_
#define _BARRADIATION_H_

#include <Element.d/Element.h>

class BarRadiation: public virtual Element {

        int nn[2];
        double *f;
public:
	BarRadiation(int*);
        ~BarRadiation();

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
        
        int              numNodes();
        int*             nodes(int * = 0);
	PrioInfo examine(int sub, MultiFront *);
	int 		getTopNumber();

        // these functions are used to get the contribution of the radiation element
        // to the rhs for linear analyses.
        bool isRadiationElement() { return true; }
        void computePressureForce(CoordSet& cs,Vector& elPressureForce,
                                  GeomState *gs = 0, int cflg = 0, double t = 0);
};

#endif
