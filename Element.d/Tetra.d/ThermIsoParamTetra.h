#ifndef _THERMISOPARAMTETRA_H_
#define _THERMISOPARAMTETRA_H_

#include <Element.d/Element.h>

class ThermIsoParamTetra: public Element {

        int order;
	int *nn;
        ThermIsoParamTetra(const ThermIsoParamTetra& e);

public:
	ThermIsoParamTetra(int,int*);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        double getMass(CoordSet&);

	Element *clone();
	void renum(int *);
        void renum(EleRenumMap&);
	void markDofs(DofSetArray &);
        int* dofs(DofSetArray &, int *p=0);
        int numDofs() { return (order*(order+1)*(order+2))/6; }
        int numNodes();
        int* nodes(int * = 0);

        PrioInfo examine(int sub, MultiFront *mf);
        int getTopNumber();
};
#endif
