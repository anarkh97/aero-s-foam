#ifndef _THERMSOPARAMHEXA_H_
#define _THERMSOPARAMHEXA_H_

#include <Element.d/Element.h>

#include <complex>
using std::complex;


class ThermIsoParamHexa: public Element {

        int order;
	int *nn;
        ThermIsoParamHexa(const ThermIsoParamHexa& e);

public:
	ThermIsoParamHexa(int,int*);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        double getMass(CoordSet&);

	Element *clone();
	void renum(int *);
        void renum(EleRenumMap&);
	void markDofs(DofSetArray &);
	int getTopNumber() {return 195;}
	int numTopNodes() {return order*order*order;}
        int* dofs(DofSetArray &, int *p=0);
        int numDofs() { return order*order*order; }
        int numNodes();
        int* nodes(int * = 0);

        PrioInfo examine(int sub, MultiFront *mf);

	double weight() {
	  return order;
	}
	
	double trueWeight() {
	  return weight();
	}
};
#endif
