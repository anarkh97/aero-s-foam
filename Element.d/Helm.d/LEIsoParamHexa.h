#ifndef _LEISOPARAMHEXA_H_
#define _LEISOPARAMHEXA_H_

#include <complex>
using std::complex;

#include <Element.d/Element.h>

class LEIsoParamHexa: public Element {

        int order;
	int *nn;
        LEIsoParamHexa(const LEIsoParamHexa& e);

public:
	LEIsoParamHexa(int,int*);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        double getMass(CoordSet&);
        void   getGravityForce(CoordSet&,double *gravity,Vector &force,
                                       int gravflg, GeomState *gs=0);


	Element *clone();
	void renum(int *);
	void markDofs(DofSetArray &);
	int getTopNumber() {return 195;}
	int numTopNodes() {return order*order*order;}
        int* dofs(DofSetArray &, int *p=0);
        int numDofs() { return 3*order*order*order; }
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
