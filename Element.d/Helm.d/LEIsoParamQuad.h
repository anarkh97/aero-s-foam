#ifndef _LEISOPARAMQUAD_H_
#define _LEISOPARAMQUAD_H_

#include <complex>
using std::complex;

#include <Element.d/Element.h>

class LEIsoParamQuad: public Element {

        int order;
	int *nn;
        LEIsoParamQuad(const LEIsoParamQuad& e);

public:
	LEIsoParamQuad(int,int*);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        double getMass(CoordSet&);

	Element *clone();
	void renum(int *);
	void markDofs(DofSetArray &);
//	int getTopNumber() {return 195;}
	int numTopNodes() {return order*order;}
        int* dofs(DofSetArray &, int *p=0);
        int numDofs() { return 2*order*order; }
        int numNodes();
        int* nodes(int * = 0);

//        PrioInfo examine(int sub, MultiFront *mf);

	double weight() {
	  return order;
	}
	
	double trueWeight() {
	  return weight();
	}
};
#endif
