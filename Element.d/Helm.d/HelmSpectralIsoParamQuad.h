#ifndef _HELMSPECTRALISOPARAMQUAD_H_
#define _HELMSPECTRALISOPARAMQUAD_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmSpectralIsoParamQuad: public HelmElement, public Element {

        int order;
	int *nn;
        HelmSpectralIsoParamQuad(const HelmSpectralIsoParamQuad& e);

public:
	HelmSpectralIsoParamQuad(int,int*);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        double getMass(CoordSet&);

	Element *clone();
	void renum(int *);
	void markDofs(DofSetArray &);
	int getTopNumber() {return 163;}
	int numTopNodes() {return order*order;}
        int* dofs(DofSetArray &, int *p=0);
        int numDofs() { return order*order; }
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
