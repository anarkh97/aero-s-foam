#ifndef _HELMISOPARAMQUAD_H_
#define _HELMISOPARAMQUAD_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmIsoParamQuad: public HelmElement, public Element {

        int order;
	int *nn;
        HelmIsoParamQuad(const HelmIsoParamQuad& e);

public:
	HelmIsoParamQuad(int,int*);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        FullSquareMatrixC stiffness(CoordSet&, complex<double> *d);
        FullSquareMatrixC massMatrix(CoordSet&, complex<double> *d);
        double getMass(CoordSet&);

	Element *clone();
	void renum(int *);
        void renum(EleRenumMap&);
	void markDofs(DofSetArray &);
//	int getTopNumber() {return 195;}
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
