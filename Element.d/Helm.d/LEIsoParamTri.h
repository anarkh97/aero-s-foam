#ifndef _LEISOPARAMTRI_H_
#define _LEISOPARAMTRI_H_

#include <complex>
using std::complex;

#include <Element.d/Element.h>

class LEIsoParamTri: public Element {

        int order;
	int *nn;
        LEIsoParamTri(const LEIsoParamTri& e);

public:
	LEIsoParamTri(int,int*);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        double getMass(CoordSet&);
        double weight(CoordSet&, double *, int);
        double weightDerivativeWRTthickness(CoordSet&, double *, int);

	Element *clone();
	void renum(int *);
        void renum(EleRenumMap&);
	void markDofs(DofSetArray &);
//	int getTopNumber() {return 195;}
	int numTopNodes() {return (order*(order+1))/2;}
        int* dofs(DofSetArray &, int *p=0);
        int numDofs() { return (order*(order+1)); }
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
