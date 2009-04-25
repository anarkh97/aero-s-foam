#ifndef _HELMISOPARAMTRI_H_
#define _HELMISOPARAMTRI_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmIsoParamTri: public HelmElement, public Element {

        int order;
	int *nn;
        HelmIsoParamTri(const HelmIsoParamTri& e);

public:
	HelmIsoParamTri(int,int*);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        FullSquareMatrixC stiffness(CoordSet&, complex<double> *d);
        FullSquareMatrixC massMatrix(CoordSet&, complex<double> *d);
        double getMass(CoordSet&);

	Element *clone();
	void renum(int *);
	void markDofs(DofSetArray &);
//	int getTopNumber() {return 195;}
	int numTopNodes() {return (order*(order+1))/2;}
        int* dofs(DofSetArray &, int *p=0);
        int numDofs() { return (order*(order+1))/2; }
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
