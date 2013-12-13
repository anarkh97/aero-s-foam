#ifndef _HELMISOPARAMHEXA_H_
#define _HELMISOPARAMHEXA_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmIsoParamHexa: public HelmElement, public Element {

        int order;
	int *nn;
        HelmIsoParamHexa(const HelmIsoParamHexa& e);

public:
	HelmIsoParamHexa(int,int*);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        FullSquareMatrixC stiffness(CoordSet&, complex<double> *d);
        FullSquareMatrixC massMatrix(CoordSet&, complex<double> *d);
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
	void addFaces(PolygonSet *pset);

        PrioInfo examine(int sub, MultiFront *mf);

};
#endif
