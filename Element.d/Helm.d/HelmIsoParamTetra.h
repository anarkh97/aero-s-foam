#ifndef _HELMISOPARAMTETRA_H_
#define _HELMISOPARAMTETRA_H_

#include <complex>
using std::complex;

#include <Element.d/Helm.d/HelmElement.h>

class HelmIsoParamTetra: public HelmElement, public Element {

        int order;
	int *nn;
        HelmIsoParamTetra(const HelmIsoParamTetra& e);

public:
	HelmIsoParamTetra(int,int*);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        FullSquareMatrixC stiffness(CoordSet&, complex<double> *d);
        FullSquareMatrixC massMatrix(CoordSet&, complex<double> *d);
        double getMass(CoordSet&);

	Element *clone();
	void renum(int *);
        void renum(EleRenumMap&);
	void markDofs(DofSetArray &);
        int getTopNumber() {return 196;}
        int numTopNodes() {return (order*(order+1)*(order+2))/6;}
        int* dofs(DofSetArray &, int *p=0);
        int numDofs() { return (order*(order+1)*(order+2))/6; }
        int numNodes();
        int* nodes(int * = 0);
	void addFaces(PolygonSet *pset);

        PrioInfo examine(int sub, MultiFront *mf);
        int nDecFaces() { return 4;}
        int getDecFace(int iFace, int *fn);


//	double weight() {
// RT: using a condensed dofs number
//           double w = 
//                      (6.0*double(order)*double(order+1) -
//                       29.0*double(order) +
//                       23.0)/6.0 ;
//           return w;
//	}
	
};
#endif
