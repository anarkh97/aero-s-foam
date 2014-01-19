#ifndef _LEISOPARAMTETRA_H_
#define _LEISOPARAMTETRA_H_

#include <complex>
using std::complex;

#include <Element.d/Element.h>

class LEIsoParamTetra: public Element {

        int order;
	int *nn;
        LEIsoParamTetra(const LEIsoParamTetra& e);

public:
	LEIsoParamTetra(int,int*);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *d, int cmflg=1);
        double getMass(CoordSet&);

	Element *clone();
	void renum(int *);
        void renum(EleRenumMap&);
	void markDofs(DofSetArray &);
        int* dofs(DofSetArray &, int *p=0);
        int numDofs() { return (order*(order+1)*(order+2))/2; }
        int numNodes();
        int* nodes(int * = 0);

        PrioInfo examine(int sub, MultiFront *mf);
        int nDecFaces() { return 4;}
        int getDecFace(int iFace, int *fn);

};
#endif
