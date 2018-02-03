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
        double getMassThicknessSensitivity(CoordSet&);

        Element *clone() override;
        void renum(int *) override;
        void renum(EleRenumMap&) override;
        void markDofs(DofSetArray &) override;
//	int getTopNumber() override {return 195;}
        int numTopNodes() {return (order*(order+1))/2;}
        int* dofs(DofSetArray &, int *p) override;
        int numDofs() const { return (order*(order+1)); }
        int numNodes() const override;
        int* nodes(int * = 0) const override;

//        PrioInfo examine(int sub, MultiFront *mf) override;

};
#endif
