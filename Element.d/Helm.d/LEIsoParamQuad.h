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

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
        FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
        double getMass(const CoordSet&) const;
        double getMassThicknessSensitivity(CoordSet&);

        Element *clone() override;
        void renum(int *) override;
        void renum(EleRenumMap&) override;
        void markDofs(DofSetArray &) const override;
//	int getTopNumber() override {return 195;}
        int numTopNodes() {return order*order;}
        int* dofs(DofSetArray &, int *p) const override;
        int numDofs() const override { return 2*order*order; }
        int numNodes() const override;
        int* nodes(int * = 0) const override;

//        PrioInfo examine(int sub, MultiFront *mf) override;

};
#endif
