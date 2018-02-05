#ifndef _LEISOPARAMHEXA_H_
#define _LEISOPARAMHEXA_H_

#include <complex>
using std::complex;

#include <Element.d/Element.h>

class LEIsoParamHexa: public Element {

        int order;
	int *nn;
        LEIsoParamHexa(const LEIsoParamHexa& e);

public:
	LEIsoParamHexa(int,int*);

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const;
        FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const;
        void aRubberStiffnessDerivs(CoordSet&, complex<double> *d, int n, double omega);
        double getMass(const CoordSet&) const;
        void   getGravityForce(CoordSet&,double *gravity,Vector &force,
                                       int gravflg, GeomState *gs=0);


	Element *clone() override;
	void renum(int *) override;
        void renum(EleRenumMap&) override;
        
	void markDofs(DofSetArray &) override;
	int getTopNumber() override {return 195;}
	int numTopNodes() {return order*order*order;}
        int* dofs(DofSetArray &, int *p) override;
        int numDofs() const override { return 3*order*order*order; }
        int numNodes() const override;
        int* nodes(int * = 0) const override;

        PrioInfo examine(int sub, MultiFront *mf) override;
        int nDecFaces() const override { return 6;}
        int getDecFace(int iFace, int *fn);

};
#endif
