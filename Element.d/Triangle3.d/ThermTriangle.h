#ifndef _THERMTRIANGLE_H_
#define _THERMTRIANGLE_H_

#include <Element.d/Element.h>

class ThermTriangle: public Element {

	int nn[3];
public:
	ThermTriangle(int*);

	Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix  stiffness(const CoordSet& cs, double *d, int flg = 1) const;
        FullSquareMatrix  massMatrix(const CoordSet& cs, double *mel, int cmflg=1) const;
        double            getMass(const CoordSet&) const;


	void markDofs(DofSetArray &) const override;
        int* dofs(DofSetArray &, int *p=0) const override;
         int numDofs() const override;

        int             numNodes() const override;
        int * nodes(int *) const override;
	int getTopNumber() override;

        void computeTemp(CoordSet&cs, State &state, double gp[2], double*res);
        void getFlFlux(double gp[2], double *flF, double *resF);
	PrioInfo examine(int sub, MultiFront *) override;
        Corotator * getCorotator(CoordSet &, double*, int, int) { return 0; }
};
#endif

