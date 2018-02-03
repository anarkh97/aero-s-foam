#ifndef _TRIANGLECONVEC_H_
#define _TRIANGLECONVEC_H_

#include <Element.d/Element.h>

class TriangleConvec: public Element {

	int nn[4];
public:
	TriangleConvec(int*);

	Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg) const override;
        FullSquareMatrix massMatrix(const CoordSet& cs,double *d, int cmflg) const override;
        double           getMass(const CoordSet&) const;

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int             numNodes() const override;
        int * nodes(int *) const override;
	int getTopNumber() override;
	PrioInfo examine(int sub, MultiFront *) override;
};

#endif
