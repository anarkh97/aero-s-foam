#ifndef _BARCONVEC_H_
#define _BARCONVEC_H_

#include <Element.d/Element.h>

class BarConvec: public Element {

        int nn[2];
public:

	BarConvec(int*);

        Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(CoordSet&,double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double           getMass(CoordSet&);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int             numNodes() const override;
        int*             nodes(int * = 0) const;
	PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() override;
};
#endif
