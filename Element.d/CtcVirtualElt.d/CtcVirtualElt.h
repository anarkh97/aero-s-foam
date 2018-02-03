#ifndef _CTCVIRTUALELT_H_
#define _CTCVIRTUALELT_H_

#include <Element.d/Element.h>

class CtcVirtualElt : public Element {

        int nn[2];

public:

	CtcVirtualElt(int,int*);

        Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) override;
	int* dofs(DofSetArray &, int *p) override;
	 int numDofs() const override;
	
	int numNodes() const override;
	int* nodes(int * = 0) const override;
	int getTopNumber() override;


	FullSquareMatrix stiffness(CoordSet&,double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
};

#endif
