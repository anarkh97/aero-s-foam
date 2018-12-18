#ifndef _CTCVIRTUALELT_H_
#define _CTCVIRTUALELT_H_

#include <Element.d/Element.h>

class CtcVirtualElt : public Element {

        int nn[2];

public:

	CtcVirtualElt(int,int*);

        Element *clone() override;

	void renum(const int *) override;
        void renum(EleRenumMap&) override;
	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	 int numDofs() const override;
	
	int numNodes() const override;
	int* nodes(int * = 0) const override;
	int getTopNumber() const override;


	FullSquareMatrix stiffness(const CoordSet&,double *kel, int flg=1) const override;
        FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
};

#endif
