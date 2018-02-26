#ifndef _BARCONVEC_H_
#define _BARCONVEC_H_

#include <Element.d/Element.h>

class BarConvec: public Element {

	int nn[2];
public:

	explicit BarConvec(int*);

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&,double *kel, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
	double getMass(const CoordSet&) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int numNodes() const override;
	int* nodes(int *) const override;
	PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() override;
};
#endif
