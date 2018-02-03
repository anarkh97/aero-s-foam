#ifndef _TORSPRING_H_
#define _TORSPRING_H_

#include <Element.d/Element.h>

class TorSpring : public Element {

	int nn[1];
public:

	explicit TorSpring(int*);

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(CoordSet& cs, double *kel, int flg=1) override;
	FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1) override;

	void markDofs(DofSetArray &) override;
	int* dofs(DofSetArray &, int *p) override;
	int numDofs() const override;

	int  numNodes() const override;
	int* nodes(int *) const override;
	bool isSafe() override {return true;}
	bool isSpring() override {return true;}
	int getTopNumber() override {return 111;}
	PrioInfo examine(int sub, MultiFront *) override;
	Corotator *getCorotator(CoordSet &, double*, int, int) override;
};
#endif
