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

	FullSquareMatrix stiffness(const CoordSet& cs, double *kel, int flg=1) const override;
	FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg=1) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p) const override;
	int numDofs() const override;

	int  numNodes() const override;
	int* nodes(int *) const override;
	bool isSafe() const override {return true;}
	bool isSpring() const override {return true;}
	int getTopNumber() override {return 111;}
	PrioInfo examine(int sub, MultiFront *) override;
	Corotator *getCorotator(CoordSet &, double*, int, int) override;
};
#endif
