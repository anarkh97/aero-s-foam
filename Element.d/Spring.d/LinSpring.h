#ifndef _LINSPRING_H_
#define _LINSPRING_H_

#include <Element.d/Element.h>

class LinSpring : public Element {

	int nn[1];
public:

	explicit LinSpring(int*);

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet& cs, double *kel, int flg) const override;
	FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg) const override;

	void markDofs(DofSetArray &) override;
	int* dofs(DofSetArray &, int *p) override;
	int numDofs() const override;

	int numNodes() const override;
	int * nodes(int *) const override;
	bool isSafe() const override {return true;}
	bool isSpring() override {return true;}
	int getTopNumber() override {return 111;}
	PrioInfo examine(int sub, MultiFront *) override;
	Corotator *getCorotator(CoordSet &, double*, int, int) override;
};
#endif
