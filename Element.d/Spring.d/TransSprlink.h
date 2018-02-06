#ifndef _TRANSSPRLINK_H_
#define _TRANSSPRLINK_H_

#include <Element.d/Element.h>

class TransSprlink : virtual public Element 
{
    int nn[2];

  public:

	explicit TransSprlink(int*);

    Element *clone() override;

    void renum(int *) override;
	void renum(EleRenumMap&) override;

    FullSquareMatrix stiffness(const CoordSet& cs, double* d, int flg = 1) const override;
    FullSquareMatrix massMatrix(const CoordSet& cs, double* mel, int cmflg=1) const override;

    void markDofs(DofSetArray&) override;
    int* dofs(DofSetArray&, int*) override;
    int numDofs() const override;

    int numNodes() const override;
    int* nodes(int*) const override;
    Corotator* getCorotator(CoordSet&, double*, int, int) override;

    int getTopNumber() override;
    bool isSafe() const override { return false; }
    bool isSpring() override { return true; }
    PrioInfo examine(int sub, MultiFront*) override;
};

#endif
