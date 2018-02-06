#ifndef _ROTNSPRLINK_H_
#define _ROTNSPRLINK_H_

#include <Element.d/Element.h>

class RotnSprlink : public Element 
{
    int nn[2];

  public:

    RotnSprlink(int*);

    Element *clone() override;

    void renum(int *) override;
        void renum(EleRenumMap&) override;

    FullSquareMatrix stiffness(const CoordSet& cs, double* kel, int flg = 1) const;
    FullSquareMatrix massMatrix(const CoordSet& cs, double* mel, int cmflg = 1) const;

    void markDofs(DofSetArray&);
    int* dofs(DofSetArray&, int* = 0);
     int numDofs() const override;

    int numNodes() const override;
    int* nodes(int* = 0) const override;
    Corotator* getCorotator(CoordSet&, double*, int, int);

    int getTopNumber() override;
    bool isSafe() const override { return false; }
    bool isSpring() { return true; }
    PrioInfo examine(int sub, MultiFront*);
	
};
#endif
