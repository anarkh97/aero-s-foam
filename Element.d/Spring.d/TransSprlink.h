#ifndef _TRANSSPRLINK_H_
#define _TRANSSPRLINK_H_

#include <Element.d/Element.h>

class TransSprlink : virtual public Element 
{
    int nn[2];

  public:

    TransSprlink(int*);

    Element *clone();

    void renum(int *);

    FullSquareMatrix stiffness(CoordSet& cs, double* d, int flg = 1);
    FullSquareMatrix massMatrix(CoordSet& cs, double* mel, int cmflg=1);

    void markDofs(DofSetArray&);
    int* dofs(DofSetArray&, int* = 0);
    int numDofs();

    int numNodes();
    int* nodes(int* = 0);
    Corotator* getCorotator(CoordSet&, double*, int, int);

    int getTopNumber();
    bool isSafe() { return false; }
    bool isSpring() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
