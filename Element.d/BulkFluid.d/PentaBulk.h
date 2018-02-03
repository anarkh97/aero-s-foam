#ifndef _PENTABULK_H_
#define _PENTABULK_H_

#include <Element.d/Element.h>

class PentaBulk: public virtual Element {

        int nn[5];
public:
        PentaBulk(int*);

        Element *clone() override;

        void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg) override;
        FullSquareMatrix massMatrix(CoordSet& cs,double *d, int cmflg) override;
        double getMass(CoordSet&);

        void  markDofs(DofSetArray &);
        int*  dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;
        int getTopNumber() override;
        PrioInfo examine(int sub, MultiFront *) override;
};

#endif
