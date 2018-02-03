#ifndef _PENTACONTACT_H_
#define _PENTACONTACT_H_

#include <Element.d/Element.h>

class PentaContact: public virtual Element {

        int nn[6];
public:
        PentaContact(int*);

        Element *clone() override;

        void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg) const override;
        FullSquareMatrix massMatrix(const CoordSet& cs,double *d, int cmflg) const override;
        double           getMass(const CoordSet&) const;

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;
        PrioInfo examine(int sub, MultiFront *) override;
        int getTopNumber() override;
};

#endif
