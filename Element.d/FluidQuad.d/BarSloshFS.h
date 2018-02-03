#ifndef _BARSLOSHFS_H_
#define _BARSLOSHFS_H_

#include <Element.d/Element.h>
#include <cstdio>

class BarSloshFS: public Element {

        int nn[2];
public:

        BarSloshFS(int*);

        Element *clone() override;

        void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(CoordSet&,double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;
        PrioInfo examine(int sub, MultiFront *) {
          fprintf(stderr,"BarSloshFS.h: PrioInfo examine is commented in Dec.d/ElemMFCheck.C.");
          return *(new PrioInfo);
        };
        int getTopNumber() override;
};

#endif
