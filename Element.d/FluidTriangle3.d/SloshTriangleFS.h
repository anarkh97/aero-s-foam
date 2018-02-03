#ifndef _SLOSHTRIANGLEFS_H_
#define _SLOSHTRIANGLEFS_H_

#include <Element.d/Element.h>
#include <cstdio>

class SloshTriangleFS: public Element {

        int nn[3];
public:
        SloshTriangleFS(int*);

        Element *clone() override;

        void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg = 1) const override;
        FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg=1) const override;
        double getArea(const CoordSet&) const;
        bool isSloshingElement() { return true; }

        void markDofs(DofSetArray &) override;
        int* dofs(DofSetArray &, int *p=0) override;
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;
        int getTopNumber() override;

        PrioInfo examine(int sub, MultiFront *) {
          fprintf(stderr,"SloshTriangleFS.h: PrioInfo examine is commented in Dec.d/ElemFSCheck.C");
          return *(new PrioInfo);
        };
};

#endif
