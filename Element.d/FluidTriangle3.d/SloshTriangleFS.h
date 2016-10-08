#ifndef _SLOSHTRIANGLEFS_H_
#define _SLOSHTRIANGLEFS_H_

#include <Element.d/Element.h>
#include <cstdio>

class SloshTriangleFS: public Element {

        int nn[3];
public:
        SloshTriangleFS(int*);

        Element *clone();

        void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg = 1);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);
        double           getArea(CoordSet&);
        bool             isSloshingElement() { return true; }

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
        int              getTopNumber();

        PrioInfo examine(int sub, MultiFront *) {
          fprintf(stderr,"SloshTriangleFS.h: PrioInfo examine is commented in Dec.d/ElemFSCheck.C");
          return *(new PrioInfo);
        };
};

#endif
