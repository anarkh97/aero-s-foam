#ifndef _SLOSHQUADGAL_H_
#define _SLOSHQUADGAL_H_

#include <Element.d/Element.h>
#include <cstdio>

class SloshQuadGal: public Element {

        int nn[4];
public:
        SloshQuadGal(int*);

        Element *clone() override;

        void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flag=1) const;
        FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const override;
        double           getMass(const CoordSet&) const;

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;
        int getTopNumber() override;
        PrioInfo examine(int sub, MultiFront *) {
          fprintf(stderr,"SloshQuadGal.h: PrioInfo examine is commented in Dec.d/ElemMFCheck.C");
          return *(new PrioInfo);
        };

        void computeSloshDisp(Vector& heatflux, CoordSet &cs, Vector& elTemp, 
                              int hflInd);
        void computeSloshDispAll(Vector& heatflux, CoordSet &cs, Vector& elTemp);
};

#endif
