#ifndef _SLOSHTETRA_H_
#define _SLOSHTETRA_H_

#include <Element.d/Element.h>
#include <cstdio>

class SloshTetra: public Element {

        int nn[4];
public:
        SloshTetra(int*);

        Element *clone() override;

        void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg=1) const;
        FullSquareMatrix massMatrix(const CoordSet&,double *mel, int cmflg=1) const;
        bool             isSloshingElement() { return true; }
        double           getMass(const CoordSet& cs) const;
        void markDofs(DofSetArray &) const override;
        int* dofs(DofSetArray &, int *p=0) const override;
         int numDofs() const override;

        int             numNodes() const override;
        int * nodes(int *) const override;
        int getTopNumber() override;
        PrioInfo examine(int sub, MultiFront *) {
          fprintf(stderr,"SloshTetra.h: PrioInfo examine is commented in Dec.d/ElemMFCheck.C\n");
          return *(new PrioInfo);
        };

private:

        double volume(double dOmega) const {return dOmega/6.0;}
        double computeMetrics(const CoordSet &, double gN[4][3]) const;
        void buildTetraStiff(double m[4][4], double gN[4][3], double dOmega) const;
        
};

#endif
