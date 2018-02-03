#ifndef _QUADRADIATION_H_
#define _QUADRADIATION_H_

#include <Element.d/Element.h>

class QuadRadiation: public virtual Element {

        int nn[4];
public:
        QuadRadiation(int*);
        ~QuadRadiation();

        Element *clone() override;

        void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg) override;
        FullSquareMatrix massMatrix(CoordSet& cs,double *d, int cmflg) override;
        double           getMass(CoordSet&);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

	Corotator* getCorotator(CoordSet &, double*, int, int);

        int numNodes() const override;
        int * nodes(int *) const override;
        int getTopNumber() override;
        PrioInfo examine(int sub, MultiFront *) override;

        bool isRadiationElement() { return true; }

        void computeTemp(CoordSet&, State &, double[2], double*); 
        void getFlFlux(double[2], double *, double *);
};

#endif
