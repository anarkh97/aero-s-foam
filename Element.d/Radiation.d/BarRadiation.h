#ifndef _BARRADIATION_H_
#define _BARRADIATION_H_

#include <Element.d/Element.h>

class BarRadiation: public virtual Element {

        int nn[2];
public:
	BarRadiation(int*);
        ~BarRadiation();

        Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(const CoordSet&,double *kel, int flg=1) const;
        FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg=1) const;
        double           getMass(const CoordSet&) const;

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        Corotator* getCorotator(CoordSet &, double*, int, int);
        
        int numNodes() const override;
        int * nodes(int *) const override;
	PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() override;

        bool isRadiationElement() { return true; }

        void computeTemp(CoordSet&, State &, double[2], double*);
        void getFlFlux(double[2], double *, double *);
};

#endif
