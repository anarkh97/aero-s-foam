#ifndef _HEVIBQUADGAL_H_
#define _HEVIBQUADGAL_H_

#include <Element.d/Element.h>
#include <cstdio>

class HEVibQuadGal: public Element {

	int nn[4];
public:
	HEVibQuadGal(int*);

	Element *clone() override;

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flag=1) const;
        FullSquareMatrix massMatrix(const CoordSet&,double *d, int cmflg=1) const;
        double           getMass(const CoordSet&) const;

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;
        int getTopNumber() override;

        bool             isHEVFluidElement()  { return true; }

	PrioInfo examine(int sub, MultiFront *)   {
           fprintf(stderr,"HEVibQuad.h: PrioInfo examine is commented in Dec.d/ElemMFCheck.C"); return *(new PrioInfo);
        };
/*
        void computeTemp(CoordSet&cs, State &state, double gp[2], double*res);
        void getFlFlux(double gp[2], double *flF, double *resF);
        void computeSloshDisp(Vector& heatflux, CoordSet &cs, Vector& elTemp, 
                               int hflInd);
        void getThermalForce(CoordSet &, Vector &, Vector &force, int, GeomState *geomState=0) { force.zero(); }
*/
};
#endif

