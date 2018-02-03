#ifndef _CONNECTED_TRI_H_
#define _CONNECTED_TRI_H_

#include	<Element.d/Element.h>

class ConnectedTri : public Element {

	int nn[4];
public:
	ConnectedTri(int*);

	void renum(int *) override;
        void renum(EleRenumMap&) override;

        FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const;
        FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg=1) const;

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
         int numDofs() const override;

        int numNodes() const override;
        int * nodes(int *) const override;

        void             computeDisp(CoordSet&, State &, const InterpPoint &,
                                     double*, GeomState *gs);
	void		 getFlLoad(CoordSet &, const InterpPoint &,
                                   double *flF, double *resF, GeomState *gs=0);
	PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() override;

};

#endif
