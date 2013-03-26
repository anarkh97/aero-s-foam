#ifndef _CONNECTED_TRI_H_
#define _CONNECTED_TRI_H_

#include	<Element.d/Element.h>

class ConnectedTri : public Element {

	int nn[4];
public:
	ConnectedTri(int*);

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

        void             computeDisp(CoordSet&, State &, const InterpPoint &,
                                     double*, GeomState *gs);
	void		 getFlLoad(CoordSet &, const InterpPoint &,
                                   double *flF, double *resF, GeomState *gs=0);
	PrioInfo examine(int sub, MultiFront *);
	int getTopNumber();

};

#endif
