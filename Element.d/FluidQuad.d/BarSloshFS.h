#ifndef _BARSLOSHFS_H_
#define _BARSLOSHFS_H_

#include <Element.d/Element.h>

class BarSloshFS: public Element {

        int nn[2];
public:

	BarSloshFS(int*);

        Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet&,double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        double           getMass(CoordSet&);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
        PrioInfo examine(int sub, MultiFront *)   {
           fprintf(stderr,"BarSloshFS.h: PrioInfo examine is commented in Dec.d/ElemMFCheck.C."); return *(new PrioInfo);
        } ;
	int 		getTopNumber();
};
#endif
