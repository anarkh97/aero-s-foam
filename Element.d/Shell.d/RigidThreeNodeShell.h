#ifndef _RIGIDTHREENODESHELL_H_
#define _RIGIDTHREENODESHELL_H_

#include	<Element.d/Element.h>

class RigidThreeNodeShell : public Element {

	int nn[4];
public:
	RigidThreeNodeShell(int*);

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet&, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

        int              numInternalNodes() {return 1;}
        void             setInternalNodes(int *in) { nn[3] = in[0]; }


	int getTopNumber();
        
        bool isShell() { return true; }
        bool isRigidElement() { return true; }
};
#endif

