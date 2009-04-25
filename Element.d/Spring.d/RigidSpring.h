#ifndef _RIGIDSPRING_H_
#define _RIGIDSPRING_H_

#include        <Element.d/Element.h>

class RigidSpring : public Element {

        int nn[3];
public:

	RigidSpring(int*);

	void renum(int *);


	FullSquareMatrix stiffness(CoordSet&, double *kel,int flg=1);
	FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);

        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
 
        int              numInternalNodes() {return 1;}
        void             setInternalNodes(int *in) { nn[2] = in[0]; }

        int  getTopNumber() { return 101; }
        bool isRigidElement() { return true; }
};
#endif
