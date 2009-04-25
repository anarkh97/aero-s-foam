#ifndef _RIGIDSOLID_H_
#define _RIGIDSOLID_H_

#include <Element.d/Element.h>

class RigidSolid: public Element {

	int nnodes;
        int *nn;
public:
	RigidSolid(int numnode, int*);

	void renum(int *);

	FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
  
        int              numInternalNodes() {return nnodes-2;}
        void             setInternalNodes(int *in) { 
                  for (int i = 0; i < nnodes-2; ++i)
                      nn[i+nnodes] = in[i];
        }

        int getTopNumber() { return 101; }
        int numTopNodes() { return 2; }
        bool isRigidElement() { return true; }
};
#endif

