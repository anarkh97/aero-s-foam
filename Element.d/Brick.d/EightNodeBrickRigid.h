#ifndef _EIGHTNODEBRICKRIGID_H_
#define _EIGHTNODEBRICKRIGID_H_

#include <Element.d/Element.h>

class EightNodeBrickRigid: public Element {

	int nn[9];
public:
	EightNodeBrickRigid(int*);

	void renum(int *);

	FullSquareMatrix stiffness(CoordSet& cs, double *d, int flg=1);
        FullSquareMatrix massMatrix(CoordSet& cs, double *mel, int cmflg=1);


	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
  
        int              numInternalNodes() {return 1;}
        void             setInternalNodes(int *in) { nn[8] = in[0]; }

        int              getTopNumber();

	//PrioInfo examine(int sub, MultiFront *mf);
       bool hasRot() { return true; }

       double weight() { return 1.0; }
       double trueWeight() { return 1.0; }
       bool isRigidElement() { return true; }
};
#endif

