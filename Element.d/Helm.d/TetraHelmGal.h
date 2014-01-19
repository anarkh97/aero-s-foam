#ifndef _TETRAHELMGAL_H_
#define _TETRAHELMGAL_H_

#include <Element.d/Helm.d/HelmElement.h>

typedef double Matrix33[3][3];

class TetraHelmGal: public HelmElement, public Element {

	int nn[4];
        static double tetra4_weights[4];
        static double tetra4_values[4][4];
        static double tetra4_derivatives[4][4][3];
public:
	TetraHelmGal(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&, double *kel, int flg);
        FullSquareMatrix acousticm(CoordSet&, double *kel);
        void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*,
                            double kappa, double *waveDir);
        FullSquareMatrix massMatrix(CoordSet&,double *mel,int cmflg=1);
	double           getMass(CoordSet& cs);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	int getTopNumber();
	PrioInfo examine(int sub, MultiFront *);
        int nDecFaces() { return 4;}
        int getDecFace(int iFace, int *fn);

private:

	void            addFaces(PolygonSet *pset);

        void computedxdxi(CoordSet &cs, int nint,
        double (*derivatives)[4][3], Matrix33 *dxdxi, double *det);
};
#endif

