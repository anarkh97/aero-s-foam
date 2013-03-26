#ifndef _TETRA10HELMGAL_H_
#define _TETRA10HELMGAL_H_

#include <Element.d/Helm.d/HelmElement.h>

typedef double Matrix33[3][3];

class Tetra10HelmGal: public HelmElement, public Element {
        
	int nn[10];
        static double tetra10_weights[15];
        static double tetra10_values[15][10];
        static double tetra10_derivatives[15][10][3];
        static double tetra10_vertex_derivatives[10][10][3];

public:
	Tetra10HelmGal(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&, double *kel, int flg);
        FullSquareMatrix acousticm(CoordSet&, double *kel);
        void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*);
        void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*,
                            double kappa, double *waveDir);
        FullSquareMatrix massMatrix(CoordSet&,double *mel,int cmflg=1);
	double           getMass(CoordSet& cs);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);

	void            addFaces(PolygonSet *pset);
private:
        void            computedxdxi(CoordSet &cs, int nint, double (*derivatives)[10][3], Matrix33 *dxdxi, double *det);
	PrioInfo examine(int sub, MultiFront *);
	int getTopNumber() {return (142);}
};
#endif
