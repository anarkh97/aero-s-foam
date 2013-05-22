#ifndef _TETRAHELMGLS_H_
#define _TETRAHELMGLS_H_

#include <Element.d/Helm.d/HelmElement.h>

class TetraHelmGLS: public HelmElement, public Element {

	int nn[4];
        double coef;
public:
	TetraHelmGLS(int*);

	Element *clone();

	void renum(int *);
        void renum(EleRenumMap&);

        FullSquareMatrix stiffness(CoordSet&, double *kel, int flg=1);
        FullSquareMatrix acousticm(CoordSet&, double *kel);
        FullSquareMatrix massMatrix(CoordSet&,double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);

	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
	int getTopNumber();
	PrioInfo examine(int sub, MultiFront *);

	void            addFaces(PolygonSet *pset);

        virtual double helmCoef() { return coef; }

private:

	double J[3][3];
	double Jinv[3][3];
	double dOmega;
	double gN[4][3];

	double TetraMass[4][4]; // Tetrahedra stiffness matrix
	double TetraStiff[4][4]; // Tetrahedra mass matrix


	//double		volume(CoordSet&);
	double		volume() {return dOmega/6.0;} 
	void		computeMetrics(CoordSet&);
	void		buildTetraMass(double m[4][4]);
	void 		buildTetraStiff(double s[4][4]);
	
};
#endif

