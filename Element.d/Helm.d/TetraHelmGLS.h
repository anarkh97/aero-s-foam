#ifndef _TETRAHELMGLS_H_
#define _TETRAHELMGLS_H_

#include <Element.d/Helm.d/HelmElement.h>

class TetraHelmGLS: public HelmElement, public Element {

	int nn[4];
	mutable double coef; // TODO Get rid of state variables!!!!!
public:
	TetraHelmGLS(int*);

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg=1) const;
	FullSquareMatrix acousticm(CoordSet&, double *kel);
	FullSquareMatrix massMatrix(const CoordSet&,double *mel, int cmflg=1) const;
	double           getMass(const CoordSet& cs) const;

	void             markDofs(DofSetArray &);
	int*             dofs(DofSetArray &, int *p=0);
	int numDofs() const override;

	int             numNodes() const override;
	int * nodes(int *) const override;
	int getTopNumber() override;
	PrioInfo examine(int sub, MultiFront *) override;
	int nDecFaces() const override { return 4;}
	int getDecFace(int iFace, int *fn);

	void            addFaces(PolygonSet *pset);

	virtual double helmCoef() { return coef; }

private:



	double TetraMass[4][4]; // Tetrahedra stiffness matrix
	double TetraStiff[4][4]; // Tetrahedra mass matrix


	//double		volume(CoordSet&);
	double volume(double dOmega) const {return dOmega/6.0;}
	double computeMetrics(const CoordSet &, double gN[4][3]) const;
	void buildTetraMass(double m[4][4], double dOmega) const;
	void buildTetraStiff(double s[4][4], double gN[4][3], double dOmega) const;

};
#endif

