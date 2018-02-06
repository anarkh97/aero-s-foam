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

	Element *clone() override;

	void renum(int *) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg) const;
	FullSquareMatrix acousticm(CoordSet&, double *kel);
	void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*);
	void getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD*,
	                    double kappa, double *waveDir);
	FullSquareMatrix massMatrix(const CoordSet&,double *mel,int cmflg=1) const;
	double getMass(const CoordSet& cs) const override;

	void markDofs(DofSetArray &) const override;
	int* dofs(DofSetArray &, int *p=0) const override;
	int numDofs() const override;

	int             numNodes() const override;
	int * nodes(int *) const override;

	void addFaces(PolygonSet *pset);
private:
	void computedxdxi(const CoordSet &cs, int nint, double (*derivatives)[10][3], Matrix33 *dxdxi, double *det) const;
	PrioInfo examine(int sub, MultiFront *) override;
	int getTopNumber() override {return (142);}
};
#endif
