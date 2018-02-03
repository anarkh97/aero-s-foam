#ifndef _QUADSOMMERBC_H_
#define _QUADSOMMERBC_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class QuadSommerBC : public SommerElement {

	int nn[4];
public:
	QuadSommerBC(int, int, int, int, Element *_el = 0, int eType = 4);

	int numNodes() const override { return 4; }

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int numDofs() const override;

	int dim() const override { return 3; }

	int *dofs(DofSetArray &, int *p) override;

	QuadSommerBC *clone() override;

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;

	FullSquareMatrix turkelMatrix(CoordSet &, double *) const override;

	FullSquareMatrix refinedSommerMatrix(CoordSet &, double *) override;

	//FullSquareMatrix surfStiffMatrix(CoordSet&, double *);
	FullSquareMatrix HSommerMatrix(CoordSet &, double *) override;
	//FullSquareMatrix HKSommerMatrix(CoordSet&, double *);

	ComplexD ffpCoef(double k) override { return {0.25 / M_PI, 0.0}; }

	void getNormal(CoordSet &, double[3]) override;

	void markDofs(DofSetArray &) override;

	void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd) override;

private:
	void getLocalCoordinates(CoordSet &, double xx[4], double yy[4], double zz[4]) const;

	void SurfaceRefinement(int nNo, double *x, double *y, double *z, double *xx, double *yy, double *zz);

	void GaussCoordinates(int Ngp, double *Pg, double *weight);
};

#endif

