#ifndef _ISOPARAMTRILINESOMMER_H_
#define _ISOPARAMTRILINESOMMER_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class IsoParamTriLineSommer : public SommerElement {

	int *nn;
	int order;
public:
	IsoParamTriLineSommer(int, int *, Element *_el = 0);

	int numNodes() const override { return order; }

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int numDofs() const { return order; }

	int numWetDofs() { return 3 * order; }

	int dim() const override { return 2; }

	int *dofs(DofSetArray &, int *p = 0);

	virtual IsoParamTriLineSommer *clone() override;

	int *wetDofs(DofSetArray &, int *p = 0);

	void neumVector(CoordSet &, ComplexVector &,
	                double, double, double, double, int pflag = 0);

	void wetInterfaceVector(CoordSet &, ComplexVector &,
	                        double, double, double, double, int, int);

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;

	GenStackFSFullMatrix<double> wetInterfaceMatrix(CoordSet &cs,
	                                                double *d);

	void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd);

	FullSquareMatrixC sommer2Matrix(CoordSet &, complex<double> *);

	FullSquareMatrix turkelMatrix(CoordSet &, double *) const override;

	ComplexD ffpCoef(double k) { return exp(ComplexD(0.0, M_PI / 4.0)) / sqrt(8.0 * M_PI * k); }

	void getNormal(const CoordSet &, double [3]) const override;

	void markDofs(DofSetArray &) override;
};

#endif
