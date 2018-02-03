#ifndef _PRESSUREELEMENT_H_
#define _PRESSUREELEMENT_H_

#include <Element.d/Sommerfeld.d/SommerElement.h>

class DofSet;

class GeomState;

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree>
class PressureElement : public SommerElement {
protected:
	int nNodes;               // number of nodes
	int *nn;                  // node numbers
	std::vector<BCond> terms;
	int nterms;
	PressureBCond *pbc;

	void addTerms(DofSet);

public:
	PressureElement(int *, PressureBCond *);

	~PressureElement();

	void renum(int *) override;

	void renum(EleRenumMap &);

	int numNodes() const override;

	int *nodes(int * = 0) const override;

	int getNode(int nd) const override { return nn[nd]; }

	const int *getNodes() const override { return nn; }
	int *getNodes() override { return nn; }

	int numDofs() const override;

	int *dofs(DofSetArray &, int * = 0) override;

	void markDofs(DofSetArray &) override;

	int findAndSetEle(CoordSet &cs, Elemset &eset, Connectivity *nodeToEle, int *eleTouch, int *eleCount, int myNum,
	                  int it) override;

	PressureBCond *getPressure() { return pbc; }

	void neumVector(CoordSet &, Vector &, int pflag = 0, GeomState * = 0, double t = 0);

	void neumVectorJacobian(CoordSet &, FullSquareMatrix &, int pflag = 0, GeomState * = 0, double t = 0);

	FullSquareMatrix sommerMatrix(CoordSet &, double *) const override;
};

#ifdef _TEMPLATE_FIX_

#include <Element.d/Sommerfeld.d/PressureElement.C>

#endif

#endif
