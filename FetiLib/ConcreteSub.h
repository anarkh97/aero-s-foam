//
// Created by Michel Lesoinne on 2/23/18.
//

#ifndef FEM_CONCRETESUB_H
#define FEM_CONCRETESUB_H

#include <Feti.d/FetiSub.h>
#include "Subdomain.h"

namespace FetiLib {

class ConcreteBaseSub : virtual public FetiBaseSub {
public:
	ConcreteBaseSub();
	const FetiInfo &getFetiInfo() const override;

	const int *getGlNodes() const override;

	int getNumUncon() const override;
	const CoordSet& getNodeSet() const override;
	double getShiftVal() const override;

	Connectivity *getNodeToNode() const override;

	DofSetArray *getDsa() const override; //!< Why do we need this?
	ConstrainedDSA *get_c_dsa() const override;
	void computeWaveNumbers() override;
	void averageMatProps() override;

private:
	FetiInfo fetiInfo;
	CoordSet coordinates;
};

template <typename Scalar>
class ConcreteSub : public ConcreteBaseSub, public FetiSub<Scalar> {
public:
	ConcreteSub(Subdomain subdomain);
};

}

#endif //FEM_CONCRETESUB_H
