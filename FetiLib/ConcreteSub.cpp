//
// Created by Michel Lesoinne on 2/23/18.
//

#include <Math.d/CuCSparse.h>
#include <Math.d/SparseSet.h>
#include "ConcreteSub.h"

namespace FetiLib {

//const FetiInfo &ConcreteBaseSub::getFetiInfo() const {
//	return <#initializer#>;
//}

void test() {
	ContreteSub<double> cs;
}

const FetiInfo &ConcreteBaseSub::getFetiInfo() const {
	return fetiInfo;
}
}