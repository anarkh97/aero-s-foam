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
	ConcreteSub<double> cs;
}

const FetiInfo &ConcreteBaseSub::getFetiInfo() const {
	return fetiInfo;
}

ConcreteBaseSub::ConcreteBaseSub() {
	subNumber =
	localSubNumber = 0;

}

void ConcreteBaseSub::computeWaveNumbers() {

}

void ConcreteBaseSub::averageMatProps() {

}

const CoordSet &ConcreteBaseSub::getNodeSet() const {
	return coordinates;
}


}