//
// Created by Michel Lesoinne on 2/1/18.
//

#include <Utils.d/Connectivity.h>
#include "Subdomain.h"
#include "DOFInfo.h"

Connectivity dofToNode(const FetiLib::DOFInfo &dofInfo) {
	using dof_t = int;
	using nt = decltype(dofInfo[0].first);
	std::vector<std::pair<size_t, nt>> dofNodeVec;
	dofNodeVec.reserve(dofInfo.size());
	for(dof_t i = 0; i < dofInfo.size(); ++i)
		dofNodeVec.emplace_back(i, dofInfo[i].first);
	return Connectivity::fromLinkRange(dofNodeVec);
}

template <typename S>
Connectivity dofToDof(const FetiLib::SparseMatrix<S> &matrix) {
	std::vector<int> pointers{matrix.outerIndexPtr(), matrix.outerIndexPtr()+matrix.outerSize()};
	std::vector<int> targets{matrix.innerIndexPtr(), matrix.innerIndexPtr()+matrix.innerSize()};
	return {static_cast<int>(pointers.size()-1), std::move(pointers), std::move(targets)};
}

namespace FetiLib {



}