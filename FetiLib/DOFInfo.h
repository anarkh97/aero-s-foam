//
// Created by Michel Lesoinne on 2/1/18.
//

#ifndef FEM_DOFINFO_H
#define FEM_DOFINFO_H

#include <vector>

namespace FetiLib {

using gl_node_t = int;
using lc_node_t = int;

enum class DOFType {
	XDisp,
	YDisp,
	ZDisp,
	XRot,
	YRot,
	ZRot
};

/** \brief DOFInfo gives forthe local node number to which the DOF belongs and the type of DOF it is. */
using DOFInfo = std::vector<std::pair<lc_node_t, DOFType>>;

}

#endif //FEM_DOFINFO_H
