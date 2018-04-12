//
// Created by Michel Lesoinne on 4/3/18.
//

#ifndef FEM_SUBIMPL_H
#define FEM_SUBIMPL_H

#include <vector>
#include <Driver.d/SComm.h>
#include "FetiLib/Types.h"

namespace FetiLib {

class SubImpl {
public:
	std::vector<global_node_index> glNodes;
	SComm scomm;
};

}


#endif //FEM_SUBIMPL_H
