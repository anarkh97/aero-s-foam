//
// Created by Michel Lesoinne on 2/5/18.
//

#ifndef FEM_OPAQUECOMMUNICATOR_H
#define FEM_OPAQUECOMMUNICATOR_H


#include <type_traits>
#include "OpaqueHandle.h"

class OpaqueCommunicator {
public:
	OpaqueCommunicator(const OpaqueHandle &handle);

private:
	OpaqueHandle handle;
};


#endif //FEM_OPAQUECOMMUNICATOR_H
