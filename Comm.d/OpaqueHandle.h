//
// Created by Michel Lesoinne on 2/5/18.
//

#ifndef FEM_OPAQUEHANDLE_H
#define FEM_OPAQUEHANDLE_H

#include <type_traits>

/** \brief Opaque handle for any implementation specific MPI handle. */
class OpaqueHandle {
public:
	template <typename T>
	OpaqueHandle(T h) { new (&handle) T{h}; }
private:
	static const size_t ls = sizeof(int);
	static const size_t la = alignof(int);
	std::aligned_storage_t<ls, la> handle;
};

using CommunicatorHandle = OpaqueHandle;
//using

#endif //FEM_OPAQUEHANDLE_H
