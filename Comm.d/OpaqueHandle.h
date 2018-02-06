//
// Created by Michel Lesoinne on 2/5/18.
//

#ifndef FEM_OPAQUEHANDLE_H
#define FEM_OPAQUEHANDLE_H

#include <memory>
#include <type_traits>

namespace com_details {
struct HandleHelper;
}

/** \brief Opaque handle for any implementation specific MPI handle. */
class OpaqueHandle {
public:
	template <typename T>
	OpaqueHandle(const T h) {
		static_assert(sizeof(T) <= ls, "OpaqueHandle cannot store such a large object.");
		unsigned char *place = reinterpret_cast<unsigned char *>(&handle);
		T *res = new (place) T(h);
	}
private:
	static const size_t ls = sizeof(void *);
	static const size_t la = alignof(void *);
	std::aligned_storage_t<ls, la> handle;
	friend struct com_details::HandleHelper;
};

using CommunicatorHandle = OpaqueHandle;
using TypeHandle = OpaqueHandle;
using OpHandle = OpaqueHandle;

extern OpHandle MaxHandle;
extern OpHandle MinHandle;
extern OpHandle SumHandle;
extern TypeHandle IntHandle;

template <typename T>
struct CommTypeTrait {
	static TypeHandle typeHandle();
};

#endif //FEM_OPAQUEHANDLE_H
