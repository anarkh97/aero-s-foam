//
// Created by Michel Lesoinne on 2/5/18.
//

#ifndef FEM_OPAQUECOMMUNICATOR_H
#define FEM_OPAQUECOMMUNICATOR_H


#include <type_traits>
#include <algorithm>
#include "OpaqueHandle.h"

class OpaqueCommunicator;

namespace com_details {

struct CommFunctions {
	int (*commSize)(const CommunicatorHandle &handle);
	int (*remoteSize)(const CommunicatorHandle &handle);
	void (*allReduce)(const void *sendbuf, void *recvbuf, int count,
	                  TypeHandle datatype, OpHandle op,
	                  const CommunicatorHandle &comm,
	                  void(*cpData)(const void *, void *, int));
	void (*barrier)(const CommunicatorHandle &);
	int (*rank)(const CommunicatorHandle &);
};

}


class OpaqueCommunicator {
public:
	OpaqueCommunicator(const OpaqueHandle &handle);

	int rank() const { return (*functions.rank)(handle); }
	int comSize() const { return (*functions.commSize)(handle); }
	int remoteSize() const { return (*functions.remoteSize)(handle); }

	void barrier() const { (*functions.barrier)(handle); }

	template <typename T>
	void allReduce(const T *sendbuf, T *recvbuf, int count,
	               TypeHandle datatype, OpHandle op) const {
		(*functions.allReduce)(sendbuf, recvbuf, count, datatype, op, handle,
		[](const void *from, void *to, int count) {
			const T*f = reinterpret_cast<const T *>(from);
			T*t = reinterpret_cast<T *>(to);
			std::copy(f, f+count, t);
		});
	}

	operator const OpaqueHandle &() const { return handle; };
private:
	OpaqueHandle handle;

	static com_details::CommFunctions functions;
};


#endif //FEM_OPAQUECOMMUNICATOR_H
