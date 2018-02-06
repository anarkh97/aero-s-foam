//
// Created by Michel Lesoinne on 2/5/18.
//

#include <mpi.h>
#include <stdexcept>
#include "OpaqueCommunicator.h"

OpaqueCommunicator::OpaqueCommunicator(const OpaqueHandle &handle) : handle(handle) {}

class mpi_except : public std::runtime_error {
public:
	explicit mpi_except(int error_code) : std::runtime_error("MPI Error") {}
};

namespace com_details {
#ifdef USE_MPI

struct HandleHelper {
	static MPI_Comm mpi_comm(const OpaqueHandle & handle) {
		return *(reinterpret_cast<const MPI_Comm *>(&handle.handle));
	}

	static MPI_Datatype mpi_type(const OpaqueHandle & handle) {
		return *(reinterpret_cast<const MPI_Datatype *>(&handle.handle));
	}

	static MPI_Op mpi_op(const OpaqueHandle & handle) {
		return *(reinterpret_cast<const MPI_Op *>(&handle.handle));
	}
};

}

namespace {
MPI_Comm mpi_comm(const OpaqueHandle &handle) {
	return com_details::HandleHelper::mpi_comm(handle);
}

}

com_details::CommFunctions OpaqueCommunicator::functions{
		.commSize = [](const OpaqueHandle &handle) {
			int size;
			int result = MPI_Comm_size(com_details::HandleHelper::mpi_comm(handle), &size);
			if (result != MPI_SUCCESS)
				throw mpi_except(result);
			return size;
		},
		.remoteSize = [](const OpaqueHandle &handle) {
			int size;
			int result = MPI_Comm_remote_size(com_details::HandleHelper::mpi_comm(handle), &size);
			if (result != MPI_SUCCESS)
				throw mpi_except(result);
			return size;
		},
		.allReduce = [](const void *sendbuf, void *recvbuf, int count, TypeHandle datatype, OpHandle op,
		                const CommunicatorHandle &comm, void(*cpData)(const void *, void *, int)) {
			int result = MPI_Allreduce(sendbuf, recvbuf, count,
			                           com_details::HandleHelper::mpi_type(datatype), com_details::HandleHelper::mpi_op(op),
			                           com_details::HandleHelper::mpi_comm(comm));
			if (result != MPI_SUCCESS)
				throw mpi_except(result);
		},
		.barrier = [](const CommunicatorHandle &comm) {
			int result = MPI_Barrier(mpi_comm(comm));
			if (result != MPI_SUCCESS)
				throw mpi_except(result);
		},
		.rank = [](const OpaqueHandle &handle) {
			int rank;
			int result = MPI_Comm_rank(com_details::HandleHelper::mpi_comm(handle), &rank);
			if (result != MPI_SUCCESS)
				throw mpi_except(result);
			return rank;
		},
};

#else

com_details::CommFunctions OpaqueCommunicator::functions{
	.commSize = [](const OpaqueHandle &handle) {
		return 1;
	},
	.remoteSize = [](const OpaqueHandle &handle) {
		return 0;
	},
	.allReduce = [](const void *sendbuf, void *recvbuf, int count, TypeHandle datatype, OpHandle op,
		                CommunicatorHandle comm, void(*cpData)(const void *, void *, int)) {
		(*cpData)(sendbuf, recvbuf, count);
	},
	.barrier = [](const CommunicatorHandle &comm) {
	},
	.rank = [](const OpaqueHandle &handle) {
		return 0;
	}
};

template <typename T>
TypeHandle CommTypeTrait<T>::typeHandle() {
	return 0;
}

template struct CommTypeTrait<int>;
template struct CommTypeTrait<double>;
template struct CommTypeTrait<std::complex<double>>;

}
#endif
