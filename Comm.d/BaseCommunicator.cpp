//
// Created by Michel Lesoinne on 2/5/18.
//

#include <mpi.h>
#include <stdexcept>
#include "BaseCommunicator.h"
#include "MPICompatTraits.h"

BaseCommunicator::BaseCommunicator(const CommunicatorHandle &handle) : handle(handle) {}

class mpi_except : public std::runtime_error {
public:
	explicit mpi_except(int error_code) : std::runtime_error("MPI Error") {}
};


using WH = OpaqueTypedHandle<HandleType::Window>;
using CH = OpaqueTypedHandle<HandleType::Communicator>;

CH make_handle(MPI_Comm com) { return CH{com}; }

WH make_handle(MPI_Win win) { return WH(win); }

void testWin(MPI_Win win) {
}

MPI_Win testWH(const WH &wh) {
	testWin(wh);
	return wh;
}

#ifdef USE_MPI


namespace {

void _allGather(const void *send_data, int send_count,
                void *recv_data, int recv_count, TypeHandle datatype,
                const CommunicatorHandle &comm) {
	MPI_Allgather(send_data, send_count, datatype,
	              recv_data, recv_count, datatype, comm);
};

WinHandle _createWindow(void *data, int size, int disp_unit, const CommunicatorHandle &comm) {
	MPI_Win winHandle;
	int err_code = MPI_Win_create(data, size, disp_unit, MPI_INFO_NULL, comm, &winHandle);
	if (err_code != MPI_SUCCESS)
		throw mpi_except(err_code);
	return {winHandle};
};

void _fence(bool openOrClose, WinHandle handle) {
	int err_code = MPI_Win_fence(openOrClose ? MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE : MPI_MODE_NOSUCCEED,
	                             handle);
	if (err_code != MPI_SUCCESS)
		throw mpi_except(err_code);
}

void _destroyWindow(WinHandle handle) {
	MPI_Win mpiHandle = handle;
	if (mpiHandle != MPI_WIN_NULL) {
		int err_code = MPI_Win_free(&mpiHandle);
		if (err_code != MPI_SUCCESS)
			throw mpi_except(err_code);
	}

}

void check(int err_code) {
	if (err_code != MPI_SUCCESS)
		throw mpi_except(err_code);
}

void _lock(bool isShared, int remoteRank, WinHandle handle) {
	check(
			MPI_Win_lock(isShared ? MPI_LOCK_SHARED : MPI_LOCK_EXCLUSIVE,  remoteRank, 0, handle)
	);
}

void _unlock(int remoteRank, WinHandle handle) {
	check(
			MPI_Win_unlock(remoteRank, handle)
	);
}

void _lockAll(WinHandle handle) {
	check(
			MPI_Win_lock_all(0,  handle)
	);
}

void _unlockAll(WinHandle handle) {
	check(
			MPI_Win_unlock_all(handle)
	);
}

void _flushRemote(int remoteRank, WinHandle handle) {
	check(
			MPI_Win_flush(remoteRank, handle)
	);
}

void _flushLocal(int remoteRank, WinHandle handle) {
	check(
			MPI_Win_flush_local(remoteRank, handle)
	);
}

void _flushRemoteAll(WinHandle handle) {
	check(
			MPI_Win_flush_all(handle)
	);
}

void _flushLocalAll(WinHandle handle) {
	check(
			MPI_Win_flush_local_all(handle)
	);
}

void _fetchAndOp(WinHandle handle, OpHandle op, const void *sourceData, TypeHandle dataType,
                   void *resData, int remoteRank, int remoteOffset) {
	check(
			MPI_Fetch_and_op(sourceData, resData, dataType, remoteRank, remoteOffset, op, handle)
	);
}

void _accumulate(WinHandle handle, OpHandle op, const void *operand, int count, TypeHandle dataType,
                   int remoteRank, int remoteOffset) {
	check(
			MPI_Accumulate(operand, count, dataType, remoteRank, remoteOffset, count, dataType, op, handle)
	);
}

void _put(WinHandle handle, const void *sourceData, int count, TypeHandle datatype,
          int remoteRank, int remoteOffset) {
	check(
			MPI_Put(sourceData, count, datatype, remoteRank, remoteOffset, count, datatype, handle)
	);
}

void _get(WinHandle handle, void *destData, int count, TypeHandle datatype,
          int remoteRank, int remoteOffset) {
	check(
			MPI_Get(destData, count, datatype, remoteRank, remoteOffset, count, datatype, handle)
	);
}

}

com_details::CommFunctions getFunctionPointers() {
	com_details::CommFunctions result;
	result.commSize = [](const CommunicatorHandle &handle) {
		int size;
		int result = MPI_Comm_size(handle, &size);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
		return size;
	};
	result.remoteSize = [](const CommunicatorHandle &handle) {
		int size;
		int result = MPI_Comm_remote_size(handle, &size);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
		return size;
	};
	result.allReduce = [](const void *sendbuf, void *recvbuf, int count, TypeHandle datatype, OpHandle op,
	                      const CommunicatorHandle &comm, void(*cpData)(const void *, void *, int)) {
		int result = MPI_Allreduce(sendbuf, recvbuf, count,
		                           datatype,
		                           op,
		                           comm);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
	};
	result.barrier = [](const CommunicatorHandle &comm) {
		int result = MPI_Barrier(comm);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
	};
	result.rank = [](const CommunicatorHandle &handle) {
		int rank;
		int result = MPI_Comm_rank(handle, &rank);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
		return rank;
	};
	result.blockingSend = [](const void *sendbuf, int count, TypeHandle datatype, int dest, int tag,
	                         const CommunicatorHandle &comm) {
		int result = MPI_Send(sendbuf, count, datatype, dest, tag, comm);
		if (result != MPI_SUCCESS)
			throw mpi_except(result);
	};
	result.blockingRec = [](void *buffer, int len, TypeHandle datatype, int tag,
	                        const CommunicatorHandle &comm) {
		RecDetails rInfo;
		MPI_Status status;
		MPI_Recv(buffer, len,
		         datatype, MPI_ANY_SOURCE, tag, comm, &status);
		MPI_Get_count(&status, datatype, &rInfo.length);
		rInfo.source = status.MPI_SOURCE;
		return rInfo;
	};
	result.allGather = &_allGather;
	result.createWindow = &_createWindow;
	result.destroyWindow = &_destroyWindow;
	result.fence = &_fence;
	result.lock = &_lock;
	result.unlock = &_unlock;
	result.lockAll = &_lockAll;
	result.unlockAll = &_unlockAll;
	result.flushRemote = &_flushRemote;
	result.flushLocal = &_flushLocal;
	result.flushRemoteAll = &_flushRemoteAll;
	result.flushLocalAll = &_flushLocalAll;
	result.fetchAndOp = &_fetchAndOp;
	result.accumulate = &_accumulate;
	result.get = &_get;
	result.put = &_put;
	return result;
}

com_details::CommFunctions BaseCommunicator::functions = getFunctionPointers();

Window BaseCommunicator::window(void *d, int nBytes, int disp_unit) const {
	MPI_Win win;
	MPI_Win_create(d, nBytes, disp_unit, MPI_INFO_NULL, handle, &win);
	return Window(win);
}

template<>
BaseCommunicator::BaseCommunicator(const MPI_Comm &c) : handle(c) {}


CommunicatorHandle getWorldComm() {
	return CommunicatorHandle(MPI_COMM_WORLD);
}

namespace com_details {

WinHandle Constants::nullWindow{MPI_WIN_NULL};

}
#else

com_details::CommFunctions BaseCommunicator::functions{
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

namespace com_details {

template <typename T>
TypeHandle CommTypeTrait<T>::typeHandle() {
	return 0;
}

template struct CommTypeTrait<int>;
template struct CommTypeTrait<double>;
template struct CommTypeTrait<std::complex<double>>;

}
#endif
