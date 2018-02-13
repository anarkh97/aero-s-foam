//
// Created by Michel Lesoinne on 2/5/18.
//

#ifndef FEM_OPAQUECOMMUNICATOR_H
#define FEM_OPAQUECOMMUNICATOR_H


#include <type_traits>
#include <algorithm>
#include "OpaqueHandle.h"

class BaseCommunicator;

struct RecDetails {
	int length;
	int source;
};

/** \brief Structure to obtain information about the underlying implementation of the standard. */
struct VersionInfo {
	int major;
	int minor;
};

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
	void (*blockingSend)(const void *sendbuf, int count, TypeHandle datatype,
	                     int dest, int tag,
	                     const CommunicatorHandle &);
	/** \brief Receive from any source
	 *
	 * @param buffer
	 * @param len
	 * @param datatype
	 * @param tag
	 * @return The information about the length received and the origin of the data.
	 */
	RecDetails (*blockingRec)(void *buffer, int len, TypeHandle datatype,
	                    int tag,
	                    const CommunicatorHandle &);
	void (*allGather)(const void *send_data, int send_count,
	                  void *recv_data, int recv_count, TypeHandle datatype,
	                  const CommunicatorHandle &);
	/** \brief Create a Window.
	 *
	 * @param data
	 * @param size
	 * @param disp_unit
	 * @return
	 */
	WinHandle (*createWindow)(void *data, int size, int disp_unit, const CommunicatorHandle &);

	void (*fence)(bool openOrClose, WinHandle handle);

	void (*destroyWindow)(WinHandle handle);

	void (*lock)(bool isShared, int remoteRank, WinHandle handle);
	void (*unlock)(int remoteRank, WinHandle handle);
	void (*flushRemote)(int remoteRank, WinHandle handle);
	void (*flushLocal)(int remoteRank, WinHandle handle);
	void (*lockAll)(WinHandle handle);
	void (*unlockAll)(WinHandle handle);
	void (*flushRemoteAll)(WinHandle handle);
	void (*flushLocalAll)(WinHandle handle);
	void (*fetchAndOp)(WinHandle handle, OpHandle op, const void *sourceData, TypeHandle datatype,
	                   void *resData, int remoteRank, int remoteOffset);
	void (*accumulate)(WinHandle handle, OpHandle op, const void *operand, int count, TypeHandle datatype,
	                   int remoteRank, int remoteOffset);
};

struct Constants {
	static WinHandle nullWindow;
};

}


class Window {
public:
	class LockGuard;
	Window(const Window &) = delete;
	Window(Window &&w);
	~Window();

	void open() const;
	void close() const;

	void sharedLock(int remoteRank) const;
	void unlock(int remoteRank) const;

	void flushRemote(int remoteRank) const;
	void flushLocal(int remoteRank) const;

	void sharedLockAll() const;
	void unlockAll() const;
	void flushRemote() const;
	void flushLocal() const;

	template <typename T>
	void fetchAndOp(OpHandle op, const T *operand, T*remoteOperandResult, int remoteRank, int remoteOffset) const;

	template <typename T>
	void accumulate(OpHandle op, const T *operand, int count, int remoteRank, int remoteOffset) const;
private:
	/// \brief Constructor of a window.
	Window(WinHandle winHandle) : winHandle(winHandle) {}
	friend class BaseCommunicator;
	WinHandle winHandle;
};

class Window::LockGuard {
	LockGuard(const Window &window, int remoteRank = -1) : window(window), remoteRank(remoteRank) {
		if(remoteRank >= 0)
			this->window.sharedLock(this->remoteRank);
		else
			this->window.sharedLockAll();
	}
	LockGuard(LockGuard &&);
	LockGuard(const LockGuard &) = delete;
	~LockGuard(){
		try {
			if(remoteRank >= 0)
				window.unlock(remoteRank);
			else
				window.unlockAll();
		} catch(...){}
	}
	void flushRemote() const {
		window.flushRemote(remoteRank);
	}
	void flushLocal() const {
		window.flushLocal();
	}
private:
	const Window &window;
	const int remoteRank; //!< Remote rank for a shared lock, -1 for an 'all' lock.
};

class BaseCommunicator {
public:
	template <typename T>
	BaseCommunicator(const T &t);

	BaseCommunicator(const CommunicatorHandle &handle);

	/** \brief Obtain the rank of the process in the communicator. */
	int rank() const { return (*functions.rank)(handle); }
	/** \brief Obtain the number of processes involved in the communicator. */
	int commSize() const { return (*functions.commSize)(handle); }
	/** \brief For heterogeneous communicator, obtain the number of process in the other side. */
	int remoteSize() const { return (*functions.remoteSize)(handle); }
	/** \brief Block until all the processes in the communicator have reached the barrier. */
	void barrier() const { (*functions.barrier)(handle); }

	Window window(void *d, int nBytes, int disp_unit) const;

	template <typename T>
	Window window(const T *data, int count) const {
		const void * dd = data;
		return window(const_cast<void *>(dd), count * sizeof(T), sizeof(T));
	}

	template <typename T>
	void allReduce(const T *sendbuf, T *recvbuf, int count,
	               OpHandle op) const {
		(*functions.allReduce)(sendbuf, recvbuf, count, CommTypeTrait<T>::typeHandle(), op, handle,
		[](const void *from, void *to, int count) {
			const T*f = reinterpret_cast<const T *>(from);
			T*t = reinterpret_cast<T *>(to);
			std::copy(f, f+count, t);
		});
	}

	template <typename T>
	void blockingSend(const void *sendbuf, int count,
	                  int dest, int tag) const {
		(*functions.blockingSend)(sendbuf, count, CommTypeTrait<T>::typeHandle(),
		                          dest, tag, handle);
	}

	template <typename T>
	RecDetails blockingRec(int tag, T *buffer, int len) const {
		return (*functions.blockingRec)(buffer, len, CommTypeTrait<T>::typeHandle(), tag, handle);
	}

	template <typename T>
	void allGather(const T *send_data, int send_count,
	               T *recv_data, int recv_count) {
		(*functions.allGather)(send_data, send_count, recv_data, recv_count,  CommTypeTrait<T>::typeHandle(), handle);
	}

	template <class Type>
	Type globalSum(Type) const;
	template <class Type>
	Type globalMax(Type) const;
	template <class Type>
	Type globalMin(Type) const;

private:
	CommunicatorHandle handle;

	friend class Window;
	static com_details::CommFunctions functions;
};

template <typename T>
T BaseCommunicator::globalSum(T t) const {
	T res;
	allReduce(&t, &res, 1, SumHandle);
	return res;
}

template <typename T>
T BaseCommunicator::globalMax(T t) const {
	T res;
	allReduce(&t, &res, 1, MaxHandle);
	return res;
}


template <typename T>
T BaseCommunicator::globalMin(T t) const {
	T res;
	allReduce(&t, &res, 1, MinHandle);
	return res;
}

inline
void Window::open() const {
	(*BaseCommunicator::functions.fence)(true, winHandle);
}

inline
void Window::close() const {
	(*BaseCommunicator::functions.fence)(false, winHandle);
}

inline Window::Window(Window &&w) : winHandle(w.winHandle) {
	w.winHandle = com_details::Constants::nullWindow;
}

inline Window::~Window() {
	// if(winHandle != com_details::Constants::nullWindow)
	(*BaseCommunicator::functions.destroyWindow)(winHandle);
}

inline
void Window::sharedLock(int remoteRank) const {
	(*BaseCommunicator::functions.lock)(true, remoteRank, winHandle);
}

inline
void Window::unlock(int remoteRank) const {
	(*BaseCommunicator::functions.unlock)(remoteRank, winHandle);
}

inline
void Window::flushRemote(int remoteRank) const {
	(*BaseCommunicator::functions.flushRemote)(remoteRank, winHandle);
}

inline
void Window::flushLocal(int remoteRank) const {
	(*BaseCommunicator::functions.flushLocal)(remoteRank, winHandle);
}

inline
void Window::sharedLockAll() const {
	(*BaseCommunicator::functions.lockAll)(winHandle);
}

inline
void Window::unlockAll() const {
	(*BaseCommunicator::functions.unlockAll)(winHandle);
}

inline
void Window::flushRemote() const {
	(*BaseCommunicator::functions.flushRemoteAll)(winHandle);
}

inline
void Window::flushLocal() const {
	(*BaseCommunicator::functions.flushLocalAll)(winHandle);
}

template<typename T>
void Window::fetchAndOp(OpHandle op, const T *operand, T *remoteOperandResult, int remoteRank, int remoteOffset) const {
	(*BaseCommunicator::functions.fetchAndOp)(winHandle, op, operand, CommTypeTrait<T>::typeHandle(),
	                                          remoteOperandResult, remoteRank, remoteOffset);
}

template<typename T>
void Window::accumulate(OpHandle op, const T *operand, int count,
                        int remoteRank, int remoteOffset) const {
	(*BaseCommunicator::functions.accumulate)(winHandle, op, operand, count, CommTypeTrait<T>::typeHandle(),
	                                          remoteRank, remoteOffset);
}

#endif //FEM_OPAQUECOMMUNICATOR_H
