#ifndef _PHELPER_H_
#define _PHELPER_H_

#include <Timers.d/GetTime.h>
#include <Threads.d/Paral.h>
#include <Timers.d/DistTimer.h>


template <typename FType>
class FunctorExecuter : public TaskDescr {
public:
	explicit FunctorExecuter(const FType &ft) : fctor(ft) {}

	FunctorExecuter(const FunctorExecuter &) = default;

	void runFor(int) override;
private:
	const FType &fctor; //!< The functor to execute.
};

template <typename FType>
void
FunctorExecuter<FType>::runFor(int i) {
	fctor(i);
}

template <typename FType>
auto makeExecuter(const FType &ftor) {
	return FunctorExecuter<FType>{ftor};
}


template <typename TA, typename TB, typename ... FArgs, typename ... Args>
void execParal(int n, TA *target, void (TB::*f)(int, FArgs ...) const, Args &&...args)
{
	auto call =[&](int i) { (static_cast<const TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};

template <typename TA, typename TB, typename ... FArgs, typename ... Args>
void execParal(int n, TA *target, void (TB::*f)(int, FArgs ...), Args &&...args)
{
	auto call =[&](int i) { (static_cast<TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};

template <typename TA, typename TB, typename ... FArgs, typename ... Args>
void timedParal(DistTimer &timer, int n, TA *target, void (TB::*f)(int, FArgs ...), Args &&...args)
{
	double initTime = getTime();
	long initMem  = threadManager->memoryUsed();
	auto call =[&](int i) { (static_cast<TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	threadManager->execTimedParal(timer, n, &fe);
	timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
};


template <typename TA, typename TB, typename ... FArgs, typename ... Args, typename X = typename std::enable_if<sizeof...(FArgs) == sizeof...(Args)>::type>
void timedParal(DistTimer &timer, int n, TA *target, void (TB::*f)(int, FArgs ...) const, Args &&...args)
{
	auto call =[&](int i) { (static_cast<const TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	double initTime = getTime();
	long initMem  = threadManager->memoryUsed();
	threadManager->execTimedParal(timer, n, &fe);
	timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
};

namespace thread_details {

template<class A, class B>
class Temp {
};

template<class T>
class Temp<T, T> {
public:
	inline static T subEval(T z, int) { return z; }
};

template<class T>
class Temp<T, T *> {
public:
	inline static T subEval(T *z, int i) { return z[i]; }
};

template<class T>
class Temp<T &, T &> {
public:
	inline static T &subEval(T &z, int) { return z; }
};

template<class T>
class Temp<T &, T *> {
public:
	inline static T &subEval(T *z, int) { return *z; }
};

template<class T>
class Temp<T &, T> {
public:
	inline static T &subEval(T &z, int) { return z; }
};

}

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApply(int n, A **target, void(B::*fct)(Args ...),PassedArgs ...pargs) {
	auto call =[&](int i) { (
		static_cast<B *>(target[i])->*fct)(thread_details::Temp<Args, PassedArgs>::subEval(pargs, i)...);
	};
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApplyToAll(int n, A *target, void(B::*fct)(Args ...),PassedArgs ...pargs) {
	auto call =[&](int i) {
		(static_cast<B &>(target[i]).*fct)(thread_details::Temp<Args, PassedArgs>::subEval(pargs, i)...);
	};
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApplyToAll(int n, A **target, void(B::*fct)(Args ...),PassedArgs ...pargs) {
	paralApply(n, target, fct, pargs...);
};

#endif
