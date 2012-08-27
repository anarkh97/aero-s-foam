#ifndef _PARAL_H_
#define _PARAL_H_

#if defined(sgi) &&  !defined(_OPENMP)
#include <ulocks.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

class ThreadLock {
#if defined(sgi) &&  !defined(_OPENMP)
      ulock_t lockV;
#endif
#if defined(_OPENMP)
     omp_lock_t lockV;
#endif
   public:
      ThreadLock();
      ~ThreadLock();
      void lock();
      void unlock();
};

/***************************************************************************
 * TaskDescr is a purely generic class definition. It should be subclassed
 * to specific task objects that cointer sufficient data to refer to the
 * work to be done, and overload the run() routine to perform the work
 ***************************************************************************/

class TaskDescr {
   public:
	virtual void run() {};
	virtual void runFor(int) {};
};

class ThreadManager;
class DistTimer;

class OneSproc {
   OneSproc *next;
   TaskDescr **allTasks;
   int numTasks;
   int step;
#if defined(sgi) &&  !defined(_OPENMP)
   usema_t *wait;
   usema_t *done;
#endif
   ThreadManager *tman;
   int myNum;

   static void run(void *);
   friend class ThreadManager;
};

/***************************************************************************
 * The ThreadManager class manages the creation, and deletion of a given
 * number of threads. Once it has been created, it can be given a list
 * of TaskDescr's that it will run on the available threads. In the default
 * case of a sequential machine, it will loop over the tasks to execute them
 *
 * NB: Because of this paradigm in which we do not know how many tasks are
 *     really running in parallel, a Task should not intend to communicate
 *     with another. It is only allowed to use Locks to insure proper
 *     update of shared data between tasks
 ***************************************************************************/
class ThreadManager {
	int numThreads;		// Number of threads
        int single;
	OneSproc *firstProc; 	// List of threads
        OneSproc *allProc;
#if defined(sgi) &&  !defined(_OPENMP)
        ulock_t sprocListLock;
        usema_t *readyProc;
        usema_t *allDone;
#endif
        DistTimer *timer;
   protected:
	static void threadStart(void *);
        void dispatch(TaskDescr *task);
   public:
	ThreadManager(int);	// Create a Manager with n threads
	~ThreadManager();
	void execTasks(int, TaskDescr **); // run the n given Tasks
        void execTasks(int, TaskDescr *);
        void execParal(int, TaskDescr **);
        void execParal(int, TaskDescr *);
        void execTimedParal(DistTimer &, int, TaskDescr **);
        void execTimedParal(DistTimer &, int, TaskDescr *);
        void memUsage();
        void memRequests();
        long getLocalMem();
	long memoryUsed();
        int numThr() { return numThreads; }
#if defined(sgi) &&  !defined(_OPENMP)
        barrier_t *getBarrier();
#endif
        friend class OneSproc;
};

extern ThreadManager *threadManager;

void paralDist( int nthr, void * (*f)(void *), void **);

#endif
