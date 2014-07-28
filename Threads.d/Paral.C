#include <cstdlib>
#include <algorithm>
#include <sys/types.h>

#if defined(sgi) &&  !defined(_OPENMP)
#include <sys/prctl.h>
#include <ulocks.h>
#endif

#include <unistd.h>

//--- UH ---
#ifndef NO_MALLOC_DOT_H
#include <malloc.h>
#endif
//--- UH ---

#ifndef WINDOWS
#include <sys/mman.h>
#endif
#include <fcntl.h>

#include <Threads.d/Paral.h>
#include <Utils.d/DistHelper.h>
#include <Timers.d/DistTimer.h>
#include <Timers.d/GetTime.h>

#if defined(sgi) &&  !defined(_OPENMP)
usptr_t * usPtr;

// isParal = 0 sequential
// isParal = 1 parallel

// NOTE: if you change this variable in a 
//       sub-directory, you can force 
//       sequential execution for that part of the code.

int isParal = 1;

int initSize=2*1024*1024;
//int growth = 16*1024*1024;
int growth = 32*1024*1024; 
//int growth =  64*1024*1024;

ulock_t allocLock = 0;

pid_t pproc = 0;
int numProc = 0;
pid_t *childProc = 0;
void **arenas;
int *maxSizes;
long *curSizes;
void *pparena = 0;
int zeroFd = -1;
#endif
long currentSizes = 0;


#if defined(sgi) &&  !defined(_OPENMP)
void *
arenaGrow(size_t size, void *)
{
 // sbrk must be sequentialized
/*double t0 = -getTime();
 ussetlock(allocLock);
 t0 += getTime();
 double t1 = -getTime();
 void *p = sbrk(size);
 t1 += getTime();
 usunsetlock(allocLock);
 fprintf(stderr,"Growth timings: %f %f\n",t0,t1);
 return p;
 return malloc(size);
*/

 void *gr =  mmap(0, size, PROT_READ|PROT_WRITE , MAP_PRIVATE, zeroFd, 0);
 if(gr == 0) { fprintf(stderr,"Out of memory\n"); exit(-1); }

// Test on touching initial memory allocation
/*
 char *mem = (char *)gr;
 int i;
 for(i=0; i < size; i += 4096)
   mem[i] =  0;
*/
 return gr;
}
#endif

bool loud = false;

void * operator new(size_t size)
{
#if defined(sgi) && !defined(_OPENMP) && defined(NEWNEW)
// fprintf(stderr,"allocated %20d bytes  Running Total %14.3f Mb \n",
//                 size, memoryUsed()/(1024.0*1024.0));

 // get my process ID
 pid_t myPid = getpid();

 if(pproc == 0 || myPid == pproc) {
   if(pparena == 0) {
      pparena = malloc(initSize);
      zeroFd = open("/dev/zero",O_RDWR);
      acreate(pparena,  initSize, 0, 0, &arenaGrow);
      amallopt(M_BLKSZ, growth, pparena);
   }
   currentSizes += size;
   return amalloc(size,pparena);
 }

 // Find my process number
 int iProc; 
 for(iProc = 0; iProc < numProc && myPid != childProc[iProc]; ++iProc)
  ;

 if(iProc == numProc) {
   fprintf(stderr,"Error\n");
   return malloc(size);
 }

 if(loud)
   fprintf(stderr,"MA %d\n", iProc);
 curSizes[iProc] += size;
 void *p = amalloc(size,arenas[iProc]); 
 return p;
#else
  //currentSizes += size;
  currentSizes = (long(size)/long(1024)) + currentSizes; // PJSA convert to KB since currentSizes is long (max 2147483647 = 2 GB, not enough)
  return malloc(size);
#endif
}

long
ThreadManager::getLocalMem()
{
#if defined(sgi) && !defined(_OPENMP) && defined(NEWNEW)
 // get my process ID
 pid_t myPid = getpid();
 if(pproc == 0 || myPid == pproc)
   return currentSizes;

 // Find my process number
 int iProc;
 for(iProc = 0; iProc < numProc && myPid != childProc[iProc]; ++iProc)
  ;
 if(iProc == numProc) {
   fprintf(stderr,"Error\n");
   return 0;
 }
 return curSizes[iProc];
#else
 return currentSizes;
#endif
}

void operator delete(void *p)
{
#if defined(sgi) && !defined(_OPENMP) && defined(NEWNEW)
 pid_t myPid = getpid();
 if(myPid == pproc || pproc==0) {
    currentSizes -= amallocblksize(p,pparena);
    afree(p,pparena);
    return;
 }
 myPid = getpid();
 // Find my process number
 int iProc; 
 for(iProc = 0; iProc < numProc && myPid != childProc[iProc]; ++iProc)
  ;
 if(iProc == numProc) {
   fprintf(stderr,"Error\n");
   free(p);
   return;
 }
 // Free using my arena...
 curSizes[iProc] -= amallocblksize(p,arenas[iProc]);
 afree(p,arenas[iProc]);
#else
 free(p);
#endif
}

void
ThreadLock::lock()
{
#if defined(sgi) &&  !defined(_OPENMP)
 ussetlock(lockV);
#endif
#if defined(_OPENMP)
  omp_set_lock(&lockV);
#endif
}

void
ThreadLock::unlock()
{
#if defined(sgi) &&  !defined(_OPENMP)
 usunsetlock(lockV);
#endif
#if defined(_OPENMP)
  omp_unset_lock(&lockV);
#endif

}


ThreadLock::ThreadLock()
{
#if defined(sgi) &&  !defined(_OPENMP)
 lockV = usnewlock(usPtr);
#endif
#if defined(_OPENMP)
  omp_init_lock(&lockV);
#endif
}

ThreadLock::~ThreadLock()
{
#if defined(sgi) &&  !defined(_OPENMP)
#endif
#if defined(_OPENMP)
  omp_destroy_lock(&lockV);
#endif
}

typedef void (*P)(void *, size_t);
ThreadManager::ThreadManager(int nThr)
{
#if defined(sgi) &&  !defined(_OPENMP)
 numThreads = nThr;
 // fprintf(stderr, "Creating %d threads\n", numThreads);
 
 usconfig(CONF_INITUSERS, numThreads);
 usPtr = usinit("/dev/zero");

 firstProc = 0;
 allDone = usnewsema(usPtr , 0);
 readyProc =  usnewsema(usPtr , 0);
 sprocListLock = usnewlock(usPtr);
 allocLock = usnewlock(usPtr);

 allProc = new OneSproc[numThreads-1];

 // Setup arenas
 numProc = numThreads-1;
 arenas = new void *[numThreads-1];
 curSizes = new long[numThreads-1];
 maxSizes = new int[numThreads-1];
 childProc = new pid_t[numThreads-1];
 int i;
 for(i = 0; i < numThreads-1; ++i) {
    OneSproc *thisProc = allProc + i;
    thisProc->done = allDone;
    thisProc->wait = usnewsema(usPtr , 0);
    if(thisProc->wait == 0) perror("Could not get a semaphore:");
    thisProc->tman = this;
    thisProc->myNum = i;
    thisProc->step = numThreads;
//    pid_t pid = sproc(OneSproc::run, PR_SALL, thisProc);
    // fprintf(stderr, "Calling sprocsp\n");
    pid_t pid = sprocsp((P)OneSproc::run, PR_SALL, thisProc,0,10000000);
    if(pid < 0) perror("Thread did not start:");

    childProc[i] = pid;
    curSizes[i] = maxSizes[i] = 0;
 }
 for(i = 0; i < numThreads-1; ++i) {
//  set up a memory arena;
    arenas[i] = malloc(initSize);
    acreate(arenas[i], initSize, 0, 0, &arenaGrow);
    amallopt(M_BLKSZ, growth, arenas[i]);
    curSizes[i] = maxSizes[i] = 0;
 }
 pproc = getpid();
 //filePrint(stderr," ... Setting Memory Arenas          ...\n");
#elif defined(_OPENMP)
 numThreads = nThr;
// fprintf(stderr, "Forcing %d threads\n",numThreads);
 omp_set_dynamic(0);
 omp_set_num_threads(numThreads);
#else
 nThr = 1;
 numThreads = nThr;
 firstProc=0;
#endif
 //fprintf(stderr, " In ThreadManager::ThreadManager(...), Creating %d threads\n", numThreads);
}

ThreadManager::~ThreadManager()
{
#if defined(sgi) &&  !defined(_OPENMP)
 int i;
 for(i = 0; i < numThreads-1; ++i)
   allProc[i].allTasks = 0;
 for(i = 0; i < numThreads-1; ++i)
   usvsema(allProc[i].wait);
 //long totSizes = 0;
 //for(i = 0; i < numThreads-1; ++i)
 // totSizes += curSizes[i];
 // fprintf(stderr,
 //  "Total memory consumption in megabytes = %10.3f\n",totSizes/(1024.*1024.));
#endif
}

long
ThreadManager::memoryUsed()
{
//HB: decommented this
#if defined(sgi) &&  !defined(_OPENMP)
 int i;
 long totSizes = currentSizes;
 for(i = 0; i < numThreads-1; ++i)
    totSizes += curSizes[i];
    return totSizes;
#else
 return currentSizes;
#endif
}

void
ThreadManager::memUsage()
{
#if defined(sgi) &&  !defined(_OPENMP)
 long totSizes = 0;

 long minMem = currentSizes;
 // long avgMem = 0;
 long maxMem = currentSizes;

 int i;
 for(i = 0; i < numThreads-1; ++i) {
    minMem = std::min( minMem, curSizes[i]);
    maxMem = std::max( maxMem, curSizes[i]);
    totSizes += curSizes[i];
 }

 minMem = std::min( minMem, currentSizes );
 maxMem = std::max( maxMem, currentSizes );

 totSizes += currentSizes;

 // avgMem = totSizes/numThreads;
#endif
}

#if defined(sgi) &&  !defined(_OPENMP)
void
OneSproc::run(void *p)
{
 OneSproc *sp = (OneSproc *) p;
 ThreadManager *tman = sp->tman;

 while(1) {
     uspsema(sp->wait);
     if(sp->allTasks == 0) return;
     long initMem;
     double initTime;
     if(tman->timer) {
       initTime = getTime();
       initMem  = threadManager->getLocalMem();
     }
     int index;
     for(index = sp->myNum; index < sp->numTasks; index+= sp->step) {
       if(tman->single) sp->allTasks[0]->runFor(index);
       else sp->allTasks[index]->run();
     }
     if(tman->timer) {
       long finalMem = threadManager->getLocalMem();
       tman->timer->addTo(sp->myNum, finalMem-initMem, getTime()-initTime);
     }
     usvsema(sp->done);
 }
}
#endif


void
ThreadManager::execTasks(int ntasks, TaskDescr **td)
{
 int i;
 for(i = 0; i < ntasks; ++i) 
   td[i]->run();
}

void
ThreadManager::execTasks(int ntasks, TaskDescr *td)
{
 int i;
 for(i = 0; i < ntasks; ++i) {
   td->runFor(i);
 }
}

void
ThreadManager::execParal(int ntasks, TaskDescr **td)
{
  timer = 0;
#if defined(sgi) &&  !defined(_OPENMP)
  int i;
  if(isParal == 0)
    for(i=0; i < ntasks; ++i)
      td[i]->run();
  else {
    for(i = 0; i < numThreads-1; ++i) {
       allProc[i].allTasks = td;
       allProc[i].numTasks = ntasks;
    }
    single = 0;
    for(i=0; i < numThreads-1; ++i)
      usvsema(allProc[i].wait);
 
    int index;
    for(index = numThreads-1; index < ntasks; index+= numThreads)
      td[index]->run();
 
    for(i=0; i < numThreads-1; ++i)
      uspsema(allDone);
  }
#else
 int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
 for(i=0; i < ntasks; ++i)
    td[i]->run();
#endif
}

void
ThreadManager::execParal(int ntasks, TaskDescr *td)
{
  timer = 0;
#if defined(sgi) &&  !defined(_OPENMP)
  int i;
  if(isParal == 0)
    for(i = 0; i < ntasks; ++i) {
      td->runFor(i);
    }
  else {
    for(i = 0; i < numThreads-1; ++i) {
      allProc[i].allTasks = &td;
      allProc[i].numTasks = ntasks;
    }
    single = 1;
    for(i=0; i < numThreads-1; ++i)
      usvsema(allProc[i].wait);

    int index;
    for(index = numThreads-1; index < ntasks; index+= numThreads)
      td->runFor(index);

    for(i=0; i < numThreads-1; ++i)
      uspsema(allDone);
  }
#else
 int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
 for(i = 0; i < ntasks; ++i)
   td->runFor(i);
#endif
}


void
ThreadManager::execTimedParal(DistTimer &thisTimer, int ntasks, TaskDescr **td)
{

    timer = &thisTimer;
#if defined(sgi) &&  !defined(_OPENMP)
    int i;
    for(i = 0; i < numThreads-1; ++i) {
      allProc[i].allTasks = td;
      allProc[i].numTasks = ntasks;
    }
    single = 0;
    for(i=0; i < numThreads-1; ++i)
      usvsema(allProc[i].wait);

    double initTime = getTime();
    long initMem  = threadManager->getLocalMem();

    int index;
    for(index = numThreads-1; index < ntasks; index+= numThreads)
      td[index]->run();

    long finalMem = threadManager->getLocalMem();
    timer->addTo(numThreads-1, finalMem-initMem, getTime()-initTime);

    for(i=0; i < numThreads-1; ++i)
      uspsema(allDone);
#else
 double initTime = getTime();
 long initMem  = threadManager->getLocalMem();
 int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
 for(i=0; i < ntasks; ++i)
    td[i]->run();

 long finalMem = threadManager->getLocalMem();
 timer->addTo(0, finalMem-initMem, getTime()-initTime);
#endif
}

void
ThreadManager::execTimedParal(DistTimer &thisTimer, int ntasks, TaskDescr *td)
{
    timer = &thisTimer;
#if defined(sgi) &&  !defined(_OPENMP)
    int i;
    for(i = 0; i < numThreads-1; ++i) {
       allProc[i].allTasks = &td;
       allProc[i].numTasks = ntasks;
    }
    single = 1;
    for(i=0; i < numThreads-1; ++i)
      usvsema(allProc[i].wait);

    double initTime = getTime();
    long initMem  = threadManager->getLocalMem();

    int index;
    for(index = numThreads-1; index < ntasks; index+= numThreads)
      td->runFor(index);

    long finalMem = threadManager->getLocalMem();
    timer->addTo(numThreads-1, finalMem-initMem, getTime()-initTime);

    for(i=0; i < numThreads-1; ++i)
      uspsema(allDone);
#else
 double initTime   = getTime();
 long initMem = threadManager->getLocalMem();

 int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1)
#endif
 for(i = 0; i < ntasks; ++i)
   td->runFor(i);

 long finalMem = threadManager->getLocalMem();
 timer->addTo(0, finalMem-initMem, getTime()-initTime);
#endif
}


#if defined(sgi) &&  !defined(_OPENMP)
barrier_t *
ThreadManager::getBarrier()
{
 return new_barrier(usPtr);
}
#endif

class SimpleExec : public TaskDescr {
   void **dt;
   void *(*f)(void *);
 public:
   SimpleExec(void *(*_f)(void *), void **d) : dt(d), f(_f) {}
   void runFor(int);
};

void
SimpleExec::runFor(int i)
{
 (*f)(dt[i]);
}

void
paralDist(int nt, void *(*_f)(void *), void **d)
{
  filePrint(stderr, "Doing Spooles in Parallel\n");
  SimpleExec se(_f, d);
  threadManager->execParal(nt, &se);
}
