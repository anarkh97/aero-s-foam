#include <cstdlib>
#include <cstdio>

#include <Driver.d/Domain.h>
#include <Timers.d/GetTime.h>
#include <Utils.d/linkfc.h>
#include "dec.h"

class ThreadManager;

#include <Threads.d/Paral.h>
#include <Utils.d/Memory.h>

extern bool estFlag;
extern bool weightOutFlag;
extern bool useFull;
extern int verboseFlag;

namespace Dec
{

// ThreadManager *threadManager;

int dec(int numProcessors,               // aka -p <number of processors>
	int numThreads,                  // aka -n <number of threads>
	int numSubdomains,               // aka -s <number of subdomains> 
	//string decompositionFileName,  // aka -d <decomposition file name>
	int topFlag                      // -1 --> do nothing
	                                 // 0   ~   aka -t
	                                 // 1   ~   aka -T
	                                 // 2   ~   aka -m
	)
{

 if(verboseFlag) filePrint(stderr, " ... Decomposing mesh with %d processors, %d threads, asking for %d subdomains ...\n",
           numProcessors,numThreads,numSubdomains);

 // threadManager = new ThreadManager(numThreads);

 double tTotal = getTime();

/*
 if(estFlag) {
   filePrint(stderr," ... Selected Option Accounts for Effects of Mid-Side Nodes ... \n");
   filePrint(stderr," Outputing Subdomain Memory Estimates in filename.memory ...\n");
 }
 else
   filePrint(stderr," ... Selected Option Ignores Effects of Mid-Side Nodes ...\n");

 if(weightOutFlag) 
   filePrint(stderr," ... Outputting Weight Distribution Statistics in filename.weight ...\n");
*/
 long decomposeMemory = -memoryUsed();

 geoSource->simpleDecomposition(numSubdomains, estFlag, weightOutFlag);

 decomposeMemory += memoryUsed();

 long totalMemory = /*readMemory + setUpDataMemory*/ + decomposeMemory;

 if(verboseFlag) filePrint(stderr," ... Total Elapsed Time and Memory for the decomposer is %14.5f sec and %14.3f Mb ...\n",
           (getTime() - tTotal)/1000.0, totalMemory/(1024.0*1024.0));

 // delete threadManager;
 return 0;
}

} // namespace Dec
