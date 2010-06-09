#ifndef _PITATIMERS_H_
#define _PITATIMERS_H_

#include <cstdio>
using std::FILE;

class ControlInfo;

class PitaTimers {

  public:
    
    int numCPU, myCPU;
 
    double buildTStime, dynOpstime, step0time, Bitime, PitaTotaltime, solveITA, exchangeProptime, 
           CkTotaltime, CkfineProjtime, CkcoarseProjtime, basetime, newCondtime, getInitCondtime,
           TStime, StoreFineGridTime, exchange_info, compNorm, newTSid, jump;
 
    FILE *fptr;
    
    PitaTimers() {
          exchange_info = 0.0;
          compNorm = 0.0;
          jump     = 0.0;
          newTSid  = 0.0;

          PitaTotaltime = 0.0;
          buildTStime   = 0.0;
          dynOpstime    = 0.0;
          step0time     = 0.0;
          Bitime        = 0.0;
          solveITA      = 0.0;

          exchangeProptime = 0.0;
          CkTotaltime      = 0.0;
          CkfineProjtime   = 0.0;
          CkcoarseProjtime = 0.0;
          basetime         = 0.0;
          newCondtime      = 0.0;
      
          getInitCondtime  = 0.0;
          
          TStime     = 0.0;
   
          StoreFineGridTime = 0.0;

          fptr = NULL;
   }
    
   ~PitaTimers() {
      if (fptr) fclose(fptr); 
   } 
   
   void openPitaTimersFile(ControlInfo &cinfo);
   void printPitaTimesFile();
};

#endif
