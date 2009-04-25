#include <Pita.d/PitaTimers.h>
#include <Driver.d/StructProp.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/Communicator.h>
#include <Comm.d/Communicator.h>
#include <Driver.d/GeoSource.h>
#include <ctime>

using namespace std;

extern Communicator *structCom;
extern GeoSource *geoSource;

//------------------------------------------------------------------
void PitaTimers::openPitaTimersFile(ControlInfo &cinfo)
{
  if (fptr == 0) {
     //snprintf is used to avoid problems of size encountered with sprintf
     int len = strlen(cinfo.checkfile)+strlen(".pitaTiming")+1;
     char *Name = new char[len];
     snprintf(Name,len,"%s%s",cinfo.checkfile,".pitaTiming");
     
     if ( (fptr=fopen(Name,"w"))==(FILE *) 0 )
       filePrint(stderr," *** ERROR: Cannot open %s.pitaTiming ...\n", cinfo.checkfile);
     fflush(fptr);
  }

}

//------------------------------------------------------------------
void PitaTimers::printPitaTimesFile()
{
#ifdef DISTRIBUTED

   numCPU = structCom->numCPUs();
   myCPU  = structCom->myID();
   
   if (myCPU==0) {
     ControlInfo *cinfo = geoSource->getCheckFileInfo();
     openPitaTimersFile(cinfo[0]);
   }
   
   double PitaTotaltime_Min, PitaTotaltime_Max, PitaTotaltime_Avg;   
   double buildTStime_Min, buildTStime_Max, buildTStime_Avg;
   double dynOpstime_Min, dynOpstime_Max, dynOpstime_Avg;
   double step0time_Min, step0time_Max, step0time_Avg;
   double Bitime_Min, Bitime_Max, Bitime_Avg;
   double solveITA_Min, solveITA_Max, solveITA_Avg;
   double exchangeProptime_Min, exchangeProptime_Max, exchangeProptime_Avg;
   double CkTotaltime_Min, CkTotaltime_Max, CkTotaltime_Avg;
   double CkfineProjtime_Min, CkfineProjtime_Max, CkfineProjtime_Avg;
   double CkcoarseProjtime_Min, CkcoarseProjtime_Max, CkcoarseProjtime_Avg; 
   double basetime_Min, basetime_Max, basetime_Avg;
   double newCondtime_Min, newCondtime_Max, newCondtime_Avg;
   
   double getInitCondtime_Min, getInitCondtime_Max, getInitCondtime_Avg;

   double TStime_Max;

   double StoreFineGridTime_Min, StoreFineGridTime_Max, StoreFineGridTime_Avg;

   double exchange_info_Max, compNorm_Max, newTSid_Max, jump_Max;

   structCom->sync();

   PitaTotaltime_Min = structCom->globalMin(PitaTotaltime);
   PitaTotaltime_Max = structCom->globalMax(PitaTotaltime);
   PitaTotaltime_Avg = structCom->globalSum(PitaTotaltime)/numCPU;

   buildTStime_Min = structCom->globalMin(buildTStime);
   buildTStime_Max = structCom->globalMax(buildTStime);
   buildTStime_Avg = structCom->globalSum(buildTStime)/numCPU;
 
   dynOpstime_Min = structCom->globalMin(dynOpstime);
   dynOpstime_Max = structCom->globalMax(dynOpstime);
   dynOpstime_Avg = structCom->globalSum(dynOpstime)/numCPU;
  
   step0time_Min = structCom->globalMin(step0time);
   step0time_Max = structCom->globalMax(step0time);
   step0time_Avg = structCom->globalSum(step0time)/numCPU;

   Bitime_Min = structCom->globalMin(Bitime);
   Bitime_Max = structCom->globalMax(Bitime);
   Bitime_Avg = structCom->globalSum(Bitime)/numCPU;
   
   solveITA_Min = structCom->globalMin(solveITA);
   solveITA_Max = structCom->globalMax(solveITA);
   solveITA_Avg = structCom->globalSum(solveITA)/numCPU;

   exchangeProptime_Min = structCom->globalMin(exchangeProptime);
   exchangeProptime_Max = structCom->globalMax(exchangeProptime);
   exchangeProptime_Avg = structCom->globalSum(exchangeProptime)/numCPU;

   CkTotaltime_Min = structCom->globalMin(CkTotaltime);
   CkTotaltime_Max = structCom->globalMax(CkTotaltime);
   CkTotaltime_Avg = structCom->globalSum(CkTotaltime)/numCPU;

   CkfineProjtime_Min = structCom->globalMin(CkfineProjtime);
   CkfineProjtime_Max = structCom->globalMax(CkfineProjtime);
   CkfineProjtime_Avg = structCom->globalSum(CkfineProjtime)/numCPU;

   CkcoarseProjtime_Min = structCom->globalMin(CkcoarseProjtime);
   CkcoarseProjtime_Max = structCom->globalMax(CkcoarseProjtime);
   CkcoarseProjtime_Avg = structCom->globalSum(CkcoarseProjtime)/numCPU;

   basetime_Min = structCom->globalMin(basetime);
   basetime_Max = structCom->globalMax(basetime);
   basetime_Avg = structCom->globalSum(basetime)/numCPU;

   newCondtime_Min = structCom->globalMin(newCondtime);
   newCondtime_Max = structCom->globalMax(newCondtime);
   newCondtime_Avg = structCom->globalSum(newCondtime)/numCPU;

   getInitCondtime_Min = structCom->globalMin(getInitCondtime);
   getInitCondtime_Max = structCom->globalMax(getInitCondtime);
   getInitCondtime_Avg = structCom->globalSum(getInitCondtime)/numCPU;
   
   TStime_Max = structCom->globalMax(TStime);
   
   StoreFineGridTime_Min = structCom->globalMin(StoreFineGridTime);
   StoreFineGridTime_Max = structCom->globalMax(StoreFineGridTime);
   StoreFineGridTime_Avg = structCom->globalSum(StoreFineGridTime)/numCPU;

   exchange_info_Max = structCom->globalMax(exchange_info);
   compNorm_Max = structCom->globalMax(compNorm);
   newTSid_Max  = structCom->globalMax(newTSid);
   jump_Max = structCom->globalMax(jump);

   structCom->sync();
   
   if (myCPU==0) {

   filePrint(fptr,"*********************************************************************\n");
   filePrint(fptr,"                    ... Pita Time Information ...                    \n"); 
   filePrint(fptr,"*********************************************************************\n\n");
  
   filePrint(fptr,"1.  PitaTotaltime                   Min          Max          Avg   \n");
   filePrint(fptr,"                                   %.3es  %.3es %.3es \n\n",
        PitaTotaltime_Min/1000.0, PitaTotaltime_Max/1000.0, PitaTotaltime_Avg/1000.0);

   filePrint(fptr,"2.  buildTStime                     Min          Max          Avg   \n");
   filePrint(fptr,"                                   %.3es  %.3es %.3es \n\n",
        buildTStime_Min/1000.0, buildTStime_Max/1000.0, buildTStime_Avg/1000.0);

   filePrint(fptr,"***********************************************************************\n\n");
   filePrint(fptr,"3.  dynOpstime                      Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es   %.3es \n\n",
        dynOpstime_Min/1000.0, dynOpstime_Max/1000.0, dynOpstime_Avg/1000.0);

   filePrint(fptr,"*********************************************************************\n\n");

   filePrint(fptr,"4.  step0time                       Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es   %.3es \n\n",
        step0time_Min/1000.0, step0time_Max/1000.0, step0time_Avg/1000.0);

   filePrint(fptr,"5.  Bitime                          Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es   %.3es \n\n",
        Bitime_Min/1000.0, Bitime_Max/1000.0, Bitime_Avg/1000.0);

   filePrint(fptr,"6.  exchangeProptime                Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es   %.3es \n\n",
        exchangeProptime_Min/1000.0, exchangeProptime_Max/1000.0, exchangeProptime_Avg/1000.0);

   filePrint(fptr,"*********************************************************************\n\n");
   filePrint(fptr,"7.  CkTotaltime                     Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es   %.3es \n\n",
        CkTotaltime_Min/1000.0, CkTotaltime_Max/1000.0, CkTotaltime_Avg/1000.0);

   filePrint(fptr,"8.  basetime                        Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es   %.3es \n\n",
        basetime_Min/1000.0, basetime_Max/1000.0, basetime_Avg/1000.0);

   filePrint(fptr,"9.  CkfineProjtime                  Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es  %.3es \n\n",
        CkfineProjtime_Min/1000.0, CkfineProjtime_Max/1000.0, CkfineProjtime_Avg/1000.0);

   filePrint(fptr,"10. CkcoarseProjtime                Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es   %.3es \n\n",
        CkcoarseProjtime_Min/1000.0, CkcoarseProjtime_Max/1000.0, CkcoarseProjtime_Avg/1000.0);
   filePrint(fptr,"**********************************************************************\n\n");

   filePrint(fptr,"11. newCondtime                     Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es   %.3es \n\n",
        newCondtime_Min/1000.0, newCondtime_Max/1000.0, newCondtime_Avg/1000.0);

   filePrint(fptr,"12. solveITA                        Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es   %.3es \n\n",
        solveITA_Min/1000.0, solveITA_Max/1000.0, solveITA_Avg/1000.0);

   filePrint(fptr,"13. StoreFineGridTime               Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es   %.3es \n\n",
        StoreFineGridTime_Min/1000.0, StoreFineGridTime_Max/1000.0, StoreFineGridTime_Avg/1000.0);

   filePrint(fptr,"*********************************************************************\n\n");

   filePrint(fptr,"14. getInitCondtime                 Min           Max           Avg   \n");
   filePrint(fptr,"                                   %.3es   %.3es   %.3es \n\n",
        getInitCondtime_Min/1000.0, getInitCondtime_Max/1000.0, getInitCondtime_Avg/1000.0);

   filePrint(fptr,"*********************************************************************\n\n");
   
   filePrint(fptr,"15. TStime                          %.3es \n\n", TStime_Max/1000.0);
   
   filePrint(fptr,"*********************************************************************\n\n");
   filePrint(fptr,"16 exchange_info_Max          %.3es   \n\n", exchange_info_Max/1000.0);
   filePrint(fptr,"17 compNorm_Max               %.3es   \n\n", compNorm_Max/1000.0);
   filePrint(fptr,"18 jump_Max                   %.3es   \n\n", jump_Max/1000.0);
   filePrint(fptr,"19 newTSid                    %.3es   \n\n", newTSid_Max/1000.0);
   filePrint(fptr,"*********************************************************************\n\n");

   if (fptr) { fclose(fptr); fptr = 0; }
 
   }

#endif   


}

