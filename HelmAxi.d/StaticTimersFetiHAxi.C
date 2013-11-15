#include <cstdio>
#include <Utils.d/dbg_alloca.h>

#include <Timers.d/Timing.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/MatrixTimers.h>
#include <Driver.d/Domain.h>
#include <Threads.d/Paral.h>
#include <Utils.d/Memory.h>
#include <Comm.d/Communicator.h>
#include <Utils.d/DistHelper.h>
#include <HelmAxi.d/FourierHelmBCs.h>
#include <Driver.d/GeoSource.h>

extern long totMemSparse;
extern long totMemSky;

extern char* message[];
extern char* yesno[];

//const double oneMegaByte = (1024.0*1024.0);
const double byteToMb    = (1.0 / oneMegaByte);

extern char* renumMessage[];
extern char* scalingMessage[];
extern char* rbmMessage[];
extern char* precMessage[];
extern char* projectMessage[];
extern char* subSolverMessage[];
extern char* precSolverMessage[];
extern char* gtgSolverMessage[];
extern char* solverMessage[];

// This function prints timers for FETI-H AXI Solver 

void
StaticTimers::printStaticTimersFetiHAxi(MatrixTimers matrixTimer, 
              double solveTime, SolverInfo& sInfo, Timings& timers,
              ControlInfo& cinfo, Domain* domain, FourierHelmBCs *glBC, 
              int MPCsize) {

 openTimingFile(cinfo);

 int numCPUs       = threadManager->numThr();
 int numSubdomains = timers.numSubdomains;

 if (f != 0) {
 
   int numCon = glBC->numDir;

   fprintf(f,"\n***********************************************************"
             "********************\n");
   fprintf(f," ... Problem Information ... \n");
   fprintf(f,"***********************************************************"
             "********************\n\n");
  
   fprintf(f,"1. Number of Nodes                         = %14d\n\n",
              domain->numNode());
  
   fprintf(f,"2. Number of Elements                      = %14d\n\n",
              domain->numElements());

   // Special case of one dof per node
 
   fprintf(f,"3. Number of Degrees of Freedom            = %14d\n",
              domain->numNode());
   fprintf(f,"         Number of Constrained Dofs        = %14d\n",numCon);
   fprintf(f,"         Number of Unconstrained Dofs      = %14d\n\n",
              domain->numNode() - numCon);
 
   fprintf(f,"4. Number of Applied Loads                 = %14d\n\n",
              glBC->numNeu);
 
   fprintf(f,"5. Number of Sommerfeld Elements           = %14d\n\n",
              glBC->numSomm);

   fprintf(f,"6. Number of Fourier Coefficients          = %14d\n\n",
              2*glBC->numModes+1);

   if (MPCsize>0) {
     fprintf(f,"7. Number of MPC                           = %14d\n\n",
                MPCsize);
  
     fprintf(f,"8. Number of Output Files                  = %14d\n\n",
                geoSource->getNumOutInfo());
  
     fprintf(f,"9. Renumbering                             = %14s\n\n",
             renumMessage[sInfo.renum]);
   }
   else { 
     fprintf(f,"7. Number of Output Files                  = %14d\n\n",
              geoSource->getNumOutInfo());
  
     fprintf(f,"8. Renumbering                             = %14s\n\n",
             renumMessage[sInfo.renum]);
   }
  
   fprintf(f,"***********************************************************"
             "********************\n");
   fprintf(f," ... Solver Information ... \n");
   fprintf(f,"***********************************************************"
             "********************\n\n");

   fprintf(f,"1. FETI-H\n");

   if (timers.numCRNs) {
     fprintf(f,"         Two-Level\n");
     switch (sInfo.getFetiInfo().nonLocalQ) {
       case 0:
         fprintf(f,"         Basic Projector\n");
         break;
       case 1:
         fprintf(f,"         Fourier Projector\n");
         break;
     }
   } else
        fprintf(f,"         One-Level\n");

   fprintf(f,"         %s", precMessage[sInfo.getFetiInfo().precno]);

   fprintf(f,"         %s", subSolverMessage[sInfo.getFetiInfo().solvertype]);
   
   if (sInfo.getFetiInfo().precno)
   fprintf(f,"         %s", precSolverMessage[sInfo.getFetiInfo().solvertype]);

   fprintf(f,"         Coarse Solver Selected            =        %s",
              gtgSolverMessage[sInfo.getFetiInfo().gtgSolver]);
 
   fprintf(f,"         Zero Pivot Tolerance for Local K  = %14.3e\n",
             sInfo.trbm);
 
   fprintf(f,"         Total Size of QtFQ                = %14d\n",
             timers.numCRNs);
 
   fprintf(f,"         Zero Pivot Tolerance for QtFQ     = %14.3e\n",
             sInfo.getFetiInfo().grbm_tol);

   fprintf(f,"         Number of Linearly Depend. Modes  = %14d\n",
             timers.numRBMs);
 
   fprintf(f,"         Maximum Number of Iterations      = %14d\n",
             sInfo.getFetiInfo().maxiter());

   fprintf(f,"         Maximum Size of Reortho. Vectors  = %14d\n",
             sInfo.getFetiInfo().maxorth());

   fprintf(f,"         Tolerance for Convergence         = %14.3e\n",
             sInfo.getFetiInfo().tolerance());
 
   fprintf(f,"\n***********************************************************"
             "********************\n");
   fprintf(f," ... Timing Statistics for %d Threads and %d Subdomains ...\n",
                  numCPUs,numSubdomains);
   fprintf(f,"***********************************************************"
             "********************\n\n");
  
 }

 //////////////////////////////////////////////////
 //          Preprocessing items                 // 
 //////////////////////////////////////////////////

 // Begin timer output

 memoryPrecond       = timers.memoryPrecond;

 TimeData assembleMin  = timers.assembleMat.getMin();
 TimeData assembleMax  = timers.assembleMat.getMax();
 TimeData assembleAvg  = timers.assembleMat.getAvg();
 TimeData assembleTot  = timers.assembleMat.getTot();

 TimeData buildRhsMin  = timers.buildRhs.getMin();
 TimeData buildRhsMax  = timers.buildRhs.getMax();
 TimeData buildRhsAvg  = timers.buildRhs.getAvg();
 TimeData buildRhsTot  = timers.buildRhs.getTot();

 long totMemSubMatrices = timers.memorySubMatrices - 16*memoryPrecond 
                               - 16*memoryK;

 long totalMemRead = matrixTimer.memoryParse + matrixTimer.memorySetUp;

 // Timer sub totals

 double subTotal[6];

 subTotal[0] = (matrixTimer.readTime + matrixTimer.readDecomp);

 subTotal[1] = preProcess + matrixTimer.setUpDataTime 
                          + corotatorTime + kelArrayTime + timeGeom
			  - matrixTimer.readDecomp
                          - matrixTimer.createDofs;

 subTotal[2] = assembleTot.time;

 subTotal[3] = buildRhsTot.time;

 if (f != 0) {

   fprintf(f,"1. Total Read Input Files              time: %14.5f s %14.3f Mb"
             "\n", subTotal[0]/1000.0,totalMemRead*byteToMb);

   fprintf(f,"         Read Mesh                     time: %14.5f s\n",
             matrixTimer.readTime/1000.0);

   fprintf(f,"         Read Mesh Partition           time: %14.5f s\n",
             matrixTimer.readDecomp/1000.0);          

   fprintf(f,"\n");

   fprintf(f,"2. Total Preprocessing                 time: %14.5f s %14.3f Mb"
             "\n",subTotal[1]/1000.0, 
                  (memoryPreProcess-matrixTimer.memoryForm-memoryRhs)*byteToMb);

   fprintf(f,"         Process Input Data            time: %14.5f s %14.3f Mb"
             "\n", matrixTimer.setUpDataTime/1000.0,
                   matrixTimer.memorySetUp*byteToMb);

   fprintf(f,"         Make Subdomains               time: %14.5f s %14.3f Mb"
             "\n", matrixTimer.makeSubDomains/1000.0,
                   matrixTimer.memorySubdomain*byteToMb);

   fprintf(f,"            Element to Node Connectivity   : %31.3f Mb\n", 
             matrixTimer.memoryElemToNode*byteToMb);

   fprintf(f,"            Sub. to Node Connectivity      : %31.3f Mb\n", 
             matrixTimer.memorySubToNode*byteToMb);

   fprintf(f,"            Node to Sub. Connectivity      : %31.3f Mb\n", 
             matrixTimer.memoryNodeToSub*byteToMb);

   if (matrixTimer.memoryMPCToNode>0) {

     fprintf(f,"            MPC to Node Connectivity       : %31.3f Mb\n", 
               matrixTimer.memoryMPCToNode*byteToMb);

     fprintf(f,"            MPC to Sub. Connectivity       : %31.3f Mb\n", 
               matrixTimer.memoryMPCToSub*byteToMb);

   }

   fprintf(f,"         Distribute BCs                time: %14.5f s %14.3f Mb"
             "\n", matrixTimer.distributeBCs/1000.0,
                   matrixTimer.memoryDistBC*byteToMb);

   if (matrixTimer.memoryMPCToNode>0)
     fprintf(f,"         Distribute MPCs               time: %14.5f s %14.3f"
             " Mb\n", matrixTimer.distributeMPCs/1000.0,
                   matrixTimer.memoryDistMPC*byteToMb);

   fprintf(f,"         Make Connectivities           time: %14.5f s %14.3f Mb"
             "\n", matrixTimer.makeConnectivity/1000.0,
                   matrixTimer.memoryConnect*byteToMb);

   fprintf(f,"         Make Interface                time: %14.5f s %14.3f Mb"
             "\n", matrixTimer.makeInterface/1000.0,
                   matrixTimer.memoryInterface*byteToMb);

   fprintf(f,"         Make Polygons                 time: %14.5f s %14.3f Mb"
             "\n", matrixTimer.makeInternalInfo/1000.0,
                   matrixTimer.memoryInternal*byteToMb);

   fprintf(f,"\n");

   fprintf(f,"3. Total Subdomain Matrices Processing time: %14.5f s %14.3f Mb"
             "\n\n", subTotal[2]/1000.0,totMemSubMatrices*byteToMb);

   fprintf(f,"4. Total Subdomain RHS Processing      time: %14.5f s %14.3f Mb"
             "\n\n", subTotal[3]/1000.0, buildRhsMax.memory*byteToMb);
   
 }

 //////////////////////////////////////////////////
 //          Solver items                        // 
 //////////////////////////////////////////////////

 double solutionTime = timers.solve + getFetiSolverTime ;

 TimeData constructMin = timers.consMatrix.getMin();
 TimeData constructMax = timers.consMatrix.getMax();
 TimeData constructAvg = timers.consMatrix.getAvg();
 TimeData constructTot = timers.consMatrix.getTot();

 double coarse1Max  = timers.coarse1;
 double coarse1Tot  = timers.coarse1;
 double coarse1Min  = timers.coarse1;

 double parfac1Max  = timers.pfactor;
 double parfac1Tot  = timers.pfactor;
 double parfac1Min  = timers.pfactor;
  
 double precondMax   = timers.precond;

 double sAndJMaximum = timers.sAndJ;

 TimeData factorMin    = timers.factorMat.getMin();
 TimeData factorMax    = timers.factorMat.getMax();
 TimeData factorAvg    = timers.factorMat.getAvg();
 TimeData factorTot    = timers.factorMat.getTot();

 TimeData sAndJMin     = timers.solveAndJump.getMin();
 TimeData sAndJMax     = timers.solveAndJump.getMax();
 TimeData sAndJAvg     = timers.solveAndJump.getAvg();
 TimeData sAndJTot     = timers.solveAndJump.getTot();

 TimeData precMin      = timers.preconditioner.getMin();
 TimeData precMax      = timers.preconditioner.getMax();
 TimeData precAvg      = timers.preconditioner.getAvg();
 TimeData precTot      = timers.preconditioner.getTot();

 TimeData orthoMin     = timers.orthogonalize.getMin();
 TimeData orthoMax     = timers.orthogonalize.getMax();
 TimeData orthoAvg     = timers.orthogonalize.getAvg();
 TimeData orthoTot     = timers.orthogonalize.getTot();

 double proj1Seq = timers.project - timers.sAndJ;

 // Timer sub totals

 subTotal[4] = solutionTime - subTotal[2] + matrixTimer.createDofs 
               + timers.coarse2 + timers.pfactor2;

 subTotal[5] = output;

 double coarseTime = coarse1Max;

 long totMemCoarse = timers.memoryGtG + memoryRhs;

 long totalMemFactor = factorTot.memory + 16*memoryK;

 long totalSolverMemory = timers.memoryFETI + timers.memorySolve
                               - totMemSubMatrices + 
                               matrixTimer.memoryForm;

 if (f != 0) {

   fprintf(f,"5. Total Solver                        time: %14.5f s %14.3f Mb"
             "\n\n",subTotal[4]/1000.0, totalSolverMemory*byteToMb);

   fprintf(f,"         Factor Subdomain Matrices     time: %14.5f s %14.3f Mb"
             "\n\n",factorMax.time/1000.0, totalMemFactor*byteToMb);

   fprintf(f,"         Total Building Coarse Pbs.    time: %14.5f s %14.3f Mb"
             "\n",(coarseTime+matrixTimer.createDofs)/1000.0, 
                  (totMemCoarse+matrixTimer.memoryForm)*byteToMb );

   fprintf(f,"               Preproc. Coarse Data    time: %14.5f s %14.3f Mb"
             "\n", matrixTimer.createDofs/1000.0,
                   matrixTimer.memoryForm*byteToMb);

   fprintf(f,"               Local Coarse Pbs.       time: %14.5f s "
             "\n\n", constructTot.time/1000.0);

   if (matrixTimer.memoryMPCToNode>0) {
     fprintf(f,"         Total Building MPC Solver     time: %14.5f s %14.3f"
             " Mb \n", (timers.coarse2+timers.pfactor2)/1000.0, 
             timers.memoryPCtFPC*byteToMb);
     fprintf(f,"               Factorization           time: %14.5f s\n",
               timers.pfactor2/1000.0);
     fprintf(f,"\n");
   }

   fprintf(f,"         Total Paral. Fac. Coarse Pbs. time: %14.5f s %14.3f Mb"
           "\n\n", parfac1Max/1000.0, timers.memoryFactor*byteToMb);

   fprintf(f,"         Total Solve Loop              time: %14.5f s %14.3f Mb"
             "\n", timers.solve/1000.0, timers.memorySolve*byteToMb);

   fprintf(f,"               Projection              time: %14.5f s %14.3f Mb"
             "\n", proj1Seq/1000.0, 
                  (timers.memoryProject1-timers.memorySAndJ)*byteToMb);

   fprintf(f,"               Precondition            time: %14.5f s %14.3f Mb"
             "\n", precondMax/1000.0, 16*memoryPrecond*byteToMb);

   fprintf(f,"               Local Solve             time: %14.5f s %14.3f Mb"
             "\n", timers.sAndJ/1000.0, 
                   timers.memorySAndJ*byteToMb);

   fprintf(f,"               Reorthogonalize         time: %14.5f s %14.3f Mb"
             "\n", orthoMax.time/1000.0, timers.memoryOSet*byteToMb);

   fprintf(f,"\n");

   fprintf(f,"6. Write Output Files                  time: %14.5f s %14.3f Mb"
           "\n",subTotal[5]/1000.0, memoryOutput*byteToMb);

   // Compute the total time spent on this simulation

   double total = 0.0;

   int i;
   for (i=0; i<6; ++i)
     total += subTotal[i];

   long totalMemSimulation = totalSolverMemory + buildRhsTot.memory +
                                  totMemSubMatrices + memoryOutput +       
                                  memoryPreProcess + totalMemRead;

   fprintf(f,"\n");

   fprintf(f,"TOTAL SIMULATION (1+2+3+4+5+6)         time: %14.5f s %14.3f Mb"
             "\n",total/1000.0, totalMemSimulation*byteToMb);

 //////////////////////////////////////////////////
 //      Output FETI Solver information          // 
 //////////////////////////////////////////////////

   if (sInfo.type == 2) {

      fprintf(f,"\n");

      fprintf(f,"***********************************************************"
             "********************\n");
      fprintf(f," ... FETI Monitoring ... \n");
      fprintf(f,"***********************************************************"
             "********************\n\n");

      fprintf(f,"1. Total Amount of Requested Memory        = %14.3f Mb\n\n",
           (memoryUsed())*byteToMb);

      fprintf(f,"2. FETI Solver Amount of Requested Memory  = %14.3f Mb\n\n",
                timers.memoryFETI*byteToMb);

      // Check whether we converged or not

      if (timers.converged == 1)
         fprintf(f,"3. Number of Iterations for Convergence    = %14d\n\n",
                   timers.numIter);
      else if (timers.converged == 0) {
              fprintf(f,"3. Stagnation Occured After a # of iter    = %14d",
                         timers.numIter);
              fprintf(f,"\n\n");
           } else {
              fprintf(f,"3. Maximum Number of Iterations reached    = %14d",
                         timers.numIter);
              fprintf(f,"\n\n");
           }

      fprintf(f,"4. Relative Primal Error Reached           = %14.3e\n\n",
                timers.iterations[0].finalPrimal);

      fprintf(f,"5. Relative Dual Error Reached             = %14.3e\n\n",
                timers.iterations[0].finalDual);

      fprintf(f,"6. Size of Coarse Problem                  = %14d %14.3f Mb"
              "\n\n", timers.numCRNs, timers.memoryGtGsky*byteToMb);

      fprintf(f,"7. Maximum dependencies in Coarse Pb.      = %14d "
              "\n\n", timers.numRBMs);

      if (sInfo.getFetiInfo().solvertype == 0)
         fprintf(f,"8. Total Memory Subdomain Skyline K        = %14.3f Mb\n\n",
                 16.0*totMemSky*byteToMb);
      else if (sInfo.getFetiInfo().solvertype == 1)
             fprintf(f,"8. Total Memory Subdomain Sparse K         = %14.3f Mb"
                       "\n\n", 16.0*totMemSparse*byteToMb);
           else
             fprintf(f,"\n"); // if we have other subdomain solvers

      fprintf(f,"\n********************************************************"
             "***********************\n");

   }

 }
 
 ////////////////////////////////////////////////////////////
 //           Detailed CPU Statistics                      //
 ////////////////////////////////////////////////////////////

 if(f) {

   fprintf(f," ... Detailed CPU Statistics (Seconds) ");
   fprintf(f,"\n***********************************************************"
             "********************\n");
   fprintf(f,"\n                                             minimum      "
             "average      maximum\n");

   double tot1MinTime = matrixTimer.readTime + matrixTimer.readDecomp;
   double tot1AvgTime = matrixTimer.readTime + matrixTimer.readDecomp;
   double tot1MaxTime = matrixTimer.readTime + matrixTimer.readDecomp;

   fprintf(f,"\n1. Total Read Input Files             : %12.4f %12.4f %12.4f"
           "\n", tot1MinTime/1000.0, tot1AvgTime/1000.0, tot1MaxTime/1000.0);
	
   double tot2MinTime = subTotal[1];
   double tot2AvgTime = subTotal[1];
   double tot2MaxTime = subTotal[1];
	 
   fprintf(f,"\n2. Total Preprocessing                : %12.4f %12.4f %12.4f"
           "\n", tot2MinTime/1000.0,tot2AvgTime/1000.0,tot2MaxTime/1000.0);
	 
   fprintf(f,"         Process Input Data           : %12.4f %12.4f %12.4f\n",
           matrixTimer.setUpDataTime/1000.0,
           matrixTimer.setUpDataTime/1000.0,
           matrixTimer.setUpDataTime/1000.0);

   fprintf(f,"         Make Subdomains              : %12.4f %12.4f %12.4f\n", 
           matrixTimer.makeSubDomains/1000.0,
           matrixTimer.makeSubDomains/1000.0,
           matrixTimer.makeSubDomains/1000.0);

   fprintf(f,"         Distribute BCs               : %12.4f %12.4f %12.4f\n", 
           matrixTimer.distributeBCs/1000.0,
           matrixTimer.distributeBCs/1000.0,
           matrixTimer.distributeBCs/1000.0);

   if (matrixTimer.memoryMPCToNode>0)
      fprintf(f,"         Distribute MPCs              : %12.4f %12.4f %12.4f"
           "\n", matrixTimer.distributeMPCs/1000.0,
           matrixTimer.distributeMPCs/1000.0,
           matrixTimer.distributeMPCs/1000.0);

   fprintf(f,"         Make Connectivities          : %12.4f %12.4f %12.4f\n", 
           matrixTimer.makeConnectivity/1000.0,
           matrixTimer.makeConnectivity/1000.0,
           matrixTimer.makeConnectivity/1000.0);

   fprintf(f,"         Make Interface               : %12.4f %12.4f %12.4f\n", 
           matrixTimer.makeInterface/1000.0,
           matrixTimer.makeInterface/1000.0,
           matrixTimer.makeInterface/1000.0);

   fprintf(f,"         Make Polygons                : %12.4f %12.4f %12.4f\n",
         matrixTimer.makeInternalInfo/1000.0,
         matrixTimer.makeInternalInfo/1000.0,
         matrixTimer.makeInternalInfo/1000.0);

   double tot3MinTime = assembleMin.time;
   double tot3AvgTime = assembleAvg.time;
   double tot3MaxTime = assembleMax.time; 

   fprintf(f,"\n3. Total Subdomain Matrices Processing: %12.4f %12.4f %12.4f"
             "\n",tot3MinTime/1000.0,tot3AvgTime/1000.0,tot3MaxTime/1000.0);
 
   double tot4MinTime = buildRhsMin.time;
   double tot4AvgTime = buildRhsAvg.time;
   double tot4MaxTime = buildRhsMax.time;

   fprintf(f,"\n4. Total Subdomain RHS Processing     : %12.4f %12.4f %12.4f"
           "\n\n",tot4MinTime/1000.0,tot4AvgTime/1000.0,tot4MaxTime/1000.0);

   // This needs to be improved to count better the parallel/sequential time

   double timeCoarseMin = coarse1Min + matrixTimer.createDofs;
   double timeCoarseTot = coarse1Tot + matrixTimer.createDofs;
   double timeCoarseMax = coarse1Max + matrixTimer.createDofs; 

   double timeParFacMin = parfac1Min;
   double timeParFacTot = parfac1Tot;
   double timeParFacMax = parfac1Max;

   double com1 = precondMax   - precMax.time;
   double com2 = sAndJMaximum - sAndJMax.time;

   double proj1Min = proj1Seq;
   double proj1Avg = proj1Seq;
   double proj1Max = proj1Seq;
 
   double solveMin =  proj1Min +
                      precMin.time + com1 + sAndJMin.time + com2 + 
                      orthoMin.time;
   double solveAvg =  proj1Avg +
                      precAvg.time + com1 + sAndJAvg.time + com2 + 
                      orthoAvg.time;
   double solveMax =  proj1Max +  
                      precMax.time + com1 + sAndJMax.time + com2 + 
                      orthoMax.time;

   double missingTime = timers.solve - solveMax;

   solveMin += missingTime;
   solveAvg += missingTime;
   solveMax += missingTime;

   double tot5MinTime   = timeCoarseMin + factorMin.time + timeParFacMin
                        + solveMin + timers.coarse2+timers.pfactor2; 

   double tot5AvgTime   = factorAvg.time + timeCoarseTot
                        + timeParFacTot  + solveAvg
                        + timers.coarse2+timers.pfactor2;

   double tot5MaxTime   = timeCoarseMax + factorMax.time + timeParFacMax
                        + solveMax + timers.coarse2+timers.pfactor2;
   
   fprintf(f,"5. Total Solver                       : %12.4f %12.4f %12.4f"
           "\n\n",tot5MinTime/1000.0, tot5AvgTime/1000.0, tot5MaxTime/1000.0);

   fprintf(f,"         Factor Subdomain Matrices    : %12.4f %12.4f %12.4f"
         "\n", factorMin.time/1000.0, factorAvg.time/1000.0,
                factorMax.time/1000.0);

   fprintf(f,"         Total Building  Coarse Pbs.  : %12.4f %12.4f %12.4f\n",
           timeCoarseMin/1000.0,timeCoarseTot/1000.0,timeCoarseMax/1000.0);

   fprintf(f,"               Preproc. Coarse Data   : %12.4f %12.4f %12.4f"
           "\n", matrixTimer.createDofs/1000.0,
           matrixTimer.createDofs/1000.0, matrixTimer.createDofs/1000.0);

   fprintf(f,"               Local Coarse Pbs.      : %12.4f %12.4f %12.4f"
           "\n", constructMin.time/1000.0, 
           constructAvg.time/1000.0, constructMax.time/1000.0);

   if (matrixTimer.memoryMPCToNode>0) {
     fprintf(f,"         Total Building MPC Solver    : %12.4f %12.4f "
             "%12.4f\n", (timers.coarse2+timers.pfactor2)/1000.0, 
             (timers.coarse2+timers.pfactor2)/1000.0,
             (timers.coarse2+timers.pfactor2)/1000.0);
   }

   fprintf(f,"         Total Paral. Fac. Coarse Pbs.: %12.4f %12.4f %12.4f"
           "\n", timeParFacMin/1000.0,timeParFacTot/1000.0,
                timeParFacMax/1000.0);

   fprintf(f,"         Total Solve loop             : %12.4f %12.4f %12.4f"
           "\n",solveMin/1000.0,solveAvg/1000.0,solveMax/1000.0);

   fprintf(f,"               Projection             : %12.4f %12.4f %12.4f"
           "\n",proj1Seq/1000.0,proj1Seq/1000.0,proj1Seq/1000.0);

   fprintf(f,"               Precondition           : %12.4f %12.4f %12.4f"
          "\n",(precMin.time+com1)/1000.0,(precAvg.time+com1)/1000.0,
          precondMax/1000.0);

   fprintf(f,"               Local Solve            : %12.4f %12.4f %12.4f"
         "\n", (sAndJMin.time+com2)/1000.0,(sAndJAvg.time+com2)/1000.0,
          sAndJMaximum/1000.0);

   fprintf(f,"               Reorthogonalize        : %12.4f %12.4f %12.4f"
           "\n", orthoMin.time/1000.0, orthoAvg.time/1000.0,
           orthoMax.time/1000.0);

   double tot6MinTime = output;
   double tot6AvgTime = output;
   double tot6MaxTime = output;

   fprintf(f,"\n6. Write Output Files                 : %12.4f %12.4f %12.4f"
          "\n", tot6MinTime/1000.0,tot6AvgTime/1000.0,tot6MaxTime/1000.0);

   double timeSimMin = tot1MinTime + tot2MinTime + tot3MinTime + tot4MinTime
                     + tot5MinTime + tot6MinTime;
   double timeSimAvg = tot1AvgTime + tot2AvgTime + tot3AvgTime + tot4AvgTime
                     + tot5AvgTime + tot6AvgTime;
   double timeSimMax = tot1MaxTime + tot2MaxTime + tot3MaxTime + tot4MaxTime
                     + tot5MaxTime + tot6MaxTime;
 
   fprintf(f,"\nTOTAL SIMULATION (1+2+3+4+5+6)        : %12.4f %12.4f %12.4f"
             "\n",timeSimMin/1000.0, timeSimAvg/1000.0,timeSimMax/1000.0);
 }

}





