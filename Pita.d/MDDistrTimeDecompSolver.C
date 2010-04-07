#ifndef _MDDISTRTIMEDECOMPSOLVER_C_
#define _MDDISTRTIMEDECOMPSOLVER_C_

#include <Pita.d/MDDistrTimeDecompSolver.h>
#include <Driver.d/Communicator.h>
#include <Comm.d/Communicator.h>
#include <cstdlib>

using std::exit;

extern Communicator *structCom;

//----------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
MDDistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::MDDistrTimeDecompSolver(ProblemDescriptor * _probDesc):DistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>(_probDesc)
{

  /* Each MPI Process has to communicate with:

        - MPI Process that have the same space subdomains but differents 
          time subdomains to exchange the vectors yi(Ti+1). 

          For these communications, use timeCom.    

        - MPI Process that have the time subdomains but differents space 
          subdomain to exchange infos used by FETI. 

          For these communications, use structCom.
           
  */

  // numSpaceMPIProc == number of MPI Process that will be used to solve
  // the problem in space.

  int rank, totalnumCPU;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnumCPU);            
                                                                                                                              
  numSpaceMPIProc = _probDesc->getDomain()->solInfo().numSpaceMPIProc;

  if ( totalnumCPU%numSpaceMPIProc!=0 ){
     fprintf(stderr," !!!!!!!! total number of CPU isn't can't be divided by numSpaceMPIProc !!!!!!! \n");
     exit(-1);
  }

  int spaceColor = rank%numSpaceMPIProc;
  int timeColor  = rank/numSpaceMPIProc;
                                                                                                                                          
  MPI_Comm comm1,comm2;
                                                                                                                                          
  MPI_Comm_split(MPI_COMM_WORLD, spaceColor+1, rank, &comm1);
  MPI_Comm_split(MPI_COMM_WORLD, timeColor+1, rank, &comm2);
                                                                                                                                          
  if (structCom) delete structCom;

  structCom = new Communicator(comm1,stderr);
  timeCom   = new Communicator(comm2,stderr);

}

//----------------------------------------------------------------------------------
template<
    class DynOps,
    class VecType,
    class PostProcessor,
    class ProblemDescriptor,
    class InfoSize>
MDDistrTimeDecompSolver<DynOps,VecType,PostProcessor,ProblemDescriptor,InfoSize>::~MDDistrTimeDecompSolver()
{
  //delete only timeCom, structCom is deleted in main.C
  if (timeCom)  delete timeCom;
}

//----------------------------------------------------------------------------------

#endif
