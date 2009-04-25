#ifdef STRUCTOPT

#include <Utils.d/dbg_alloca.h>
#include <boost/scoped_array.hpp>
#include <stdio.h>

#include <Math.d/mathUtility.h>

#include <Structopt.d/Optfilter.h>
#include <Structopt.d/Structopt_sd.h>

#ifdef USE_MPI
#include <Comm.d/Communicator.h>
extern Communicator *structCom;
#endif

//-----------------------------------------------------------------------------

void Optfilter::print(int* critList, int numCrit,FILE * optunitout)
{

   fprintf(optunitout,"\n\tFilter for gradients of criteria :");
   int i;

    for (i=0;i<numCrit;i++) 
       if (critList[i] != 0 ) fprintf(optunitout," %d ( %d ) ",i+1,critList[i]);

    fprintf(optunitout," \n\n");

    fprintf(optunitout,"\t(1) scaling for maximum  (-1) scaling for minimum\n");
   
    fprintf(optunitout,"\tFiter method ................... : ");
    
    switch(typ)
    {
       case 0:
         fprintf(optunitout,"Maute\n");
         break;
       case 1:
         fprintf(optunitout,"Sigmund without variable scaling \n");
         break;
       case 2:
         fprintf(optunitout,"Sigmund with variable scaling\n");
         break;
    }

    fprintf(optunitout," \n\n");

    fprintf(optunitout,"\tFiter radius ................... : %e\n",refRadius);    
    fprintf(optunitout,"\tFiter power (min)............... : %e\n",minExp);    
    fprintf(optunitout,"\tFiter power (max)............... : %e\n",maxExp);    
    fprintf(optunitout,"\tNumber of adaption steps........ : %d\n",static_cast<int>(maxCount));    

    fprintf(optunitout,"\tScaling of gradients............ : ");

    switch(sclTyp)
    {
       case 0:
         fprintf(optunitout,"no scaling\n");
         break;
       case 1:
         fprintf(optunitout,"scaling for L2 norm of gradients \n");
         break;
       case 2:
         fprintf(optunitout,"scaling for maximum/minium value\n");
         break;
    }

    fprintf(optunitout,"\n");

    fprintf(optunitout,"\tNumber of groups................ : %d\n\n",numGroups);    

    for (i=0;i<numGroups;i++) 
         fprintf(optunitout,"\tgroup %d:  first=%d   last=%d\n",i+1,
	                                smoGroups[i][0]+1,smoGroups[i][1]);
    fprintf(optunitout,"\n");

}

//-----------------------------------------------------------------------------

Optfilter::Optfilter(int _typ, int _sclTyp, double& _refRadius, 
                     double& _maxCount, double& _minExp, double& _maxExp,
		     int _numGroups, int listsize, int *list) 
{

   refRadius = _refRadius;
   typ       = _typ;
   sclTyp    = _sclTyp;

   Weights   = 0;
   Connected = 0;
   numNeighb = 0;

   counter  = 0;
   maxCount = _maxCount;

   minExp   = _minExp;
   maxExp   = _maxExp;

   numGroups = 0;
   
   if ( _numGroups ) {
   
     numGroups = _numGroups;

     if ( listsize != (2*numGroups) ) 
       fprintf(stderr," .... WARNING: check Filter grouping\n");

     smoGroups = new int*[numGroups];

     int i;
     for (i=0;i<numGroups;i++) {
       smoGroups[i]    = new int[2];
       smoGroups[i][0] = list[i*2];
       smoGroups[i][1] = list[i*2+1]+1;
     }
   }

} 

//-----------------------------------------------------------------------------
void Optfilter::project(double* vector,int xFlag)
{ 
  double expFac = std::min(minExp + pow(counter/maxCount,4.0) * (maxExp-minExp),
			   maxExp);
   if(typ == 0 )
     { filePrint(stderr," ... gradient exponent : %e\n",expFac); }

   int myStart, myEnd;
#ifdef USE_MPI
   const int locSize = static_cast<int>(ceil(size/structCom->numCPUs()));
   myStart = structCom->myID()*locSize;
   myEnd   = std::min((structCom->myID()+1)*locSize, size);
#else
   myStart = 0; myEnd = size;
#endif
   boost::scoped_array<double> work((myEnd>myStart)? new double[myEnd-myStart] : 0);

   switch(typ) 
     {
     case 0:
       for(int i=myStart; i<myEnd; ++i)
	 {
	   work[i-myStart] = 0;
	   double scale   = 0;
	   for(int j=0; j<numNeighb[i]; ++j)
	     { scale +=  pow(Weights[i][j],expFac); }
	   for(int j=0; j<numNeighb[i]; ++j)
	     {
	       int nn = Connected[i][j];
	       work[i-myStart] += pow(Weights[i][j],expFac)*vector[nn]/scale;
	     }
	 }
       break;
     case 1:
       for(int i=myStart; i<myEnd; ++i) 
	 {
	   work[i-myStart] = 0;
	   for(int j=0; j<numNeighb[i]; ++j)
	     {
	       int nn = Connected[i][j];
	       work[i-myStart] += Weights[i][j]*vector[nn];
	     }
	 }
       break;
     case 2:
       for(int i=myStart; i<myEnd; ++i) 
	 {
	   work[i-myStart] = 0;
	   for(int j=0; j<numNeighb[i]; j++)
	     {
	       int nn = Connected[i][j];
	       work[i-myStart] += pow((variables[nn]/variables[i]),expFac)*Weights[i][j]*
		 vector[nn];
	     }
	 }
       break;
     default:
       assert(0);
       break;
     }
   // filter each group individually
   // variable not belonging to any group are not modified
   for(int ig=0; ig<numGroups; ++ig) 
     {
       double orgNorm=0.0;
       double newNorm=0.0;
       
       double orgMax=-1.0e20;
       double orgMin= 1.0e20;
       double newMax=-1.0e20;
       double newMin= 1.0e20;
       double orgSum= 0;
       double newSum= 0;

       double scale=1.0;
       double shift=0.0;
       
       // get group     
       int min2 = std::max(smoGroups[ig][0], myStart);
       int max2 = std::min(smoGroups[ig][1], myEnd);
       // scale gradients
       switch (sclTyp) 
       {
       case 0:
	 break;
       case 1:
         for(int i=min2; i<max2; ++i) 
	   { 
	     orgNorm +=vector[i]*vector[i];
	     newNorm +=work[i-myStart]*work[i-myStart]; 
	   }
#ifdef USE_MPI
	 orgNorm = structCom->globalSum(orgNorm);	 
	 newNorm = structCom->globalSum(newNorm);
#endif	 
         if (newNorm > 0 ) { scale = sqrt(orgNorm/newNorm); }
         break;
       case 2:
         for(int i=min2; i<max2; ++i) 
	   { 
	     orgMax  = std::max(vector[i],orgMax);
	     orgMin  = std::min(vector[i],orgMin);
	     orgSum += vector[i];
	     newMax  = std::max(work[i-myStart],  newMax);
	     newMin  = std::min(work[i-myStart],  newMin);
	     newSum +=work[i-myStart]; 
	   }
#ifdef USE_MPI
	 orgMax = structCom->globalMax(orgMax);
	 orgMin = structCom->globalMin(orgMin);
	 orgSum = structCom->globalSum(orgSum);
	 newMax = structCom->globalMax(newMax);
	 newMin = structCom->globalMin(newMin);
	 newSum = structCom->globalSum(newSum);
#endif

         orgSum /= size;
         newSum /= size;

         if( abs(newSum-newMax) > 0 && xFlag > 0 ) 
         {
           scale = (orgSum - orgMax)/(newSum - newMax);
           shift = orgMax - scale*newMax;
         }
         if( abs(newSum-newMin) > 0 && xFlag < 0 ) 
         {
	   scale = (orgSum - orgMin)/(newSum - newMin);
	   shift = orgMin - scale*newMin;
         }
         break;
       default:
	 assert(0);
	 break;
       }

       for(int i=min2; i<max2; ++i) { vector[i] = scale*work[i-myStart] + shift; }

       filePrint(stderr," ... filtering group %d\n",ig+1);
       filePrint(stderr," ... gradient scaling  by: %e\n",scale);
       filePrint(stderr," ... gradient shifting by: %e\n",shift);
     }

   // bcast vector
#ifdef USE_MPI
   for(int iCpu=0; iCpu<structCom->numCPUs(); ++iCpu)
     {
       int iStart = iCpu*locSize;
       int iEnd   = std::min((iCpu+1)*locSize, size);
       int iSend  = std::max(0, iEnd-iStart);
       if(iStart >= size) { break; }
       MPI_Bcast(vector+iStart,
		 iSend, MPI_DOUBLE, iCpu,
		 *(structCom->getCommunicator()));
     }
#endif
   return;
}

//-----------------------------------------------------------------------------

void Optfilter::initialize(Structopt *structopt,double* _var, int _size)
{
   const double refRadius2 = refRadius*refRadius;

   size      = _size;
   variables = _var;

   Weights   = new double*[size];
   Connected = new int*[size];
   numNeighb = new int[size];

   // check grouping information
   
   if (!numGroups) {
     numGroups       = 1;
     smoGroups       = new int*[1];
     smoGroups[0]    = new int[2];
     smoGroups[0][0] = 0;
     smoGroups[0][1] = size;
   } 
   
   int ierror=0;
   
   if ( smoGroups[0][0] < 0 )              ierror=1;
   if ( smoGroups[numGroups-1][1] > size ) ierror=1;

   if (ierror) {
     fprintf(stderr," ERROR: Filter Grouping information incorrect (min/max)\n");
     exit(-1);
   }

   // get influence areas

   double *totmas     = new double[size]; 
   double **midPoints = new double*[size];

   for(int n1=0;n1<size;n1++) midPoints[n1] = new double[3];

   structopt->getMidPointMass(totmas,midPoints);

   // build filter operator
   boost::scoped_array<double> workWeight (new double[size]);
   boost::scoped_array<int>    workConnect(new int   [size]);

   int myStart, myEnd;
#ifdef USE_MPI
   const int locSize = static_cast<int>(ceil(size/structCom->numCPUs()));
   myStart = structCom->myID()*locSize;
   myEnd   = std::min((structCom->myID()+1)*locSize, size);
#else
   myStart = 0; myEnd = size;
#endif

   for(int n1=myStart;n1<myEnd;++n1) 
     {
       numNeighb[n1]  = 0;
       double scale   = 0;     
       int  icount    = 0;
       
       const Node& nodeA = midPoints[n1];
       // get group     
       int min2 = -1;
       int max2 = -1;       
       for(int ig=0; ig<numGroups; ++ig) 
	 {
	   if(n1 >= smoGroups[ig][0] && n1 < smoGroups[ig][1]) 
	     {
	       min2=smoGroups[ig][0];
	       max2=smoGroups[ig][1];
	       break;
	     }
	 }
     
       // if no group: no smoothing     
       if (min2<0) 
	 {
	   min2=n1;
	   max2=n1+1;
	 }
       for(int n2=min2; n2<max2; ++n2) 
	 {
	   const Node& nodeB = midPoints[n2];
	   const double distance2 = nodeB.distance2(nodeA);
	   if(distance2 <= refRadius2) 
	     {
	       const double distance = sqrt(distance2);
	       workConnect[icount] = n2;	    
	       switch (typ) 
		 {	       
		 case 0:
		   workWeight [icount] = 1.0-distance/refRadius;
		   scale = 1.0;
		   break;
		 case 1:
		 case 2:
		   workWeight [icount] = refRadius-distance;
		   scale              += refRadius-distance;
		 default:
		   assert(0);
		 }
	       ++icount;
	     }
	 }     
       Weights[n1]   = new double[icount];
       Connected[n1] = new int[icount];
       numNeighb[n1] = icount;
       for(int i=0;i<icount;i++)
	 {
	   Weights[n1][i]   = workWeight[i]/scale;
	   Connected[n1][i] = workConnect[i];
	 }
     }
#ifdef USE_MPI
   // exchange numNeighb, Weights, Connected between subdomains
   for(int iCpu=0; iCpu<structCom->numCPUs(); ++iCpu)
     {
       int iStart = iCpu*locSize;
       int iEnd   = std::min((iCpu+1)*locSize, size);
       int iSend  = std::max(0, iEnd-iStart);
       if(iStart >= size) { break; }
       MPI_Bcast(numNeighb+iStart,
		 iSend, MPI_INT, iCpu,
		 *(structCom->getCommunicator()));
     }
   for(int iCpu=0; iCpu<structCom->numCPUs(); ++iCpu)
     {
       int iStart = iCpu*locSize;
       int iEnd   = std::min((iCpu+1)*locSize, size);
       if(iStart >= size) { break; }
       for(int n1=iStart; n1<iEnd; ++n1)
	 {
	   if(structCom->myID() != iCpu)
	     {
	       // allocate memory
	       Weights[n1]   = new double[numNeighb[n1]];
	       Connected[n1] = new int[numNeighb[n1]];
	     }
	   MPI_Bcast(Weights[n1],
		     numNeighb[n1], MPI_DOUBLE, iCpu,
		     *(structCom->getCommunicator()));	   
	   MPI_Bcast(Connected[n1],
		     numNeighb[n1], MPI_INT, iCpu,
		     *(structCom->getCommunicator()));	   
	 }
     }
#endif   
   for(int n1=0; n1<size; ++n1) delete [] midPoints[n1];   
   delete [] totmas;
   delete [] midPoints;
}

#endif
