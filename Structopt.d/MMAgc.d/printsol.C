#ifdef STRUCTOPT

#include <stdio.h>

#include "mma.h"

void MMAgc::printsol(int iter)
{
   int i,j;
   
   fprintf(ofile,"\n\n");
   fprintf(ofile,"ITERATION %5d\n",iter);
   fprintf(ofile,"\n");

   fprintf(ofile,"VARIABLES:	  ");
   fprintf(ofile,"VALUE 	   	");
   fprintf(ofile,"LOWER ASSYMPOTE  	");
   fprintf(ofile,"UPPER ASSYMPOTE  	");
   fprintf(ofile,"LOWER LOCAL BOUND	");
   fprintf(ofile,"UPPER LOCAL BOUND\n");

   for (i=0;i<numvar;i++) {
     fprintf(ofile,"VARIABLE %5d  %20.10e  %20.10e  %20.10e  %20.10e  %20.10e\n",
                    i+1,xval[i],low[i],upp[i],alfa[i],beta[i]);
   }

   fprintf(ofile,"\n");
   fprintf(ofile,"OBJECTIVE :         %20.10e   %20.10e\n",f0val,raa0);

   for(j=0;j<numcon;j++) {
     fprintf(ofile,"CONSTRAINT: %5d   %20.10e   %20.10e\n",j+1,-1.0*fval[j],raa[j]);
   }

   fprintf(ofile,"\n");
   
   fflush(ofile);
}

#endif
