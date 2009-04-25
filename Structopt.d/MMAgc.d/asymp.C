#ifdef STRUCTOPT

#include "mma.h"

void MMAgc:: asymp(int iter)
{

  //  update assymptotes and consistency variables
  
  int i,j;
  
  // update raa
  
  if ( iter > 1.5 ) {
    raa0 = max(raa0eps,(0.1*raa0));
    for (j=0;j<numcon;j++) 
      raa[j]  = max(raa0eps,0.1*raa[j]);
  }  

  // update assymptotes
   
  if ( iter < 2.5 ) {
    for(i=0;i<numvar;i++) {
      low[i] = xval[i] - saFac*(xmax[i]-xmin[i]);
      upp[i] = xval[i] + saFac*(xmax[i]-xmin[i]);
    }
  }
  else {
    for(i=0;i<numvar;i++) {
  
      double factor = 1.0;
      double xxx    = (xval[i]-xold1[i])*(xold1[i]-xold2[i]);
  
      if ( xxx > 0 )  factor = scFac;
      if ( xxx < 0 )  factor = sbFac;
  
      low[i] = xval[i] - factor*(xold1[i] - low[i]);
      upp[i] = xval[i] + factor*(upp[i]   - xold1[i]);
  
      double lowmin = xval[i] - 10.0*(xmax[i]-xmin[i]);
      double lowmax = xval[i] - 0.01*(xmax[i]-xmin[i]);
      double uppmin = xval[i] + 0.01*(xmax[i]-xmin[i]);
      double uppmax = xval[i] + 10.0*(xmax[i]-xmin[i]);
  
      low[i] = max(low[i],lowmin);
      low[i] = min(low[i],lowmax);
      upp[i] = min(upp[i],uppmax);
      upp[i] = max(upp[i],uppmin);
    }
  }
} 

#endif
