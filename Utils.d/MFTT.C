#include <Utils.d/MFTT.h>
#include <iostream>

MFTTData::MFTTData()
{
 maxval = 16;
 np     = 0;
 curp   = 0;
 time   = new double[maxval];
 value  = new double[maxval];
}

MFTTData::MFTTData(int _id)
{
 maxval = 16;
 np     = 0;
 curp   = 0;
 time   = new double[maxval];
 value  = new double[maxval];
 id = _id;
}

void
MFTTData::add(double t, double v)
{
 if(np == maxval) {
   int n_maxval = 2*maxval;
   double *n_time  = new double[n_maxval];
   double *n_value = new double[n_maxval];
   int i;
   for(i=0; i < maxval; ++i) {
     n_time[i] = time[i];
     n_value[i] = value[i];
   }
  delete[] time;
  delete[] value;
  time   = n_time;
  value  = n_value;
  maxval = n_maxval;
 }
 time[np] = t;
 value[np] = v;
 np++;
}

double
MFTTData::getVal(double t)
{
 // This function returns zero if t < t_min or t > t_max
 
 // np = total number of points
 if (np) {

 // interpolate and make sure we deal with the case 
 // curp = current point
 // curp==np-1 and curp=-1 correctly

 if(time[curp] > t) // Reverse the search
   while(curp >= 0 && time[curp] > t) curp--;
 else
   while(curp < np-1 && time[curp+1] < t) curp++;

   if (curp < 0 )
      return 0;
   else if (curp == np - 1)
      return 0;
   else {
      double v1 = value[curp], v2 = value[curp+1];
      double t1 =  time[curp], t2 =  time[curp+1];
      return  v1 + (v2-v1) / (t2 - t1) * ( t - t1);
   }
 }
 else
   return 1;
}

double
MFTTData::getValAlt(double t)
{
 // This function returns the value corresponding to t_min if t < t_min,
 // and the value corresponding to t_max if t > t_max

 // np = total number of points

 if (np) {
   curp = np/2; // starting point of search

   // interpolate and make sure we deal with the case
   // curp = current point
   // curp==np-1 and curp=-1 correctly

   if(time[curp] > t) // Reverse the search
     while(curp >= 0 && time[curp] > t) curp--;
   else
     while(curp < np-1 && time[curp+1] < t) curp++;

   if(curp < 0 ) 
     return value[0];
   else if (curp == np - 1) 
     return value[curp];
   else {
     double v1 = value[curp], v2 = value[curp+1];
     double t1 =  time[curp], t2 =  time[curp+1];
     if(t2 == t1) return v1;
     else return  v1 + (v2-v1) / (t2 - t1) * ( t - t1);
   }
 }
 else 
   return 0.0;
}

void 
MFTTData::getValAndSlopeAlt(double t, double *v, double *s)
{
 // This function returns the value and slope corresponding to t_min if t < t_min,
 // and the value and slope corresponding to t_max if t > t_max

 // np = total number of points

 if (np) {
   curp = np/2; // starting point of search

   // interpolate and make sure we deal with the case
   // curp = current point
   // curp==np-1 and curp=-1 correctly

   if(time[curp] > t) // Reverse the search
     while(curp >= 0 && time[curp] > t) curp--;
   else
     while(curp < np-1 && time[curp+1] < t) curp++;

   if(curp < 0 ) {
     *v = value[0];
     *s = 0.0;
   }
   else if (curp == np - 1) {
     *v = value[curp];
     *s = 0.0;
   }
   else {
     double v1 = value[curp], v2 = value[curp+1];
     double t1 =  time[curp], t2 =  time[curp+1];
     *v = v1 + (v2-v1) / (t2 - t1) * ( t - t1);
     *s = (v2-v1) / (t2 - t1);
   }
 }
 else { *v = 0.0; *s = 0.0; }
}

