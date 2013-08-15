#ifndef _MFTT_H_
#define _MFTT_H_

// Mechanical Force Time Table data
// also used for young's modulus vs temperature table
class MFTTData {
 public:
     int     id;
     int     np;  // number of points
     double* time;
     double* value;
     int maxval;
     int curp; // Because calls to getVal will usually be
               // for neighboring values of t
   public:
     MFTTData();
     MFTTData(int);
     void add(double, double);
     double getVal(double);
     double getValAlt(double);
     void getValAndSlopeAlt(double t, double *v, double *s);
     int getID() { return id; }
     int getNumPoints() { return np; }
     double getT(int i) { return time[i]; }
     double getV(int i) { return value[i]; }
};


#endif
