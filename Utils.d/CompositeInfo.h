#ifndef _COMPOSITE_INFO_H_
#define _COMPOSITE_INFO_H_

#include <Utils.d/resize_array.h>

// Coefficient Data class for composites
class CoefData {
    public:
      double c[6][6];
    public:
      //CoefData() {};
      //CoefData(double *d);
      void zero();
      void setCoef(int,int,double);
      void setCoef(double d[36]);
      double *values () { return (double *)c; }
};


// composite layer information
class LayInfo {
   public:
     int numLayers;
     int maxLayer;
     int type;
   private:
     ResizeArray<int> matids; // PJSA 3-30-05
   public:
     double (*data)[9];
     double (*grad)[9];

     LayInfo(int _type);
 
     void add(int k, double *d, int m = -1); // PJSA 3-30-05 added m (layer material id)
     void setAllLayers(int, double (*)[9]);
     void setGrad();
     void zeroGrad();

     int nLayers() { return numLayers; }
     int getType() { return type; }

     double *values()  { return (double *) data; }
     double *gradval() { return (double *) grad; }

     void setLayerMaterialProperties(int k, double *d);
     int getLayerMaterialId(int k);
};

// composite layer material information
class LayMat {
  public:
    int m;
    double data[7]; // data[] = { E1  E2  nu12  G12  mu1,12 mu2,12 rho }
    LayMat(int _m, double *d) { m = _m; for(int i=0; i<7; ++i) data[i] = d[i]; }
};

#endif
