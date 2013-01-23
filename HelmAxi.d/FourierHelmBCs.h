#ifndef _FOURIER_HELMBCS_H_
#define _FOURIER_HELMBCS_H_

#include <Utils.d/resize_array.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>

struct ComplexFDBC {
// Complex Fourier BC of type c*e^(i*k*(x*dx+y*dy+z*dz))
  int nnum;
};

inline ComplexFDBC
FDBC(int n)
 { ComplexFDBC r; 
   r.nnum = n; 
   return r;
 }

struct ComplexFNBC {
// Complex Fourier BC of type c*e^(i*k*(x*dx+y*dy+z*dz))
  int nnum1, nnum2, nnum3;
};

inline ComplexFNBC
FNBC(int n1, int n2)
 { ComplexFNBC r;
    r.nnum1 = n1; r.nnum2 = n2; r.nnum3 = -1;
   return r;
 }

inline ComplexFNBC
FNBC(int n1, int n2, int n3)
 { ComplexFNBC r;
    r.nnum1 = n1; r.nnum2 = n2; r.nnum3 = n3;
   return r;
 }

class FourierHelmBCs {

   public : 

   DComplex cstant;
   double dirX, dirY, dirZ;

   int numDir;
   ResizeArray<ComplexFDBC> dirBC;

   int numNeu;
   ResizeArray<ComplexFNBC> neuBC;

   int numSomm;
   int somType;
   double surfR, surfZ;
   ResizeArray<SommerElement *> somBC;

   int numModes;
   int numSlices;

   FourierHelmBCs();
   void addDirichlet(ComplexFDBC &boundary);
   void addNeuman(ComplexFNBC &boundary);
   void addSommer(SommerElement *ele);
   void setModes(int mode);
   void setSlices(int slice);
   void setSurf(double aR, double aZ);
   void setSomType(int type);
   void setDir(double dx, double dy, double dz);
   void setConst(DComplex C);
   DComplex Dir2Four(int mode,int pos,double k,double x,double y);
   DComplex Neu2Four(int mode,int pos,double k,double x,double y, 
                     double nr, double nz);
   void IntNeu2(int mode, int pos, double k, double *x, double *y, 
               double nx, double ny, DComplex *T);
   void IntNeu3(int mode, int pos, double k, double *x, double *y, 
               double nx, double ny, DComplex *T);
  
};

extern FourierHelmBCs *fourHelmBC;

#endif
