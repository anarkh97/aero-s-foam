#ifndef _FILTEROPT_H_
#define _FILTEROPT_H_

#ifdef STRUCTOPT

#include <Structopt.d/Structopt_sd.h>
#include <Math.d/Vector.h>

class Optfilter {

   private:

     int    size;
     int    typ;            // typ 0: Thesis Maute
                            // typ 1: Sigmund/Petersson wo. variable scaling
                            // typ 2: Sigmund/Petersson original

     int    sclTyp;         // scaling of filtered gradient 
                            // typ 0: no scaling;
                            // typ 1: by L2-norm;
                            // typ 2: by max. absolute value

     double refRadius;

     double counter;
     double maxCount;

     double minExp;
     double maxExp;
     
     int    numGroups;

     int    **  smoGroups;

     int    *   numNeighb;
     int    **  Connected;
     double **  Weights;

     double *   variables;



   public:

    Optfilter(int,int,double&,double&,double&,double&,int,int,int*);

    void initialize(Structopt*,double*,int);
    void project   (double*,int);
    void print     (int*,int,FILE*);

    void stepCounter() { counter++; }
};

#endif

#endif
