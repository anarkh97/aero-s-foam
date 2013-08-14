#ifndef _MATNL_COROTATOR_H_
#define _MATNL_COROTATOR_H_

#include <Corotational.d/Corotator.h>

class MatNLElement;
class NLMaterial;

class MatNLCorotator : public Corotator {

     MatNLElement *ele;
     bool own;

   public:
     MatNLCorotator(MatNLElement *, bool own = true);
     ~MatNLCorotator(); 

     void getStiffAndForce(GeomState &curState, CoordSet &cs,
                           FullSquareMatrix &elk, double *f, double dt, double t)
      { getStiffAndForce((GeomState *) NULL, curState, cs, elk, f, dt, t); }

     void getInternalForce(GeomState &curState, CoordSet &cs,
                           FullSquareMatrix &elk, double *f, double dt, double t)
      { getInternalForce((GeomState *) NULL, curState, cs, elk, f, dt, t); }

     void getStiffAndForce(GeomState *refState, GeomState &curState, CoordSet &cs,
                           FullSquareMatrix &elk, double *f, double dt, double t);

     void getInternalForce(GeomState *refState, GeomState &curState, CoordSet &cs,
                           FullSquareMatrix &elk, double *f, double dt, double t);

     virtual void extractDeformations(GeomState &geomState, CoordSet &cs,
                                     double *vld, int &nlflag);

     void getNLVonMises(Vector& stress, Vector& weight, GeomState &,
                        GeomState *, CoordSet &, int strIndex, int surface = 0,
                        double *ndTemps = 0, double ylayer = 0,
                        double zlayer = 0, int avgnum = 0, int measure = -1);

     void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0);

     int getNumGaussPoints();
};

#endif
