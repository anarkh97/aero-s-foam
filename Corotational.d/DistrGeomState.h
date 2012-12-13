#ifndef _DISTR_GEOM_STATE_H_
#define _DISTR_GEOM_STATE_H_

template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
class GeomState;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;

#include <map>
#include <vector>

class DistrGeomState {
   private:
     GeomState **gs;	// pointer to an array of GeomStates
     int numSub;

  public:
     std::map<std::pair<int,int>, double> *mu; // lagrange multipliers for the contact surfaces
     std::vector<double> *lambda;              // lagrange multipliers for all other constraints

     // Constructor
     DistrGeomState(DecDomain *domain);
     // Copy Constructor
     DistrGeomState(const DistrGeomState &);
     // Destructor
     ~DistrGeomState();

     // return the ith subdomain's GeomState
     GeomState* operator[](int i) const { return gs[i]; }
     GeomState* getSubGeomState(int i) { return gs[i]; }

     // Update the GeomStates
     void update(DistrVector &v);
     void setVelocity(DistrVector &, DistrVector &);

// The following functions are necessary to implement NL dynamics and
// the arclength method
     int getTotalNumElemStates();
     void midpoint_step_update(DistrVector &veloc_n, DistrVector &accel_n, double &delta, DistrGeomState &ss,
                               double beta, double gamma, double alphaf, double alpham, bool zeroRot);
     void get_inc_displacement(DistrVector &inc_Vec, DistrGeomState &ss, bool zeroRot);
     void rotateVec(DistrVector &vec);
     void get_tot_displacement(DistrVector &totVec);
     void interp(double, DistrGeomState &, DistrGeomState &);
     void diff(DistrGeomState &unp, DistrVector &un);
     void print() { };

     //HB
     DistrGeomState &operator=(DistrGeomState &unp);     
     void subCopy(int isub, DistrGeomState &unp);

     int getNumSub() const { return numSub; }

  private:
     void subStep_update(int isub, DistrVector &veloc_n, DistrVector &accel_n,
                         double &delta, DistrGeomState &ss,
                         double beta, double gamma, double alphaf, double alpham, bool zeroRot);
     void subInc_get(int isub, DistrVector &inc_Vec, DistrGeomState &ss, bool zeroRot);
     void subTot_get(int isub, DistrVector &totVec);
     void subInterp(int isub, double&, DistrGeomState &, DistrGeomState &);
     void subDiff(int isub, DistrGeomState &unp, DistrVector &un);
     void subUpdate(int isub, DistrVector &v);
     void subSetVelocity(int isub, DistrVector &v, DistrVector &a);
     void makeSubGeomStates(int isub, DecDomain *domain);
     void subCopyConstructor(int isub, const DistrGeomState &g2);
     void subRotateVec(int isub, DistrVector &vec);
};

#endif
