#ifndef _TEMP_STATE_H_
#define _TEMP_STATE_H_

class DofSetArray;
class CoordSet;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
class ControlLawInfo;
class BCond;
class Node;
class NodeState;
#include <Corotational.d/GeomState.h>


class TemperatureState : public GeomState {
   private:
//     int numnodes;	// number of nodes
//     NodeState *ns;	// node state (x,y,z position and rotation tensor)
     int (*loc)[1];	// dof location array
//     double refCG[3];   // reference CG
//     double gRot[3][3]; // Global Rotation Matrix
//     const CoordSet &X0;
//     int    numReal;    // number of 'real' nodes
//     bool  *flag; 	// signifies if node is connected to element

   public:

     // Default Constructor
     TemperatureState();

     // Constructor
     //TemperatureState(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs) : GeomState(dsa, cdsa, cs) { }
     TemperatureState(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs);

     // Copy Constructor
     TemperatureState(const TemperatureState &);

     ~TemperatureState() { }

/*
     NodeState & operator[](int i)  { return ns[i]; }
     NodeState *getNodeState() { return ns; }

     // int getLocation(int inode, int dof) { return (loc[inode][dof]-1); }
     int numNodes() { return numnodes; }

     void getPositions(double *positions);
     void getRotations(double *rotations);

     void setPositions(double *positions);
     void setRotations(double *rotations);
*/
     void update(const Vector &);
     void updatePrescribedDisplacement(BCond *dbc, int numDirichlet, 
                                       double delta = 1.0);
/*
     void updatePrescribedDisplacement(double *userDefinedDisplacement,
                                       ControlLawInfo* claw,
                                       CoordSet &cs );
*/
     void midpoint_step_update(Vector &veloc_n, Vector &accel_n, double delta, GeomState &ss,
                               double beta, double gamma, double alphaf, double alpham);
     void get_inc_displacement(Vector &inc_Vec, GeomState &ss, bool zeroRot = true);
/*
     void zeroRotDofs(Vector &vec);
     void interp(double, const TemperatureState &, const TemperatureState &);
     void diff(const TemperatureState &unp, Vector &un);
     
     TemperatureState &operator=(const TemperatureState &); // copy geometric data -- Assume similar geomState objects

     void print();
     void printNode(int nodeNumber);

     void computeGlobalRotation();
     void getGlobalRot(double R[3][3]);
     void computeCG(double cg[3]);
     void computeRotMat(double angle[3], double mat[3][3]);
     void solve(double [3][3], double [3]);
     void computeRotGradAndJac(double cg [3], 
			       double  grad[3], double jac[3][3]);
     void rotate(double mat[3][3], double vec[3]);
*/
};

#endif
