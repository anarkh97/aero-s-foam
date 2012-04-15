#ifndef _GEOM_STATE_H_
#define _GEOM_STATE_H_

#include <map>
#include <vector>

class DofSetArray;
class CoordSet;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
class ControlLawInfo;
class BCond;
class Node;
class Elemset;

class NodeState {
//  private:
//    DirectionCosineMatrix* dRdr;        // vector of the partial derivatives of R wrt thetax, thetay and thetaz
//    DirectionCosineMatrix* d2Rdr2;      // vector of the second parital derivatives of R wrt thetax, thetay and thetaz
  public:
    double x,   y,  z;			// x,y,z coordinates
    double d[6], v[6], a[6];	        // x,y,z velocities and accelerations
    double R[3][3];      	        // Rotation Tensor
    void operator=(const NodeState &);
    double diff(const Node &un, int dof);
    NodeState() { for(int i = 0; i < 6; ++ i) d[i] = v[i] = a[i] = 0; }
};

class ElemState {
  public:
    int numInternalStates;
    double *internalStates;
    ElemState() : numInternalStates(0), internalStates(0) {}
    ~ElemState() { if(internalStates) delete [] internalStates; }
    void operator=(const ElemState &);
};

class GeomState {
  protected:
     std::vector<NodeState> ns; // node state (x,y,z position and rotation tensor)
     int numnodes;	// number of nodes
     std::vector<std::vector<int> > loc; // dof location array
     double refCG[3];   // reference CG
     double gRot[3][3]; // Global Rotation Matrix
     const CoordSet &X0;
     int    numReal;    // number of 'real' nodes
     std::vector<int> flag; // signifies if node is connected to element
     int numelems;
     ElemState *es;
     std::map<int,int> emap;
     std::map<std::pair<int,int>,int> multiplier_nodes;

   public:
     // Default Constructor
     GeomState();
     GeomState(CoordSet &cs);

     // Constructor
     GeomState(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs, Elemset *elems = 0);

     // Copy Constructor
     GeomState(const GeomState &);

     virtual ~GeomState();

     NodeState & operator[](int i)  { return ns[i]; }
     const NodeState & operator[](int i) const { return ns[i]; }
     void clearMultiplierNodes();
     void resizeLocAndFlag(DofSetArray &cdsa);

     double * getElemState(int glNum) { return (numelems > 0) ? es[emap[glNum]].internalStates : 0; }
     int getNumElemStates(int glNum) { return (numelems > 0) ? es[emap[glNum]].numInternalStates : 0; }
     int getTotalNumElemStates();

     // int getLocation(int inode, int dof) { return (loc[inode][dof]-1); }
     int numNodes() const { return numnodes; }

     void getPositions(double *positions);
     void getRotations(double *rotations);
     void getElemStates(double *elemStates);

     void setPositions(double *positions);
     void setRotations(double *rotations);
     void setElemStates(double *elemStates);

     void extract(double *p);

     virtual void update(const Vector &);
     virtual void explicitUpdate(CoordSet &cs, const Vector &v);
     virtual void setVelocity(const Vector &, const Vector &, const Vector &);
     virtual void updatePrescribedDisplacement(BCond *dbc, int numDirichlet, 
                                       double delta);
     void updatePrescribedDisplacement(BCond *dbc, int numDirichlet,
                                       CoordSet &cs);
     void updatePrescribedDisplacement(double *userDefinedDisplacement,
                                       ControlLawInfo* claw,
                                       CoordSet &cs);
     virtual void midpoint_step_update(Vector &veloc_n, Vector &accel_n, double delta, GeomState &ss,
                                       double beta, double gamma, double alphaf, double alpham,
                                       bool zeroRot);
     virtual void get_inc_displacement(Vector &inc_Vec, GeomState &ss, bool zeroRot);
     virtual void get_tot_displacement(Vector &totVec);
     void zeroRotDofs(Vector &vec);
     void interp(double, const GeomState &, const GeomState &);
     void diff(const GeomState &unp, Vector &un);
     void diff1(const GeomState &un, Vector &vD, int inode);
     
     GeomState &operator=(const GeomState &); // copy geometric data -- Assume similar geomState objects

     void print();
     void printNode(int nodeNumber);

     // these functions are used by the spring corotator
     void computeGlobalRotation();
     void getGlobalRot(double R[3][3]);
     void computeCG(double cg[3]);
     void computeRotMat(double angle[3], double mat[3][3]);
     void solve(double [3][3], double [3]);
     void computeRotGradAndJac(double cg [3], 
			       double  grad[3], double jac[3][3]);
     void rotate(double mat[3][3], double vec[3]);
     void setNewmarkParameters(double _beta, double _gamma, double _alpham, double _alphaf);

     void addMultiplierNode(std::pair<int,int> &lmpc_id, double value);
     double getMultiplier(std::pair<int,int> &lmpc_id);
};

#endif
