#ifndef _GEOM_STATE_H_
#define _GEOM_STATE_H_

class DofSetArray;
class CoordSet;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
class ControlLawInfo;
class BCond;
class Node;


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


class GeomState {
  public:
     NodeState *ns;     // node state (x,y,z position and rotation tensor)
  protected:
     int numnodes;	// number of nodes
     int (*loc)[6];	// dof location array
     double refCG[3];   // reference CG
     double gRot[3][3]; // Global Rotation Matrix
     const CoordSet &X0;
     int    numReal;    // number of 'real' nodes
     bool  *flag; 	// signifies if node is connected to element
   public:

     // Default Constructor
     GeomState();
     GeomState(CoordSet &cs);

     // Constructor
     GeomState(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs);

     // Copy Constructor
     GeomState(const GeomState &);

     virtual ~GeomState();

     NodeState & operator[](int i)  { return ns[i]; }
     const NodeState & operator[](int i) const { return ns[i]; }
     NodeState *getNodeState() { return ns; }

     // int getLocation(int inode, int dof) { return (loc[inode][dof]-1); }
     int numNodes() { return numnodes; }

     void getPositions(double *positions);
     void getRotations(double *rotations);

     void setPositions(double *positions);
     void setRotations(double *rotations);

     void extract(double *p);

     virtual void update(const Vector &);
     virtual void explicitUpdate(CoordSet &cs, const Vector &v);
     virtual void setVelocity(const Vector &, const Vector &, const Vector &);
     virtual void updatePrescribedDisplacement(BCond *dbc, int numDirichlet, 
                                       double delta = 1.0);
     void updatePrescribedDisplacement(double *userDefinedDisplacement,
                                       ControlLawInfo* claw,
                                       CoordSet &cs );
     virtual void midpoint_step_update(Vector &veloc_n, Vector &accel_n, double delta, GeomState &ss,
                                       double beta, double gamma, double alphaf, double alpham);
     virtual void get_inc_displacement(Vector &inc_Vec, GeomState &ss, bool zeroRot = true);
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
};

#endif
