/* interface for feti03. evolved from feti98_interface

   This interface is written assuming one subdomain per processor. In the
   future it will be expanded to cover multiple subdomains per processor.

   Subdomain Variables Passed:
     Node X,Y,Z
     Element Connectivity
     A in sparse format (CSR format, real and imaginary parts passed separately)
     Force Vector (passed in solve)
     list of boundary conditions (Node #, nodal_dof)
     List of Node communicators
     dofs per node
     local to global node map
     resultant displacement (returned in solve)
     parameters and flags

   Processor information  (for future use, but required)
     pointer into subdomains owned on processor
     subdomain number list on this processor
     (for one subdomain/proc, these are {0,1,2,3..n} and {0,1,2,..n-1} )


   FetiH will be used to solve the global equation
     A * displacement = force
     where all variables are complex.


   Calling FetiH is performed in three calls
     CU_Feti *solver=new CU_Feti(...);   // constructor

     solver->Solve( force_vect_real, force_vect_imag, disp_r, disp_imag)
     delete solver;
   Of course, the solver may be called multiple times with different rhs.

   Because the LHS matrices typically change many times in an analysis,
   we add the call,
     solver->ResetLHS(...)


Example

  12--------13--------14--------15--------16--------17
  |         |          |         |         |         |
  |         |          |         |         |         |
  |   0     |     1    |   2     |  3      |  4      |
  |         |          |         |         |         |
  |         |          |         |         |         |
  6---------7----------8---------9---------10-------11
  |         |          |         |         |         |
  |         |          |         |         |         |

  |   5     |     6    |   7     |  8      |  9      |
  |         |          |         |         |         |
  0---------1----------2---------3---------4---------5

  assume 4 subdomains.
    0 - contains elements 0 1 2
    1 - contains elements 3 4
    2 - contains elements 5 6 7
    3 - contains elements 8 9

  assume 2 processors
    0 - contains subdomains 0,1  (elems 0:4)
    1 - contains subdomains 2,3  (elems 5:9)

............................................
Example Using Global numbering
      proc=0, subdomain=0               proc=0, sub=1

  12--------13--------14--------15  15--------16--------17
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  |   0     |     1    |   2     |   |  3      |  4      |
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  6---------7----------8---------9   9---------10-------11

      proc=1, subdomain=2               proc=1, sub=3
  6---------7----------8---------9   9---------10-------11
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  |   5     |     6    |   7     |   |  8      |  9      |
  |         |          |         |   |         |         |
  0---------1----------2---------3   3---------4---------5

............................................
Example Using Local numbering
      proc=0, subdomain=0               proc=0, sub=1

  4---------5----------6---------7   3---------4---------5
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  |   0     |     1    |   2     |   |    0    |   1     |
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  0---------1----------2---------3   0---------1---------2

      proc=1, subdomain=0               proc=1, sub=1
  4---------5----------6---------7   3---------4---------5
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  |   0     |     1    |   2     |   |    0    |   1     |
  |         |          |         |   |         |         |
  0---------1----------2---------3   0---------1---------2

...................

**** in what follows, only the data for subdomain 0 are shown

 Data for subdomain 0
  (note: this contains global nodes 6 7 8 9 12 13 14 15)
  X[]={ 0.0 1.0 2.0 3.0 0.0 1.0 2.0 3.0 }
  Y[]={ 1.0 1.0 1.0 1.0 2.0 2.0 2.0 2.0 }

 Element Connectivity:
  This is determined using local node numbering. Two integer
  vectors are passed. elemptr points to the offset in the 
  elem_conn array. For this example,
  int elem_ptr[]={0 4 8 12 }
  int elem_conn[]={ 0 1 5 4  1 2 6 5   2 3 7 6 }
  The local connectivity is illustrated below
      4------5------6------7
      |      |      |      |
      |      |      |      |
      |      |      |      |
      0------1------2------3

 The complex LHS matrix, A, will be represented in harwell boeing
  format, in compressed sparse row. It is not assumed to be symmetric.
  Two separated matrices are passed with the real and imaginary
  parts included. Each matrix has three vectors, A_ai, A_aj and A_aa.
  Fortran indexing is used in A_ai and A_aj. 

  A_ai[nrows]=nnz+1
  int nrows, A_ai[nrows+1], A_aj[nnz];
  double A_aa[nnz];
  int neq=nrows;

 Force Vector
  the right hand side vector corresponding to the stiffness matrix, and
  dimensioned to neq. The force vector will be provided to feti
  properly averaged at interface nodes, i.e. if an interface degree of 
  freedom is shared by two subdomains, then the force on the interface 
  will be divided by two.
  The force can be thought of as composed of body forces (such as 
  those applied by pressures or gravity), and of applied forces. The
  applied forces must be divided by the number of nodes sharing the
  force so that the total applied force is correct. Body forces need
  not be so divided because the total body force is the sum of body
  forces within the subdomains.
  The force vector is provided in real and imaginary component vectors.

 A dof map for A and F must also be provided. The standard numbering
  would have node 0 having the first degrees of freedom in the
  map (as determined by dofs_per_node below). The dimension of
  the map is neq. 

 boundary conditions:
  int ndofs_deleted
  int bc_node_ids[]
  unsigned short bc_node_dof[]
  both vectors above are dimensioned at "ndofs_deleted". The first
  lists the nodes where constraints are applied. The second, the

  nodal dof (a number between 0:5) where the bc is applied. In our
  example, if the left side were clamped,
    ndofs_deleted = 12

    bc_node_ids[]={ 0 0 0 0 0 0 4 4 4 4 4 4 } (local numbers used)
    bc_node_dof[]={ 0 1 2 3 4 5 0 1 2 3 4 5 }
  NOTE: there are two rigid body mode methods (algebraic and geometric).
  The geometric method is preferred because of robustness in identification
  of rigid body modes. The geometric method requires that feti apply the
  boundary conditions.
  The algebraic method is less stable since the rigid body modes
  must be determined from zero pivots in the stiffness matrix. If the
  algebraic method is used, the rows/columns associated with boundary
  conditions must be removed from the the stiffness matrix. If
  that is done then ndofs_deleted=0 and the other vectors are empty.

  A summary of the boundary condition options and the sizes of associated
  data is shown below. In this table, neq is the number of equations to
  solve. Nbc is the number of boundary conditions.

  RBM Method           size(A) ndofs_deleted
  ------------------------------------------
  algebraic            neq       0
  geometric            neq       nbc



 Node Communicators
    int numads      :       Number of adjacent subdomains
 
    int*           Size numads
    subdad          :       Subdomain adjacent Lists

    int*           Size numads + 1
    com_ptr          :       Pointer array to com_lis
                             com_ptr[numads] marks the end of com_lis

    int*           Size com_ptr[numads + 1]
    com_lis          :       Shared nodes definition. A list of
                             global node numbers that are shared
                             with the neighbors

   Once again take the subdomain 0 as an example:
    numads = 3
    subdad = [1,2,3]
    com_ptr = [0,2,6,7]
    com_lis = [9,15,6,7,8,9,9]


 The dofs/node are determined using a bitwise OR of the active dofs
  in the unconstrained element. For example, if only X translation is
  available "1" is the value. "2" for Y, "4" for Z, 8 for RotX, etc. The
  map may be modified in feti.
  dofmap_on_node[]={63 63 63 63 63 63 63 63} 
      (because these are shells, with 6 dofs/node and 2^6=64)

 local to global_map[]={ 6 7 8 9 12 13 14 15 }


 subdomains/processor
   int subs_ptr_proc[nproc+1]={ 0, 2, 4 }

   int  subs_per_proc[]= { 0, 1,     2, 3 }

   int local_subdomain_index
 NOTE: this subdomain information is for future use. Currently the interface
 requires one subdomain/processor. In that case
    subs_ptr_proc={0, 1, 2, .. nproc}
    subs_per_proc={0, 1, 2, .. nproc-1}
    local_subdomain_index=0

 Computed result
   double displacement[neq]. returned in call to solve.

 Parameters and Flags
   Because the parameters and flags are subject to change, and to provide
   default values for them, a C++ object containing the class is passed.
   The object is included below. It is expected that this class will change
   in the future, but it is also expected that all such changes can be
   included in this interface file.

 The solver returns:
   0 if successful
   positive integers for warnings
   negative integers for failure


*/
#ifndef SALINAS_FETIH_INTERFACE
#define SALINAS_FETIH_INTERFACE



// provide definition of MPI_Comm
#ifdef SALINAS_MPI
#include <mpi.h>
#else

#ifndef MPI_Comm
#define MPI_Comm unsigned int
#endif

#endif


class FetiParams;
class FetiValues;
class MpcLocal;  // defined externally
template<class Scalar> class GenSandiaDomain;
template<class Scalar> class GenDistrVector;
template<class Scalar> class GenFetiDPSolver;

enum face_type {QUAD4_FACE, QUAD8_FACE, TRIA3_FACE, TRIA6_FACE,
                TRIA4_FACE, LINE2_FACE, POINT_FACE};

template <class Scalar> 
class CU_Feti  {

private: // declare, but don't implement
  CU_Feti();
  CU_Feti(const CU_Feti&);
  CU_Feti& operator=(const CU_Feti&);

  GenSandiaDomain<Scalar> *sdomain;
  GenDistrVector<Scalar> *distF, *distD;
  GenFetiDPSolver<Scalar> *fSolver;

public:
  CU_Feti(
    const double *X,
    const double *Y,
    const double *Z,
    int nnodes_this_subdomain,
    int *elem_ptr,
    int *elem_conn,
    int num_elem_this_subdomain,

    int neq,
    const int *map,

    int ndofs_deleted,
    const int *bc_node_ids,
    int *bc_node_dof,

    MPI_Comm *salinas_communicator,
    int numads,
    int* subdad,
    int* com_ptr,
    int* com_lis,

    unsigned short *dofmap_on_node, 
    int *local_to_global_node_map,  

    FetiParams* params,

    int num_threads = 0
  );

  ~CU_Feti();
             
  int Solve(
    Scalar *Force,
    Scalar *U_real,
    FetiValues &returnvals
  );

  int Factor();  // only needs to be called by salinas for eigen

  void zeroLHS();
  int addToLHS(
    const Scalar multiplier,
    const double *A_aa,   // LHS as real
    const int *A_aj,
    const int *A_ai,
    const bool isK = true   // isK = false for M and C terms in direct frf analysis
  );   // returns 0 if ok

  void setMPCs(
    int NumberMpc,   // number of MPC equations
    MpcLocal* MpcVector
  );

  // the following two functions are provided to permit passing the null space
  // back from feti to Salinas.
  int getNumRBMs();           // returns the number of rbms
  void getRBMs(Scalar *rbms); // fills the rbm vectors. memory must
                              // be allocated before calling (neq)

  // Additions for support of retrieval of MPC information
  // The forces of constraint may be retrieved for each constraint equation.
  // Note that the number of MPC equations has already been passed into
  // FETI. It is redundant (but not a conflict) to provide a function to
  // retrieve that number. Note also that NumberMpc is the number of local
  // constraint equations.

  int NumberMpcConstraints();
  // returns the number of local constraint equations or zero if there are 
  // none. Note: this is the total number of constraints, and includes
  // both real and imaginary constraints

  void MpcForces(double *lambdas);
  // this function returns the constraint forces. It returns the values
  // determined after the last solve. If Solve() has not been
  // called previously, it is undefined.
  // Memory management for lambdas is performed in the calling routine. The
  // dimension of lambdas is NumberMpcConstraints().

  MpcLocal* GetLocalMpcs();
  int  updateLocalMpcs(MpcLocal **_mpc);

  void ConstraintProduct(
    int num_vect, 
    const double* R[], 
    Scalar** V,     
    int trans=1
  );
  // with trans nonzero,
  //   this function determines the product Q'*R 
  //   where Q is the constraint matrix, and
  //   R is a pointer to num_vect vectors, each of length neq (the dim of A).
  //   Results are stored in the vectors V. The data will be stored in V
  //   so that the first neq terms are the first vector, etc.
  //   memory for V must be managed in the calling routine. (num_vect*neq)
  // with trans == 0
  //   this function determines the product Q*R 
  //   where Q is the constraint matrix, and
  //   R is a pointer to num_vect vectors, each of length NumberMpcConstraints()
  //   Results are stored in the vectors V. The data will be stored in V
  //   so that the first NumberMpcConstraints()q terms are the first vector, etc.
  //   memory for V must be managed in the calling routine.
  //
  // If no Solve has been performed, the results of this operation are
  // not defined.

/* The following represents the proposed interface for contact.
   Four calls are represented.
     1. setContactSurfaces
     2. setContactNormal( vector of normals )
     3. setContactGap( vector of Gaps )
     4. updateContactDisp( array of displacements)
     5. getContactForce( array of forces )

  One must first call setContactSurfaces. Calls to the other routines
  without first calling this routine are an error. This call
  initializes the contact data. It also sets the normals and initial
  gaps. The initial displacement is set to zero. A call to Solve() may
  directly follow the setContactSurfaces, in other words, the contact
  data is complete and set up once this call has been made.

  setContactNormal may be called anytime after data is initiated with
  setContactSurfaces. setContactNormal resets the normal vectors for
  contact faces.

  setContactGap sets the initial gap. It may be called any time after
  setContactSurfaces has established the contact data structures. 
  Typically it would be called at the start of each time step for a
  nonlinear dynamics analysis.

  updateContactDisp() is used to provide an update of the current
  displacement. Typically it would be called in each iteration of a
  newton iteration. The data passed is the displacement vector on
  every degree of freedom in the model. It represents the total
  displacement incurred since the last call to setContactGap was
  called.

  PSEUDO EXAMPLE:
   in the main routines we make a call like the following.

   setContactSurfaces(...);
   for (int tstep=0; tstep<maxstep; tstep++){
      NewtonLoop(tstep);
      ... etc.
      if ( normals_have_changed)
         setContactNormal(...);
   }

within the Newton loop...

      double du[neq];
      double disp_update[neq];
      zero(disp_update,neq);
      setContactGap(...);
      while (! converged ){
          update_tangent_matrices();
	  update_residual();
	  updateContactDisp(disp_update);
	  
	  Solve(internal_residual,du);
	  for (int j=0;j<neq;j++)
	     disp_update[j] += du[j];
      }

*/



//--------------------------------------------
  void setContactSurfaces(
     // this function is called for each subdomain to set all the
     // potential contact surfaces.
     int numContactNeigh,            // number of neighbouring subdomains
                                     // sharing potential contact interfaces
                                     // with this_subdomain
     int *contactNeighb,             // subdomains with potential contact
                                     // interfaces to this_subdomain
     int *neighbPtr,                 // indexes (of fytpe and face_ptr) that
                                     // represent the start of data for
                                     // each contactNeighb. dimension is
                                     // numContactNeighb+1, last value is
                                     // the total number_of_faces
     enum face_type *ftype,          // face types of potential contact interfaces
     int *face_ptr,                  // indexes (of face_conn) that represent
                                     // the start of data for each face. dimension
                                     // is number_of_faces + 1, last value is
                                     // the cumulative number of interface nodes
                                     // on all faces
     int *face_conn,                 // nodes on each face
     double (*face_normal)[3],       // face normals
     double (*gap_vector)[3],        // inital gap vectors.
     int *glbl_id,                   // unique global contact numbering
     double *normal_tol,             // PJSA 9-10-2007 normal tolerance for ACME search (dimension of distance)
     double *tangent_tol             // PJSA 9-10-2007 tangential tolerance for ACME search (dimension of distance)
     // note: both subdomains sharing a surface should have the
     // same tolerances for that surface
  );

//--------------------------------------------
   int setContactNormal(
      // called on each subdomain to set the normal of each contact face.
      // returns 0 if no error.
	 const double (*face_normal)[3]   // face normals
	 // the vector is of length number_of_faces=neighbPtr[numContactNeighb]
    );

//--------------------------------------------
   int setContactGap( 
      // called on each subdomain to set the contact gap vector. Note
      // that gaps are really scalars. It is assumed that a dot
      // product will be performed inside the solver.
      // returns 0 for no error.

         //const double (*gap_vectors)[3]   // vector of gaps
         const double *disp
	 // the gap_vector is of dimension number of nodes in contact. 
    );

//--------------------------------------------
    int updateContactDisp(
      // called on each subdomain to set the contact displacement used
      // in the calculation of contact gaps.
      // returns 0 for no error.

	  const double *disp_updated
	  // disp_updated is of dimension neq (the number of dofs on
	  // the subdomain). It is the displacement updated since the
	  // last call to setContactGap().
     );


//--------------------------------------------
     int getContactForce(
	  // called on each subdomain to get the contact forces
	  // applied to each degree of freedom.
	  // returns 0 for no error

	   double *contact_forces
	   // contact_forces is of dimension neq (on each
	   // subdomain). It is the force computed from the lagrange
	   // multipliers computed at the previous solve.
      );

//--------------------------------------------
  void setTiedSurfaces(
     // this function is called for each subdomain to set all the
     // `tied' surfaces
     int numTiedNeigh,            // number of neighbouring subdomains
                                  // sharing potential tied surfaces
                                  // with this_subdomain
     int *tiedNeighb,             // subdomains with tied
                                  // interfaces to this_subdomain
     int *neighbPtr,              // indexes (of fytpe and face_ptr) that
                                  // represent the start of data for
                                  // each contactNeighb. dimension is
                                  // numTiedNeighb+1, last value is
                                  // the total number_of_faces
     enum face_type *faceEl_type, // face types of potential contact interfaces
     int *faceEl_ptr,             // indexes (of face_conn) that represent
                                  // the start of data for each face. dimension
                                  // is number_of_faces + 1, last value is
                                  // the cumulative number of interface nodes
                                  // on all faces
     int *faceEl_conn,            // nodes on each face
     double (*faceEl_normal)[3],  // face normals
     double (*gap_vector)[3],     // gap vector
     int *glbl_id,                // unique global surface numbering
     double *normal_tol,          // normal tolerance for ACME search (dimension of distance)
     double *tangent_tol          // tangential tolerance for ACME search (dimension of distance)
     // note: both subdomains sharing a surface should have the
     // same tolerances for that surface
);


};  // end of class definition 

#ifdef _TEMPLATE_FIX_
  #include <Interface.d/CU_Feti.C>
#endif

#endif

/* history

   date       changes
   8/3/05     significant changes to the contact interface with minor
              changes to the tied surface interface.

*/
