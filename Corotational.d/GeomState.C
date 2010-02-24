#include <Driver.d/Domain.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Utils.d/dofset.h>
#include <Element.d/Element.h>
#include <Math.d/SparseMatrix.h>
#include <Driver.d/GeoSource.h>

GeomState::GeomState(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs)
 : X0(cs)
/****************************************************************
 *
 *  Purpose: determine geometric state of nodal coordinates
 *
 *  Input:
 *     DofSetArray : Constrained Degree of freedom set array 
 *     CoordSet    : Coordinate set 
 *
 *  Output:
 *     ns   :  node state updated  
 *
 *  Coded by: Michel Lesoinne and Teymour Manzouri
 ***************************************************************/
{
  numnodes = dsa.numNodes();		// Number of nodes
  ns       = new NodeState[numnodes];	// Array of Node States
  loc      = new int[numnodes][6];	// dof locations	
  flag	   = new bool[numnodes];     // flag for node to element connectivity
  numReal = 0;

  int i;
  for(i=0; i<numnodes; ++i) {

    // Store location of each degree of freedom
    loc[i][0] = cdsa.locate( i, DofSet::Xdisp );
    loc[i][1] = cdsa.locate( i, DofSet::Ydisp );
    loc[i][2] = cdsa.locate( i, DofSet::Zdisp );
    loc[i][3] = cdsa.locate( i, DofSet::Xrot  );
    loc[i][4] = cdsa.locate( i, DofSet::Yrot  );
    loc[i][5] = cdsa.locate( i, DofSet::Zrot  );
 
    // Get Node i from the Coordinate (Node) set
    Node *node_i = cs[i];

    if (node_i)  {

      // Set the ith node's coordinates
      ns[i].x = node_i->x;
      ns[i].y = node_i->y;
      ns[i].z = node_i->z;
 
      // Set ith node's rotation tensor equal to identity
      ns[i].R[0][0] = 1.0;
      ns[i].R[0][1] = 0.0;
      ns[i].R[0][2] = 0.0;
      ns[i].R[1][0] = 0.0;
      ns[i].R[1][1] = 1.0;
      ns[i].R[1][2] = 0.0;
      ns[i].R[2][0] = 0.0;
      ns[i].R[2][1] = 0.0;
      ns[i].R[2][2] = 1.0;

      if(dsa[i].list() != 0)  {
        flag[i] = true;
        numReal++;
      }
      else 
        flag[i] = false;

    }
    else  {
      // HB: TEMPORARY BAD FIX FOR DEALING WITH LAGRANGE MULTIPLIERS (RIGID BAR) !!!
      ns[i].x = 0.0;
      ns[i].y = 0.0;
      ns[i].z = 0.0;
      // Set ith node's rotation tensor equal to identity
      ns[i].R[0][0] = 1.0;
      ns[i].R[0][1] = 0.0;
      ns[i].R[0][2] = 0.0;
      ns[i].R[1][0] = 0.0;
      ns[i].R[1][1] = 1.0;
      ns[i].R[1][2] = 0.0;
      ns[i].R[2][0] = 0.0;
      ns[i].R[2][1] = 0.0;
      ns[i].R[2][2] = 1.0;

      loc[i][0] = cdsa.locate( i, DofSet::Lagrange );
      flag[i] = false;
      //flag[i] = true;
    }
    
  }

  // Initialize Global Rotation Matrix to Identity
  double zeroRot[3] = {0.0, 0.0, 0.0};
  computeRotMat(zeroRot, gRot);
  computeCG(refCG);

}

CoordSet emptyCoord;

GeomState::GeomState() : ns(NULL), numnodes(0), loc(NULL), X0(emptyCoord), numReal(0), flag(NULL) 
{
}

GeomState::GeomState(CoordSet &cs) : ns(NULL), numnodes(0), loc(NULL), X0(cs), numReal(0), flag(NULL) 
{
}

GeomState::~GeomState() {
  delete[] flag;
  delete[] loc;
  delete[] ns;
}

void
GeomState::print()
{
/*
 // Prints nodal coordinates and associated rotation tensor
 int i;
 for(i=0; i<numnodes; ++i) {
   fprintf(stderr,"inode\tx\ty\tz\n");
   fprintf(stderr,"#%d\t%e\t%e\t%e\n",i,ns[i].x,ns[i].y,ns[i].z);
   fprintf(stderr,"Rotation Tensor\n");
   fprintf(stderr,"% e % e % e\n",ns[i].R[0][0],ns[i].R[0][1],ns[i].R[0][2]);
   fprintf(stderr,"% e % e % e\n",ns[i].R[1][0],ns[i].R[1][1],ns[i].R[1][2]);
   fprintf(stderr,"% e % e % e\n",ns[i].R[2][0],ns[i].R[2][1],ns[i].R[2][2]);
 }
*/
 for(int i=0; i<numnodes; ++i) {
   if(loc[i][0] >= 0) cerr << ns[i].x << " ";
   if(loc[i][1] >= 0) cerr << ns[i].y << " ";
   if(loc[i][2] >= 0) cerr << ns[i].z << " ";
 }
 cerr << endl;
}

void
GeomState::printNode(int i)
{
   fprintf(stderr,"inode\tx\ty\tz\n");
   fprintf(stderr,"#%d\t%e\t%e\t%e\n",i,ns[i].x,ns[i].y,ns[i].z);
   fprintf(stderr,"Rotation Tensor\n");
   fprintf(stderr,"% e % e % e\n",ns[i].R[0][0],ns[i].R[0][1],ns[i].R[0][2]);
   fprintf(stderr,"% e % e % e\n",ns[i].R[1][0],ns[i].R[1][1],ns[i].R[1][2]);
   fprintf(stderr,"% e % e % e\n",ns[i].R[2][0],ns[i].R[2][1],ns[i].R[2][2]);
}

GeomState &
GeomState::operator=(const GeomState &g2)
{
  int i,j;
  for(i=0; i<numnodes; ++i)
    ns[i] = g2.ns[i];

  for(i=0; i<3; ++i)
    for(j=0; j<3; ++j)
      gRot[i][j] = g2.gRot[i][j];

  return *this;
}

GeomState::GeomState(const GeomState &g2) : X0(g2.X0)
{
  // Copy number of nodes
  numnodes = g2.numnodes;

  // Allocate memory for node states & dof locations
  ns  = new NodeState[numnodes];
  loc = new int[numnodes][6];

  // flag for node to element connectivity
  flag     = new bool[numnodes];     

  numReal = g2.numReal;

  // Copy dof locations
  int i;
  for(i = 0; i < numnodes; ++i) {
    loc[i][0] = g2.loc[i][0];
    loc[i][1] = g2.loc[i][1];
    loc[i][2] = g2.loc[i][2];
    loc[i][3] = g2.loc[i][3];
    loc[i][4] = g2.loc[i][4];
    loc[i][5] = g2.loc[i][5];
  }

  // Copy node states
  for(i = 0; i < numnodes; ++i) {
    ns[i]  = g2.ns[i];
    flag[i]= g2.flag[i];
  }
 
  // Initialize Global Rotation Matrix & CG position // HB
  refCG[0] = g2.refCG[0];
  refCG[1] = g2.refCG[1];
  refCG[2] = g2.refCG[2];
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      gRot[i][j] = g2.gRot[i][j];
}

void
NodeState::operator=(const NodeState &node)
{
 // Set x, y, and z coordinate values
 this->x = node.x;
 this->y = node.y;
 this->z = node.z;

 // Set rotation tensor
 this->R[0][0] = node.R[0][0];
 this->R[0][1] = node.R[0][1];
 this->R[0][2] = node.R[0][2];
 this->R[1][0] = node.R[1][0];
 this->R[1][1] = node.R[1][1];
 this->R[1][2] = node.R[1][2];
 this->R[2][0] = node.R[2][0];
 this->R[2][1] = node.R[2][1];
 this->R[2][2] = node.R[2][2];
}

void
GeomState::update(const Vector &v)
{
 //if(v*v == 0.0) return; // XXXX
 // v = incremental displacement vector

 double dtheta[3];

 int i;
 for(i=0; i<numnodes; ++i) {

     // Set incremental translational displacements

     double dx = (loc[i][0] >= 0) ? v[loc[i][0]] : 0.0;
     double dy = (loc[i][1] >= 0) ? v[loc[i][1]] : 0.0;
     double dz = (loc[i][2] >= 0) ? v[loc[i][2]] : 0.0;

     // Set incremental rotations

     dtheta[0] = (loc[i][3] >= 0) ? v[loc[i][3]] : 0.0;
     dtheta[1] = (loc[i][4] >= 0) ? v[loc[i][4]] : 0.0;
     dtheta[2] = (loc[i][5] >= 0) ? v[loc[i][5]] : 0.0;

     // Increment total translational displacements

     ns[i].x += dx;
     ns[i].y += dy;
     ns[i].z += dz;

     // Increment rotation tensor R = R(dtheta)Ra
     //if(dtheta[0] == 0.0 && dtheta[1] == 0.0 && dtheta[2] == 0.0) continue; // XXXX
     inc_rottensor( dtheta, ns[i].R );
   }

 computeGlobalRotation();
}

void
GeomState::setVelocity(const Vector &v)
{
 // v = velocity vector

 int i;
 for(i=0; i<numnodes; ++i) {

     // Set Velocity

     double Vx = (loc[i][0] >= 0) ? v[loc[i][0]] : 0.0;
     double Vy = (loc[i][1] >= 0) ? v[loc[i][1]] : 0.0;
     double Vz = (loc[i][2] >= 0) ? v[loc[i][2]] : 0.0;

     ns[i].vx = Vx;
     ns[i].vy = Vy;
     ns[i].vz = Vz;

   }
}

void
GeomState::midpoint_step_update(Vector &vel_n, Vector &acc_n, double delta, GeomState &ss,
                                double beta, double gamma, double alphaf, double alpham)
{
 // note: delta = dt/2
 double coef = 1/(1-alphaf);
 double dcoef = gamma/(2*delta*beta);
 double vcoef = (1-(1-alphaf)*gamma/beta)-alphaf;
 double acoef = 2*delta*(1-alphaf)*(2*beta-gamma)/(2*beta);
 double rcoef = alphaf/(1-alphaf);

 // x,y,z velocity step update (velocity^{n+1} = 1/(1-alphaf)*(velocity^{n+1/2} - alphaf*velocity^n)
 // note: we are computing translational velocity^{n+1/2} locally
 //       as velocity^{n+1/2} = gamma/(dt*beta)*(d^{n+1-alphaf}-d^n) + (1-(1-alphaf)*gamma/beta)*v^n + dt*(1-alphaf)(2*beta-gamma)/(2*beta)*a^n
 for(int i = 0; i < numnodes; ++i) {
   if(loc[i][0] >= 0)
     vel_n[loc[i][0]] = coef*(dcoef*(ns[i].x - ss.ns[i].x) + vcoef*vel_n[loc[i][0]] + acoef*acc_n[loc[i][0]]);

   if(loc[i][1] >= 0)
     vel_n[loc[i][1]] = coef*(dcoef*(ns[i].y - ss.ns[i].y) + vcoef*vel_n[loc[i][1]] + acoef*acc_n[loc[i][1]]);

   if(loc[i][2] >= 0)
     vel_n[loc[i][2]] = coef*(dcoef*(ns[i].z - ss.ns[i].z) + vcoef*vel_n[loc[i][2]] + acoef*acc_n[loc[i][2]]);
 }

  // Update step translational displacements
  for(int i = 0; i < numnodes; ++i) {
    ss.ns[i].x = ns[i].x = coef*(ns[i].x - alphaf*ss.ns[i].x);
    ss.ns[i].y = ns[i].y = coef*(ns[i].y - alphaf*ss.ns[i].y);
    ss.ns[i].z = ns[i].z = coef*(ns[i].z - alphaf*ss.ns[i].z);
  }

  // Update step rotational tensor
  double result[3][3], result2[3][3], rotVec[3];
  for(int i = 0; i < numnodes; ++i) {
    if(alphaf == 0.0) continue; // nothing to do in this case

    mat_mult_mat(ss.ns[i].R, ns[i].R, result2, 1);
    if(alphaf != 0.5) {
      mat_to_vec(result2, rotVec);
      rotVec[0] *= rcoef;
      rotVec[1] *= rcoef;
      rotVec[2] *= rcoef;
      vec_to_mat(rotVec, result2);
    }
    mat_mult_mat(ns[i].R, result2, result, 0);

    for(int j = 0; j < 3; ++j)
      for(int k = 0; k < 3; ++k)
        ss.ns[i].R[j][k] = ns[i].R[j][k] = result[j][k];
  }
}

void
GeomState::interp(double alphaf, const GeomState &gs_n, const GeomState &gs_nplus1)
{
  double coef = 1.0 - alphaf;

  // Update displacements: d_{n+1-alphaf} = (1-alphaf)*d_{n+1} + alphaf*d_n
  int inode;
  for(inode=0; inode<numnodes; ++inode) {
    ns[inode].x = alphaf*gs_n.ns[inode].x + coef*gs_nplus1.ns[inode].x;
    ns[inode].y = alphaf*gs_n.ns[inode].y + coef*gs_nplus1.ns[inode].y;
    ns[inode].z = alphaf*gs_n.ns[inode].z + coef*gs_nplus1.ns[inode].z;
  }

  // Update rotations: note this is done using SLERP, interpolating from R_{n+1} to R_n with the incremental rotation defined as a multiplication from the left:
  // R_n = R(l)*R_{n+1} hence R(l) = R_n*R_{n+1}^T and R_{n+1-alphaf} = R(alphaf*l)*R_{n+1}, where l is a rotation "vector" (axis/angle representation)
  // typically SLERP is done using the incremental rotation defined as a multiplication from the right:
  // R_n = R_{n+1}*R(r) hence R(r) = R_{n+1}^T*R_n and R_{n+1-alphaf} = R_{n+1}*R(alphaf*r)
  double rotMat[3][3], rotVec[3];
  for(inode=0; inode<numnodes; ++inode) {
     mat_mult_mat(gs_n.ns[inode].R, gs_nplus1.ns[inode].R, rotMat, 2);
     mat_to_vec(rotMat, rotVec);
     rotVec[0] *= alphaf;
     rotVec[1] *= alphaf;
     rotVec[2] *= alphaf;
     vec_to_mat(rotVec, rotMat);
     mat_mult_mat(rotMat, gs_nplus1.ns[inode].R, ns[inode].R, 0);
  }
}

void
GeomState::diff(const GeomState &un, Vector &vD)
{
  double vec[3], dR[3][3];

  // Loop over all of nodes
  int inode;
  for(inode=0; inode<numnodes; ++inode) {

    // Perform translational dof difference
    if(loc[inode][0] >= 0) vD[loc[inode][0]] = ns[inode].x-un.ns[inode].x;
    if(loc[inode][1] >= 0) vD[loc[inode][1]] = ns[inode].y-un.ns[inode].y;
    if(loc[inode][2] >= 0) vD[loc[inode][2]] = ns[inode].z-un.ns[inode].z;

    // Perform rotational dof difference
    mat_mult_mat(ns[inode].R, un.ns[inode].R, dR, 2);
    mat_to_vec(dR, vec);

    if(loc[inode][3] >= 0) vD[loc[inode][3]] = vec[0];
    if(loc[inode][4] >= 0) vD[loc[inode][4]] = vec[1];
    if(loc[inode][5] >= 0) vD[loc[inode][5]] = vec[2]; 
  }
}

void
GeomState::diff1(const GeomState &un, Vector &vD, int inode)
{
  double vec[3], dR[3][3];

  // Perform translational dof difference
  if(loc[inode][0] >= 0) vD[0] = ns[inode].x-un.ns[inode].x;
  if(loc[inode][1] >= 0) vD[1] = ns[inode].y-un.ns[inode].y;
  if(loc[inode][2] >= 0) vD[2] = ns[inode].z-un.ns[inode].z;

  // Perform rotational dof difference
  mat_mult_mat(ns[inode].R, un.ns[inode].R, dR, 2);
  mat_to_vec(dR, vec);

  if(loc[inode][3] >= 0) vD[3] = vec[0];
  if(loc[inode][4] >= 0) vD[4] = vec[1];
  if(loc[inode][5] >= 0) vD[5] = vec[2];
}

void
GeomState::get_inc_displacement(Vector &incVec, GeomState &ss, bool zeroRot)
{
  // Update incremental translational displacements and rotations
  int inode;
  for(inode=0; inode<numnodes; ++inode) {

    // Update incremental translational displacements
    if(loc[inode][0] >= 0) incVec[loc[inode][0]] = ns[inode].x - ss.ns[inode].x;
    if(loc[inode][1] >= 0) incVec[loc[inode][1]] = ns[inode].y - ss.ns[inode].y;
    if(loc[inode][2] >= 0) incVec[loc[inode][2]] = ns[inode].z - ss.ns[inode].z;
    
    if(loc[inode][3] >= 0 || loc[inode][4] >= 0 || loc[inode][5] >= 0) {
      if(zeroRot) {
        // Set rotational displacements equal to zero.
        if(loc[inode][3] >= 0) incVec[loc[inode][3]] = 0.0;
        if(loc[inode][4] >= 0) incVec[loc[inode][4]] = 0.0;
        if(loc[inode][5] >= 0) incVec[loc[inode][5]] = 0.0;
      }
      else {
        double R[3][3], vec[3];
        mat_mult_mat( ns[inode].R, ss.ns[inode].R, R, 2 ); // LEFT
        mat_to_vec( R, vec );
        if( loc[inode][3] >= 0 ) incVec[loc[inode][3]] = vec[0];
        if( loc[inode][4] >= 0 ) incVec[loc[inode][4]] = vec[1];
        if( loc[inode][5] >= 0 ) incVec[loc[inode][5]] = vec[2];
      }
    } 
  }

}

void
GeomState::zeroRotDofs(Vector& vec)
{
  for(int inode = 0; inode < numnodes; ++inode) {
  // Set rotational displacements equal to zero.
    if(loc[inode][3] >= 0) vec[loc[inode][3]] = 0.0;
    if(loc[inode][4] >= 0) vec[loc[inode][4]] = 0.0;
    if(loc[inode][5] >= 0) vec[loc[inode][5]] = 0.0;
  }
}

// updating prescribed displacements for nonlinear statics
// i.e. non-zero displacement boundary conditions

void
GeomState::updatePrescribedDisplacement(BCond* dbc, int numDirichlet,
                                        double delta)
{
  // allocate space to store rotational prescribed dofs
  double (*dth)[3] = new double[numnodes][3];

  // initialize to zero, rotational prescribed dofs
  int i;
  for(i=0; i<numnodes; ++i)
    dth[i][0] = dth[i][1] = dth[i][2] = 0.0;

  for(i=0; i<numDirichlet; ++i) {

    int nodeNumber = dbc[i].nnum;
    int dofNumber  = dbc[i].dofnum;
    

    // we multiply the total prescribed value by delta which
    // is a parameter prescribed by the user in the input file
    // it is a load control parameter. By default, it is set to 1.0
    // which is applying the full prescribed displacement all at once.

    double prescribedValue = delta*dbc[i].val;

    // if prescribed value is zero, we do nothing.
    if( prescribedValue == 0.0 ) continue;

    switch(dofNumber) {
        case 0:
                ns[nodeNumber].x += prescribedValue;
                break;
        case 1:
                ns[nodeNumber].y += prescribedValue;
                break;
        case 2:
                ns[nodeNumber].z += prescribedValue;
                break;
        case 3:
                dth[nodeNumber][0] = prescribedValue;
                break;
        case 4:
                dth[nodeNumber][1] = prescribedValue;
                break;
        case 5:
                dth[nodeNumber][2] = prescribedValue;
                break;
        default:
                break;
    }

  }

  // Take care of rotational degrees of freedom
  for(i=0; i<numnodes; ++i) {
    if(dth[i][0] == 0.0 && dth[i][1] == 0.0 && dth[i][2] == 0.0) continue;
    // inc_rottensor( dth[i], ns[i].R );
    form_rottensor( dth[i], ns[i].R );
  }

  delete [] dth;
}

// update prescribed displacements for nonlinear dynamics
// i.e. non-zero displacement boundary conditions prescribed with
// USDD which are time dependent user defined displacements
void
GeomState::updatePrescribedDisplacement(double *v, ControlLawInfo *claw,
                                        CoordSet &cs )
{
  if(claw->numUserDisp == 0) return;
  // allocate space to store rotational prescribed dofs
  double (*dth)[3] = new double[claw->numUserDisp][3];

  // initialize to zero, rotational prescribed dofs
  int i;
  for(i=0; i<claw->numUserDisp; ++i)
    dth[i][0] = dth[i][1] = dth[i][2] = 0.0;

  for(i=0; i<claw->numUserDisp; ++i) {
  
    int nodeNumber = claw->userDisp[i].nnum;
    int dofNumber  = claw->userDisp[i].dofnum;
    if (nodeNumber < 0 || nodeNumber >= numnodes)  {
      fprintf(stderr, "Bad cntrl law for node number %d(%d) dof(%d)\n", nodeNumber+1,numnodes, claw->userDisp[i].dofnum);
    }
    
    switch(dofNumber) {
    	case 0: 
		ns[nodeNumber].x = ( cs[nodeNumber]->x + v[i] );
		break;
	case 1:
		ns[nodeNumber].y = ( cs[nodeNumber]->y + v[i] );
		break;
	case 2:
		ns[nodeNumber].z = ( cs[nodeNumber]->z + v[i] );
		break;
	case 3:
	        dth[i][0] = v[i];
		break;
	case 4:
	        dth[i][1] = v[i];
		break;
	case 5:
	        dth[i][2] = v[i];
		break;
	default:
		break;
    }

  }

  // Take care of rotational degrees of freedom for
  // the prescribed displacements
  for(i=0; i<claw->numUserDisp; ++i) {
    if(dth[i][0] == 0.0 &&
       dth[i][1] == 0.0 &&
       dth[i][2] == 0.0) continue;
    //inc_rottensor( dth[i], ns[claw->userDisp[i].nnum].R );
    form_rottensor( dth[i], ns[claw->userDisp[i].nnum].R );
  }

 // check whether this is the appropriate method of
 // deleting this data structure!
 delete [] dth;
}

/*
void
GeomState::extractR(CoordSet &cs, Vector &rigid)
{
 Node &node1 = cs.getNode(n1);

 int i;
 for(i=0; i<numnodes; ++i) {
   Node &nd = cs.getNode(i);
   rigid[loc[i][0]] = ns[i].x - nd.x;
   rigid[loc[i][0]] = ns[i].y - nd.y;
   rigid[loc[i][0]] = ns[i].z - nd.z;
 }
   
}
*/

void
GeomState::getPositions(double *positions)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   positions[3*i+0] = ns[i].x;
   positions[3*i+1] = ns[i].y;
   positions[3*i+2] = ns[i].z;
 }
}

void
GeomState::getRotations(double *rotations)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   rotations[9*i+0] = ns[i].R[0][0];
   rotations[9*i+1] = ns[i].R[0][1];
   rotations[9*i+2] = ns[i].R[0][2];
   rotations[9*i+3] = ns[i].R[1][0];
   rotations[9*i+4] = ns[i].R[1][1];
   rotations[9*i+5] = ns[i].R[1][2];
   rotations[9*i+6] = ns[i].R[2][0];
   rotations[9*i+7] = ns[i].R[2][1];
   rotations[9*i+8] = ns[i].R[2][2];
 }

}

void
GeomState::setPositions(double *positions)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   ns[i].x = positions[3*i+0];
   ns[i].y = positions[3*i+1];
   ns[i].z = positions[3*i+2];
 }
}

void
GeomState::setRotations(double *rotations)
{
 int i;
 for(i=0; i<numnodes; ++i) {
   ns[i].R[0][0] = rotations[9*i+0];
   ns[i].R[0][1] = rotations[9*i+1];
   ns[i].R[0][2] = rotations[9*i+2];
   ns[i].R[1][0] = rotations[9*i+3];
   ns[i].R[1][1] = rotations[9*i+4];
   ns[i].R[1][2] = rotations[9*i+5];
   ns[i].R[2][0] = rotations[9*i+6];
   ns[i].R[2][1] = rotations[9*i+7];
   ns[i].R[2][2] = rotations[9*i+8];
 }

}

void GeomState::computeGlobalRotation() 
{
  double cg[3];
  computeCG(cg);

  double deltaRot[3][3];  // incremental rotation
  // initialize deltaRot to Identity Matrix
  double zeroRot[3] = {0,0,0};
  computeRotMat(zeroRot, deltaRot);

  int iter = 3;  // number of iterations
  double jac[3][3], grad[3];  //minimization gradient and jacobian

  for (int it = 0; it < iter; it++)  {
    int i,j;
    // compute minimization gradients and jacobians
    computeRotGradAndJac(cg, grad, jac);

    // rotation results come back in grad
    solve(jac, grad);

    // update deltaRot
    computeRotMat(grad, deltaRot);

    // update global rotation matrix
    double R[3][3];
    mat_mult_mat(deltaRot, gRot, R, 0);
  
    for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
        gRot[i][j] = R[i][j];
  }
} 

void GeomState::computeRotMat(double *angle, double mat[3][3])
{
  // trig functions of angles
  double c1 = cos(angle[0]);
  double s1 = sin(angle[0]);
  double c2 = cos(angle[1]);
  double s2 = sin(angle[1]);
  double c3 = cos(angle[2]);
  double s3 = sin(angle[2]);


  /* Compute rotation matrix
     computed as R1.R2.R3
     where R1 is rotation about z
           R2 is rotation about y
           R3 is rotation about x
  mat[0][0] = c1*c2;
  mat[0][1] = c1*s2*s3 - c3*s1;
  mat[0][2] = c1*c3*s2 + s1*s3;
  
  mat[1][0] = c2*s1;
  mat[1][1] = c1*c3+s1*s2*s3;
  mat[1][2] = c3*s1*s2 - c1*s3;
  
  mat[2][0] = -s2;
  mat[2][1] = c2*s3;
  mat[2][2] = c2*c3; */
  

  /* computed as R1.R2.R3
     where R1 is rotation about x
           R2 is rotation about y
           R3 is rotation about z
  */ 
  mat[0][0] = c2*c3;
  mat[0][1] = -c2*s3;
  mat[0][2] = s2;

  mat[1][0] = c3*s1*s2 + c1*s3;
  mat[1][1] = c1*c3 - s1*s2*s3;
  mat[1][2] = -c2*s1;

  mat[2][0] = s1*s3 - c1*c3*s2;
  mat[2][1] = c3*s1 + c1*s2*s3;
  mat[2][2] = c1*c2;
}

void GeomState::solve(double m[3][3], double v[3]) 
{
  int i,j,k;

  for (i = 0; i < 2; i++)
    for (j = i+1; j < 3; ++j) {
      if (m[i][i] == 0.0)  {
        for (k = 0; k < 3; k++)  m[i][k] = 0;
        m[i][i] = 1.0;
        v[i] = 0;
      }
      double coef = m[j][i]/m[i][i];
      for (k = i+1; k < 3; ++k)
        m[j][k] -= coef*m[i][k];
      v[j] -= coef*v[i];
    }
//NOT VERY CLEAN
  if (m[2][2] == 0.0)  {
    for (k = 0; k < 2; k++)  m[2][k] = 0;
      m[2][2] = 1.0;
      v[2] = 0;
    }

  for (i=2; i >= 0; i--) {
    for (j=2; j > i; j--)
      v[i] -= m[i][j] * v[j];
    v[i] /= m[i][i];

  }
}

void GeomState::computeRotGradAndJac(double cg[3], double grad[3], double jac[3][3]) 
{
  // init grad and jac
  int i,j;
  for (i = 0; i < 3; i++)  {
    grad[i] = 0.0;
    for (j = 0; j < 3; j++)
      jac[i][j] = 0.0;
  }

  // rotate local vectors using R(n-1)
  for (i = 0; i < numnodes; i++) {
    if (!flag[i]) continue;

    double rd[3];
    rd[0] = X0[i]->x - refCG[0]; 
    rd[1] = X0[i]->y - refCG[1]; 
    rd[2] = X0[i]->z - refCG[2]; 
    rotate(gRot, rd);

    // compute freq. used values
    double eVec[3];  // Deformed X - rotated loc. vec - deformed cg
    eVec[0] = ns[i].x - cg[0] - rd[0];
    eVec[1] = ns[i].y - cg[1] - rd[1];
    eVec[2] = ns[i].z - cg[2] - rd[2];

    // compute gradients 
    // dg = 2*evec dot rotation gradient
    grad[0] += 2*(eVec[2]*rd[1] - eVec[1]*rd[2]);
    grad[1] += 2*(eVec[0]*rd[2] - eVec[2]*rd[0]);
    grad[2] += 2*(eVec[1]*rd[0] - eVec[0]*rd[1]);

    // compute jacobian
    jac[0][0] += 2 * (rd[1]*eVec[1]
                   +  rd[2]*eVec[2]
                   +  rd[2]*rd[2]
 		   +  rd[1]*rd[1]);
  
    jac[0][1] += -2 * (rd[0]*eVec[1]
                    +  rd[1]*rd[0]);

    jac[0][2] += -2 * (rd[0]*eVec[2]
                    +  rd[0]*rd[2]);

    jac[1][1] +=  2 * (rd[0]*eVec[0]
                    +  rd[2]*eVec[2]
                    +  rd[2]*rd[2]
		    +  rd[0]*rd[0]);

    jac[1][2] += -2 * (rd[1]*eVec[2]
                    +  rd[2]*rd[1]);

    jac[2][2] +=  2 * (rd[0]*eVec[0]
                    +  rd[1]*eVec[1]
		    +  rd[1]*rd[1]
  		    +  rd[0]*rd[0]);
  }

  // Fill in symmetric terms
  jac[1][0] = jac[0][1];
  jac[2][0] = jac[0][2];
  jac[2][1] = jac[1][2];
}

void GeomState::computeCG(double cg[3])
{
  // init cg
  cg[0] = 0.0;
  cg[1] = 0.0;
  cg[2] = 0.0;

  int i;
  for ( i = 0; i < numnodes; ++i)  {
    if (flag[i])  {
      cg[0] += ns[i].x;
      cg[1] += ns[i].y;
      cg[2] += ns[i].z;
    }
  }

  double invTotNd = 1.0 / numReal;
  for (i = 0; i < 3; ++i)
    cg[i] *= invTotNd;
}

void GeomState::rotate(double R[3][3], double v[3]) 
{
  double c[3];
  for (int j = 0; j < 3; j++)
    c[j] = R[j][0]*v[0] +
           R[j][1]*v[1] +
           R[j][2]*v[2];
 
  v[0] = c[0];
  v[1] = c[1];
  v[2] = c[2];
}

void GeomState::getGlobalRot(double R[3][3]) 
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      R[i][j] = gRot[i][j];

}

double
NodeState::diff(const Node &un, int dof)
{
  double delta = 0.0;

  // Loop over all of nodes
  switch(dof) {
    case 0:
    {
      delta = x-un.x;
      break;
    }
    case 1:
    {
      delta = y-un.y;
      break;
    }
    case 2:
    {
      delta = z-un.z;
      break;
    }
    case 3:
    case 4:
    case 5:
    {
      double vec[3];
      mat_to_vec(R, vec);
      delta = vec[dof-3];
      break;
    }
  }
  return delta;
}

#include <Corotational.d/TemperatureState.h>
TemperatureState::TemperatureState(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs)
 : GeomState(cs)
{
  numnodes = dsa.numNodes();            // Number of nodes
  ns       = new NodeState[numnodes];   // Array of Node States
  loc      = new int[numnodes][1];      // dof locations        
  flag     = new bool[numnodes];     // flag for node to element connectivity
  numReal = 0;

  int i;
  for(i=0; i<numnodes; ++i) {

    // Store location of each degree of freedom
    loc[i][0] = cdsa.locate( i, DofSet::Temp );

    // Get Node i from the Coordinate (Node) set
    Node *node_i = cs[i];

    if (node_i)  {

      ns[i].x = 0;
      if(dsa[i].list() != 0)  {
        flag[i] = true;
        numReal++;
      }
      else
        flag[i] = false;

    }
    else  {
      // HB: TEMPORARY BAD FIX FOR DEALING WITH LAGRANGE MULTIPLIERS (RIGID BAR) !!!
      ns[i].x = 0.0;
      loc[i][0] = cdsa.locate( i, DofSet::Lagrange );
      flag[i] = false;
      //flag[i] = true;
    }

  }

  // Initialize Global Rotation Matrix to Identity
  //double zeroRot[3] = {0.0, 0.0, 0.0};
  //computeRotMat(zeroRot, gRot);
  //computeCG(refCG);

}

TemperatureState::TemperatureState(const TemperatureState &g2) : GeomState((CoordSet &) g2.X0)
{
  // Copy number of nodes
  numnodes = g2.numnodes;

  // Allocate memory for node states & dof locations
  ns  = new NodeState[numnodes];
  loc = new int[numnodes][1];

  // flag for node to element connectivity
  flag     = new bool[numnodes];

  numReal = g2.numReal;

  // Copy dof locations
  int i;
  for(i = 0; i < numnodes; ++i) {
    loc[i][0] = g2.loc[i][0];
  }

  // Copy node states
  for(i = 0; i < numnodes; ++i) {
    ns[i]  = g2.ns[i];
    flag[i]= g2.flag[i];
  }
}


void
TemperatureState::update(const Vector &v)
{
 // v = incremental displacement vector

 int i;
 for(i=0; i<numnodes; ++i) {

     // Set incremental translational displacements
     double dx = (loc[i][0] >= 0) ? v[loc[i][0]] : 0.0;

     // Increment total translational displacements
     ns[i].x += dx;

   }
}

void
TemperatureState::updatePrescribedDisplacement(BCond* dbc, int numDirichlet,
                                        double delta)
{
  // allocate space to store rotational prescribed dofs
  double (*dth)[3] = new double[numnodes][3];

  // initialize to zero, rotational prescribed dofs
  int i;
  for(i=0; i<numnodes; ++i)
    dth[i][0] = dth[i][1] = dth[i][2] = 0.0;

  for(i=0; i<numDirichlet; ++i) {

    int nodeNumber = dbc[i].nnum;
    int dofNumber  = dbc[i].dofnum;


    // we multiply the total prescribed value by delta which
    // is a parameter prescribed by the user in the input file
    // it is a load control parameter. By default, it is set to 1.0
    // which is applying the full prescribed displacement all at once.

    double prescribedValue = delta*dbc[i].val;

    // if prescribed value is zero, we do nothing.
    if( prescribedValue == 0.0 ) continue;

    switch(dofNumber) {
        case 6:
                ns[nodeNumber].x += prescribedValue;
                break;
        default:
                break;
    }

  }
  delete [] dth;
}

void
TemperatureState::get_inc_displacement(Vector &incVec, GeomState &ss, bool zeroRot)
{
  int inode;
  for(inode=0; inode<numnodes; ++inode) {
    // Update incremental temperature
    if(loc[inode][0] >= 0) incVec[loc[inode][0]] = ns[inode].x - ss.ns[inode].x;
  }
}

void
TemperatureState::midpoint_step_update(Vector &vel_n, Vector &accel_n, double delta, GeomState &ss,
                                       double beta, double gamma, double alphaf, double alpham)
{
 // XXXX THIS HASN'T BEEN UPDATED FOR GENERALIZED ALPHA YET
 // Update incremental temperatures at end of step:
 double coef = 2.0/delta;

 // x,y,z velocity step update (velocity_n+1 = 2.0*velocity_n+1/2 - velocity_n)
 // note: we are computing translational velocity_n+1/2 locally
 //       as velocity_n+1/2 = (ns[i] - ss.ns[i])/delta = 2/dt*(u_n+1/2 - u_n)

 int i;
 for(i=0; i<numnodes; ++i) {
   if(loc[i][0] >= 0)
     vel_n[loc[i][0]] = coef*(ns[i].x - ss.ns[i].x) - vel_n[loc[i][0]];
 }

 // Update step translational displacements
 int inode;
 for(inode=0; inode<numnodes; ++inode) {
   ns[inode].x    = 2.0*ns[inode].x - ss.ns[inode].x;
   ss.ns[inode].x = ns[inode].x;
 }
}

