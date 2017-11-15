#include "../Driver.d/Domain.h"
#include "../Driver.d/GeoSource.h"

// Dummy declaration for global sfem, etc
Sfem *sfem;
Communicator *heatStructCom;
Communicator *fluidCom;

// Dummy declaration for a global domain, to get an executable.
Domain *domain;
std::unique_ptr<GenSubDomainFactory<double> >   subDomainFactory(new GenSubDomainFactory<double>());
std::unique_ptr<GenSubDomainFactory<DComplex> > subDomainFactoryC(new GenSubDomainFactory<DComplex>());;

// Dummy declaration for a global geoSource
GeoSource *geoSource;

long totMemSky       = 0;
long totMemSparse    = 0;

bool nosa=false;   // no simulated anealing

#include <Interface.d/CU_Feti.h>
#include <complex>

void instantiate()
{
  double *X = 0;
  double *Y = 0;
  double *Z = 0;
  int nnodes_this_subdomain = 0;
  int *elem_ptr = 0;
  int *elem_conn = 0;
  int num_elem_this_subdomain = 0;
  int _neq = 0;
  int *_map = 0;
  int ndofs_deleted = 0;
  int *bc_node_ids = 0;
  int *bc_node_dof = 0;
  MPI_Comm *salinas_communicator = 0;
  int numads = 0;
  int *subdad = 0;
  int *com_ptr = 0;
  int *com_lis = 0;
  unsigned short *dofmap_on_node = 0;
  int *local_to_global_node_map = 0;
  //int *subs_ptr_proc = 0;
  //int *subs_per_proc = 0;
  //int local_subdomain_index = 0;
  FetiParams *params = 0;
  //int body = 0;

  CU_Feti<double> *dfeti = new CU_Feti<double>(X, Y, Z, nnodes_this_subdomain,
                                          elem_ptr, elem_conn, num_elem_this_subdomain,
                                          _neq, _map, ndofs_deleted, bc_node_ids,
                                          bc_node_dof, salinas_communicator, numads,
                                          subdad, com_ptr, com_lis, dofmap_on_node,
                                          local_to_global_node_map, params);

  int *A_aj = 0;
  int *A_ai = 0;
  //int localSubNumber = 0;
  // dfeti->setMatrixPattern(A_aj, A_ai, localSubNumber);

  double *Force = 0;
  double *U_real = 0;
  FetiValues returnvals;
  dfeti->Solve(Force, U_real, returnvals);

  dfeti->zeroLHS();

  double multiplier = 0;
  double *A_aa = 0;
  dfeti->addToLHS(multiplier, A_aa, A_aj, A_ai);

  int NumberMpc = 0;
  MpcLocal *MpcVector = 0;
  dfeti->setMPCs(NumberMpc, MpcVector);

  int nMpc = dfeti->NumberMpcConstraints();

  double *lambdas = 0;
  dfeti->MpcForces(lambdas);

  int num_vect = 0;
  const double *R[1];
  double **V = 0;
  int trans = 0;
  dfeti->ConstraintProduct(num_vect, R, V, trans);

  int numRBM = dfeti->getNumRBMs();
  double *rbms = 0;
  dfeti->getRBMs(rbms);

  int numContactNeigh = 0;
  int *contactNeighb = 0;
  int *neighbPtr = 0;
  enum face_type *ftype = 0;
  int *face_ptr = 0;
  int *face_conn = 0;
  double (*face_normal)[3] = 0;
  double (*gap_vector)[3] = 0;
  int *glbl_id = 0;
  const double *disp_updated = 0;
  double *contact_forces = 0;
  dfeti->setContactSurfaces(numContactNeigh, contactNeighb, neighbPtr, ftype,
                            face_ptr, face_conn, face_normal, gap_vector, glbl_id, 0, 0);
  dfeti->setContactNormal(face_normal);
  //dfeti->setContactGap(gap_vector);
  dfeti->setContactGap(disp_updated);
  dfeti->updateContactDisp(disp_updated);
  dfeti->getContactForce(contact_forces);
  int numTiedNeigh = 0;
  int *tiedNeighb = 0;
  enum face_type *faceEl_type = 0;
  int *faceEl_ptr = 0;
  int *faceEl_conn = 0;
  double (*faceEl_normal)[3] = 0;
  dfeti->setTiedSurfaces(numTiedNeigh, tiedNeighb, neighbPtr, faceEl_type, 
                         faceEl_ptr, faceEl_conn, faceEl_normal, gap_vector, glbl_id, 0, 0);
  nMpc = dfeti->updateLocalMpcs(&MpcVector);
  delete dfeti;


  CU_Feti<complex<double> > *cfeti = new CU_Feti<complex<double> >(X, Y, Z, nnodes_this_subdomain,
                                elem_ptr, elem_conn, num_elem_this_subdomain,
                                _neq, _map, ndofs_deleted, bc_node_ids,
                                bc_node_dof, salinas_communicator, numads,
                                subdad, com_ptr, com_lis, dofmap_on_node,
                                local_to_global_node_map, params);

  // cfeti.setMatrixPattern(A_aj, A_ai, localSubNumber);

  complex<double> *c_Force = 0;
  complex<double> *c_U_real = 0;
  cfeti->Solve(c_Force, c_U_real, returnvals);

  cfeti->zeroLHS();

  complex<double> c_multiplier = 0;
  cfeti->addToLHS(multiplier, A_aa, A_aj, A_ai);

  cfeti->setMPCs(NumberMpc, MpcVector);

  nMpc = cfeti->NumberMpcConstraints();

  cfeti->MpcForces(lambdas);

  complex<double> **c_V = 0;
  cfeti->ConstraintProduct(num_vect, R, c_V, trans);

  numRBM = cfeti->getNumRBMs();
  complex<double> *c_rbms = 0;
  cfeti->getRBMs(c_rbms);

  enum face_type *c_ftype = 0;
  cfeti->setContactSurfaces(numContactNeigh, contactNeighb, neighbPtr, c_ftype,
                            face_ptr, face_conn, face_normal, gap_vector, glbl_id, 0, 0);
  cfeti->setContactNormal(face_normal);
  //cfeti->setContactGap(gap_vector);
  cfeti->setContactGap(disp_updated);
  cfeti->updateContactDisp(disp_updated);
  cfeti->getContactForce(contact_forces);

  enum face_type *c_faceEl_type = 0;
  cfeti->setTiedSurfaces(numTiedNeigh, tiedNeighb, neighbPtr, c_faceEl_type,
                         faceEl_ptr, faceEl_conn, faceEl_normal, gap_vector, glbl_id, 0, 0);
  nMpc = cfeti->updateLocalMpcs(&MpcVector);
  delete cfeti;

}


int Main(int argc, char** argv) 
{
  
  int myrank = -1;
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  cerr << "myrank = " << myrank << endl;
#endif

  bool xxx = false;
  if(xxx) instantiate();

  // hardcoded test
  double X[8] = { 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0 };
  double Y[8] = { 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0 };
  double Z[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  int nnodes_this_subdomain = 8;
  int elem_ptr[4] = { 0, 4, 8, 12 };
  int elem_conn[12] = { 0, 1, 5, 4,  1, 2, 6, 5,  2, 3, 7, 6 };
  int num_elem_this_subdomain = 3;
  int _neq = 0;
  int *_map = 0;
  int ndofs_deleted = 12;
  int bc_node_ids[12] = { 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4 };
  int bc_node_dof[12] = { 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5 };
  MPI_Comm *salinas_communicator = 0;
  int numads = 1;
  int subdad[1];
  int com_ptr[2] = { 0, 8 };
  int com_lis[8];
  unsigned short dofmap_on_node[8] = { 63, 63, 63, 63, 63, 63, 63, 63 };
  int local_to_global_node_map[8];
  //int body;

 if(myrank == 0) {
  subdad[0] = 1;
  for(int i=0; i<8; ++i) {  
    com_lis[i] = i; 
    local_to_global_node_map[i] = i;
  }
  //body = 0;
 }
 else { // myrank = 1
  subdad[0] = 0;
  for(int i=0; i<8; ++i) {  
    com_lis[i] = i+8;
    local_to_global_node_map[i] = i+8;
  }
  //body = 1;
 } 

  //int subs_ptr_proc[3] = { 0, 1, 2 };
  //int subs_per_proc[2] = { 0, 1 };
  //int local_subdomain_index = 0;
  FetiParams *params = new FetiParams();

  CU_Feti<double> *dfeti = new CU_Feti<double>(X, Y, Z, nnodes_this_subdomain,
                                elem_ptr, elem_conn, num_elem_this_subdomain,
                                _neq, _map, ndofs_deleted, bc_node_ids,
                                bc_node_dof, salinas_communicator, numads,
                                subdad, com_ptr, com_lis, dofmap_on_node,
                                local_to_global_node_map, params);

  int numTiedNeigh = 1;
  int neighbPtr[2] =  { 0, 3 };
  enum face_type faceEl_type[3] = { QUAD4_FACE, QUAD4_FACE, QUAD4_FACE };
  int faceEl_ptr[4] = { 0, 4, 8, 12 };
  double gap_vector[3][3] = { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } };

 if(myrank == 0) {
  int tiedNeighb[1] = { 1 };
  int faceEl_conn[12] = { 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6 };
  double faceEl_normal[3][3] = { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } };
  int glbl_id[3] = { 0, 0, 0 };
  dfeti->setTiedSurfaces(numTiedNeigh, tiedNeighb, neighbPtr, faceEl_type,
                         faceEl_ptr, faceEl_conn, faceEl_normal, gap_vector, glbl_id, 0, 0);
 }
 else { // myrank = 1
  int tiedNeighb[1] = { 0 };
  int faceEl_conn[12] = { 0, 4, 5, 1, 1, 5, 6, 2, 2, 6, 7, 3 };
  double faceEl_normal[3][3] = { { 0.0, 0.0, -1.0 }, { 0.0, 0.0, -1.0 }, { 0.0, 0.0, -1.0 } };
  int glbl_id[3] = { 1, 1, 1 };
  dfeti->setTiedSurfaces(numTiedNeigh, tiedNeighb, neighbPtr, faceEl_type,
                         faceEl_ptr, faceEl_conn, faceEl_normal, gap_vector, glbl_id, 0, 0);
 }

 int NumberMpc = 1;
 MpcLocal *MpcVector = new MpcLocal[1];
 if(myrank == 0) {
   // 1st MPC on CPU 0: x components of global nodes 1 and 2 are equal
   MpcVector[0].NumEntries(2);
   MpcVector[0].NumEntriesGlobal(2);
   MpcVector[0].LocalId(0,0); MpcVector[0].LocalId(1,1);
   MpcVector[0].Coef(0,1.0); MpcVector[0].Coef(1,-1.0);
   MpcVector[0].NodalDof(0,0); MpcVector[0].NodalDof(1,0);
   MpcVector[0].G(0.0);
 }
 else { // myrank = 1
   MpcVector[0].NumEntries(0);
   MpcVector[0].NumEntriesGlobal(2);
 }
 dfeti->setMPCs(NumberMpc, MpcVector);  // this can be called before or after setTiedSurfaces

 dfeti->zeroLHS();  

#ifdef USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
