#include <stdlib.h>
#include <Utils.d/dbg_alloca.h>
#include <map>

#include <Comm.d/Communicator.h>
#include <Feti.d/Feti.h>
#include <Utils.d/Memory.h>
#include <Interface.d/MpcLocal.h>
#include <Element.d/Element.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include <Utils.d/dofset.h>

extern SysCom *scom;
Communicator *structCom;
int salinasFlag = 1;

#ifdef TFLOP
extern "C" int heap_info(int*, int*, int*, int*);
#endif

template<class Scalar>
GenSandiaDomain<Scalar>::GenSandiaDomain(Domain *d, const double *_x, const double *_y, const double *_z,
                                         int _num_nodes, int *_elem_ptr, int *_elem_conn, int _num_elems,
                                         int _neq, const int *_map, int _ndofs_deleted, const int *_bc_node_ids,
                                         int *_bc_node_dof, int _numads, int *_procad, int *_com_ptr, int *_com_lis,
                                         unsigned short *_dofmap_on_node, int *_global_node_nums, FetiParams *params)
 : GenDecDomain<Scalar>(d)
{
  x = _x; y = _y; z = _z;
  num_nodes = _num_nodes;
  elem_ptr = _elem_ptr;
  elem_conn = _elem_conn;
  num_elems = _num_elems;
  neq = _neq;
  map = _map;
  ndofs_deleted = _ndofs_deleted;
  bc_node_ids = _bc_node_ids;
  bc_node_dof = _bc_node_dof;
  numads = _numads;
  procad = _procad;
  com_ptr = _com_ptr;
  com_lis = _com_lis;
  dofmap_on_node = _dofmap_on_node;
  global_node_nums = _global_node_nums;

  initialize();

  this->numCPU = scom->numCPUs();
  this->myCPU = scom->myID();
  
  // copy data from FetiParams
  problem_type = params->Solution_type(); // 0 = static, 1 = eigen, 2 = dynamic, 3 = directfrf
                                          // 4 = eigen for modal acceleration (need to compute rbms)
  computeRbms = (params->Nrbms() > 0) ? true : false;
  //computeRbms = (params->Nrbms() > 0 || problem_type == 0) ? true : false; // PJSA 9-20-07
  setFetiInfo(params);
  setSolInfo(params);
  mortar_type = params->Mortar_type();
  have_contact = false;
}


template<class Scalar>
void
GenSandiaDomain<Scalar>::initialize()
{
  finfo = 0; 
  numberMpc = 0; numberCtc = 0; 
  sysMatrix = 0; sparseK = 0;
  fetiSolver = 0;
  local_node_nums = 0;
  initial_disp = 0;
  //dx = dy = dz = 0;
}

template<class Scalar>
GenSandiaDomain<Scalar>::~GenSandiaDomain()
{
  if(this->cpuToSub) { delete this->cpuToSub; this->cpuToSub = 0; this->localSubToGl = 0; }
  if(sysMatrix) { /*for(int i=0; i<this->numSub; ++i) if(sysMatrix) delete sysMatrix;*/ delete [] sysMatrix; }
  if(sparseK) { /*for(int i=0; i<this->numSub; ++i) if(sparseK) delete sparseK;*/ delete [] sparseK; }
  if(local_node_nums) { delete [] local_node_nums; local_node_nums = 0; }
  if(initial_disp) { delete [] initial_disp; initial_disp = 0; }
  //if(dx) delete [] dx; if(dy) delete [] dy; if(dz) delete [] dz;
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::setFetiInfo(FetiParams *params)
{
  finfo = new FetiInfo();
  finfo->version = FetiInfo::fetidp;
  domain->solInfo().type = 2;
  if(params->Verbose_flag() != 0) verboseFlag = 1;
  finfo->printNumber = (verboseFlag) ? 1 : -1;
  if((params->Verbose_flag() & (1 << 4)) != 0) finfo->contactPrintFlag = 3;
  //finfo->c_normalize = false; // XXXX DEBUG CONTACT
  finfo->maxit = params->Max_iterations();
  finfo->tol = params->Solver_tol();
  if(params->Max_N_orthog_vecs() > 0)
    finfo->maxortho = params->Max_N_orthog_vecs();
  else
    finfo->maxortho = -params->Max_N_orthog_vecs();
  finfo->nonLocalQ = params->Projector();
  finfo->grbm_tol  = params->Grbm_tol();
  finfo->crbm_tol  = params->Crbm_tol();
  finfo->nlPrecFlg = params->NlPrecFlg();
  finfo->cct_tol = params->Mpc_method().cct_tol;
  finfo->rebuildcct = int(params->Contact_params().rebuild_cct);
  finfo->geometric_gap = params->Contact_params().geometric_gap;
  finfo->initial_lambda = params->Initial_lambda();
  finfo->stagnation_tol = params->Stagnation_tol();
  finfo->useMRHS = (params->Multiple_rhs()) ? 1 : 0;
  finfo->rescalef = false; // this seems to work better

  if(params->Corner_dimensionality() == 3)
    finfo->corners = FetiInfo::noEndCorners3;
  else finfo->corners = FetiInfo::noEndCorners6;

  if((scom->myID() == 0)) {
    fprintf(stdout,"-----------------------------------------------"
                   "------------------------------\n");
    if(params->Solution_type() == 0)
      fprintf(stdout," Solution Type                     = Statics\n");
    else if(params->Solution_type() == 1 || params->Solution_type() == 4)
      fprintf(stdout," Solution Type                     = Eigen\n");
    else if(params->Solution_type() == 2)
      fprintf(stdout," Solution Type                     = Implicit Transient\n");
    else if(params->Solution_type() == 3)
      fprintf(stdout," Solution Type                     = Direct Frequency Response\n");
    fprintf(stdout," Number of CPUs                    = %d\n",this->numCPU);
    fprintf(stdout," Maximum Number of Iterations      = %d\n",finfo->maxit);
    fprintf(stdout," Maximum Size of Reortho. Vectors  = %d\n",finfo->maxortho);
    fprintf(stdout," Tolerance for Convergence         = %e\n",finfo->tol);
    fprintf(stdout," Tolerance for Kcc                 = %e\n",finfo->crbm_tol);
    fprintf(stdout," Tolerance for GtG                 = %e\n",finfo->grbm_tol);
    fprintf(stdout," Mechanical Tolerance              = %e\n",params->Rbm_tol_mech());
    fprintf(stdout," SVD Tolerance                     = %e\n",params->Rbm_tol_svd());
    fprintf(stdout," Verbose Flag                      = %d\n",params->Verbose_flag());
    fprintf(stdout," Multiple LHS acceleration         = %d\n",finfo->nlPrecFlg);
    fprintf(stdout," Multiple RHS acceleration         = %d\n",finfo->useMRHS);
    fprintf(stdout," Corner Dimensionality             = %d\n",params->Corner_dimensionality());
  }

  // DPH overrides
  if(params->Feti_level() == 4 ||
     ((params->Feti_level() == 2) && (params->Shift() > 0.0))) { // PJSA 2-17-05 shifted eigen
    if((params->Feti_level() == 2) && (scom->myID() == 0))
      cerr << "*** Shifted Eigen analysis. Activating FETI-DPH. Shift = " << params->Shift() << " *** \n";
    finfo->dph_flag = true;
    if(scom->myID() == 0) cerr << "*** Setting Kcc tolerance to 1.0e-12 \n";
    finfo->crbm_tol  = 1.0e-12;
  }
  else {
    finfo->dph_flag = false;
  }

  switch(params->Corner_augmentation()) {
    default:
    case no_augmentation:
      finfo->augment = FetiInfo::none;
      finfo->rbmType = FetiInfo::None;
      finfo->nGs = 0;
      if(scom->myID() == 0) fprintf(stdout," Corner Augmentation               = None\n");
      break;
    case subdomain_augmentation:
      cerr << " *** WARNING: subdomain_augmentation is not supported \n";
      break;
    case edge_augmentation:
      finfo->augment = FetiInfo::Edges;
      if(scom->myID() == 0) fprintf(stdout," Corner Augmentation               = Edge\n");
      break;
    case weighted_edge_augmentation:
      finfo->augment = FetiInfo::WeightedEdges;
      if(scom->myID() == 0) fprintf(stdout," Corner Augmentation               = Weighted Edge\n");
      break;
  }

  if(params->Corner_augmentation() != no_augmentation) {
    switch(params->Corner_aug_rbm_type()) {
      default:
      case translation:
        finfo->rbmType = FetiInfo::translation;
        finfo->nGs = 3;
        if(scom->myID() == 0) fprintf(stdout," Corner Augmentation RBM Type      = Translation\n");
        break;
      case all:
        finfo->rbmType = FetiInfo::all;
        finfo->nGs = 6;
        if(scom->myID() == 0) fprintf(stdout," Corner Augmentation RBM Type      = All\n");
        break;
      case none:
        finfo->rbmType = FetiInfo::None;
        finfo->nGs = 0;
        if(scom->myID() == 0) fprintf(stdout," Corner Augmentation RBM Type      = None\n");
        break;
    }
  }

  if(finfo->dph_flag) {
    int numdir = params->NumWaveDirections();
    if(numdir > 0) {
      if(finfo->augment == FetiInfo::none) finfo->augment = FetiInfo::Edges;
      finfo->numdir = numdir;
      finfo->orthotol = params->Wave_ortho_tol();
      finfo->waveMethod = FetiInfo::averageK;
      finfo->spaceDimension = params->Space_dimension();
      if(scom->myID() == 0) {
         fprintf(stdout," Wave Augmentation                 = On\n");
         fprintf(stdout," Number of Wave Directions         = %d\n", numdir);
         fprintf(stdout," Space Dimension                   = %d\n", params->Space_dimension());
         fprintf(stdout," Wave Ortho Tolerance              = %f\n",params->Wave_ortho_tol());
      }
    }
    else {
      if(scom->myID() == 0) fprintf(stdout," Wave Augmentation                 = Off\n");
    }
  }

  switch(params->Outerloop_solver()) {
    default:
    case cg:
      finfo->outerloop = FetiInfo::CG;
      if(scom->myID() == 0) fprintf(stdout," Outerloop Solver                  = CG\n");
      break;
    case gmres:
      finfo->outerloop = FetiInfo::GMRES;
      if(scom->myID() == 0) fprintf(stdout," Outerloop Solver                  = GMRES\n");
      break;
    case gcr:
      finfo->outerloop = FetiInfo::GCR;
      if(scom->myID() == 0) fprintf(stdout," Outerloop Solver                  = GCR\n");
      break;
    case cgal:
      finfo->outerloop = FetiInfo::CGAL;
      finfo->expansion = 2; // make this the default for CGAL
      finfo->alphabar_cntl = params->Contact_params().alphabar;
      finfo->maxinnerit = params->Max_iterations();
      finfo->maxouterit = (params->Max_iterations()/10 < 10) ? 10 : params->Max_iterations()/10;
      if(scom->myID() == 0) fprintf(stdout," Outerloop Solver                  = CGAL\n");
      break;
  }

  switch(params->Local_solver()) {
    default:
    case Skyline:
      finfo->solvertype = FetiInfo::skyline;
      if(scom->myID() == 0) fprintf(stdout," Local Solver                      = Skyline (Sloan) \n");
      break;
    case Sparse:
      finfo->solvertype = FetiInfo::sparse;
      if(scom->myID() == 0) fprintf(stdout," Local Solver                      = Sparse (Esmond) \n");
      break;
    case Spooles:
      finfo->solvertype = FetiInfo::spooles;
      if(scom->myID() == 0) fprintf(stdout," Local Solver                      = Spooles\n");
      break;

  }

  switch(params->Coarse_solver()) {
    default:
    case Skyline:
      finfo->gtgSolver = FetiInfo::skyline;
      if(scom->myID() == 0) fprintf(stdout," Coarse Solver                     = Skyline (Sloan) \n");
      break;
    case Sparse:
      finfo->gtgSolver = FetiInfo::sparse;
      if(scom->myID() == 0) fprintf(stdout," Coarse Solver                     = Sparse (Esmond) \n");
      break;
    case Spooles:
      finfo->gtgSolver = FetiInfo::spooles;
      if(scom->myID() == 0) fprintf(stdout," Coarse Solver                     = Spooles\n");
      break;
  }

  if((scom->myID() == 0) && ((params->Local_solver() == Spooles) || (params->Coarse_solver() == Spooles))) {
    if(params->Pivoting() == 0)
      fprintf(stdout," Spooles Pivoting                  = Off \n");
    else
      fprintf(stdout," Spooles Pivoting                  = On \n");
  }

  switch(params->Precondition()) {
    default:
    case PRECOND_DIRICHLET:
      finfo->precno = FetiInfo::dirichlet;
      if(scom->myID() == 0) fprintf(stdout," Preconditioner                    = Dirichlet\n");
      break;
    case PRECOND_LUMPED:
      finfo->precno = FetiInfo::lumped;
      if(scom->myID() == 0) fprintf(stdout," Preconditioner                    = Lumped\n");
      break;
    case 0:
      finfo->precno = FetiInfo::noPrec;
      if(scom->myID() == 0) fprintf(stdout," Preconditioner                    = None\n");
      break;
  }

  switch(params->Weighting()) {
    default:
    case Stiffness:
      finfo->scaling = FetiInfo::kscaling;
      if(scom->myID() == 0) fprintf(stdout," Preconditioner Weighting          = Stiffness\n");
      break;
    case Topological:
      finfo->scaling = FetiInfo::tscaling;
      if(scom->myID() == 0) fprintf(stdout," Preconditioner Weighting          = Topological\n");
      break;
  }

  int rbmFlag = 0;
  switch(params->Rbm_method()) {
    default:
    case RBM_ALGEBRAIC:
      rbmFlag = 0;
      if(scom->myID() == 0) fprintf(stdout," RBM Method                        = Algebraic\n");
      break;
    case RBM_GEOMETRIC:
      rbmFlag = 1;
      if(scom->myID() == 0) fprintf(stdout," RBM Method                        = Geometric\n");
      break;
  }
  domain->solInfo().rbmflg = rbmFlag;

  switch(params->Mpc_method().method) {
    default:
    case MPC_Method::Dual:
      finfo->mpcflag = 1;
      if(scom->myID() == 0) fprintf(stdout," MPC Method                        = Dual\n");
 
      switch(params->Mpc_method().submethod) {
        default:
        case MPC_Method::Auto:
          finfo->mpc_precno = FetiInfo::autoSelectCCt;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner                = Auto Select\n");
          break;
        case MPC_Method::Full:
          finfo->mpc_precno = FetiInfo::globalCCt;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner                = Full\n");
          break;
        case MPC_Method::None:
          finfo->mpc_precno = FetiInfo::noMpcPrec;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner                = None\n");
          break;
        case MPC_Method::Diag:
          finfo->mpc_precno = FetiInfo::diagCCt;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner                = Diagonal\n");
          break;
        case MPC_Method::BlockDiag:
          finfo->mpc_precno = FetiInfo::blockDiagCCt;
          finfo->mpc_block = FetiInfo::topoBlock;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner                = Block Diagonal\n");
          break;
        case MPC_Method::PerSub:
          finfo->mpc_precno = FetiInfo::blockDiagCCt;
          finfo->mpc_block = FetiInfo::subBlock;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner                = SubBlock Diagonal\n");
          break;
        case MPC_Method::PerFace:
          finfo->mpc_precno = FetiInfo::blockDiagCCt;
          finfo->mpc_block = FetiInfo::mortarBlock;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner                = EdgeBlock Diagonal\n");
          break;
      }
 
      switch(params->Mpc_method().solver) {
        default:
        case Skyline:
          finfo->cctSolver = FetiInfo::skyline;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner Solver         = Skyline\n");
          break;
        case Sparse:
          finfo->cctSolver = FetiInfo::sparse;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner Solver         = Sparse\n");
          break;
        case Spooles:
          finfo->cctSolver = FetiInfo::spooles;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner Solver         = Spooles\n");
          break;
      }

      switch(params->Mpc_method().weighting) {
        default:
        case Topological:
          finfo->mpc_scaling = FetiInfo::tscaling;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner Weighting      = Topological\n");
          break;
        case Stiffness:
          finfo->mpc_scaling = FetiInfo::kscaling;
          if(scom->myID() == 0) fprintf(stdout," MPC Preconditioner Weighting      = Stiffness\n");
          break;
      }

      if(scom->myID() == 0)     fprintf(stdout," Tolerance for MPC Preconditioner  = %e\n",finfo->cct_tol);

      break;
    case MPC_Method::Primal:
      finfo->mpcflag = 2;
      if(scom->myID() == 0) fprintf(stdout," MPC Method:                       = Primal\n");
      break;
  }

  // contact parameters
  if(finfo->outerloop == 0 || finfo->outerloop == 3) { // projected conjugate gradient or augmented lagrangian conjugate gradient
    finfo->gamma = params->Contact_params().gamma;
    if(scom->myID() == 0)     fprintf(stdout," Proportioning Tolerance           = %e\n",finfo->gamma);
    finfo->equi_tol = params->Contact_params().eq_tol;
    if(scom->myID() == 0)     fprintf(stdout," Equality Constraint Tolerance     = %e\n",finfo->equi_tol);
    finfo->iequ_tol = params->Contact_params().ieq_tol;
    if(scom->myID() == 0)     fprintf(stdout," Inequality Constraint Tolerance   = %e\n",finfo->iequ_tol);
  }
  if(finfo->outerloop == 3) { // augmented lagrangian conjugate gradient
    finfo->rho_cntl = params->Contact_params().rho;
    if(scom->myID() == 0)     fprintf(stdout," Initial Penalty Parameter         = %e\n",finfo->rho_cntl);
    finfo->beta = params->Contact_params().beta;
    if(scom->myID() == 0)     fprintf(stdout," Penalty Update Factor             = %e\n",finfo->beta);
    finfo->M = params->Contact_params().M;
    if(scom->myID() == 0)     fprintf(stdout," Adaptive Precision Control (M)    = %e\n",finfo->M);
    finfo->eta = params->Contact_params().eta;
    if(scom->myID() == 0)     fprintf(stdout," Adaptive Precision Control (eta)  = %e\n",finfo->eta);
    if(params->Contact_params().cgal_prec > 0) {
      finfo->cgal_prec = true;
      if(scom->myID() == 0)     fprintf(stdout," CGAL Preconditioner               = On\n");
    }
    else {
      finfo->cgal_prec = false;
      if(scom->myID() == 0)     fprintf(stdout," CGAL Preconditioner               = Off\n");
    }
  }

  if((scom->myID() == 0)) {
    fprintf(stdout,"-----------------------------------------------"
                   "------------------------------\n");
  }
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::setSolInfo(FetiParams *params)
{ 
  domain->solInfo().renum = 1;    // Sloan renumbering for skyline solvers
  domain->solInfo().sparse_renum = 0;  // Esmond renumbering for sparse solvers
  domain->solInfo().pivot = params->Pivoting(); // pivoting for spooles solvers
  domain->solInfo().getFetiInfo() = *finfo;
  delete finfo; finfo = &domain->solInfo().getFetiInfo(); // make sure we have pointer to the same object as domain->solInfo
  domain->solInfo().trbm   = params->Rbm_tol_mech();
  domain->solInfo().tolsvd = params->Rbm_tol_svd();
} 

template<class Scalar>
void
GenSandiaDomain<Scalar>::makeSubD()
{
 // determine number of components on this cpu
  Connectivity *elemToNode = new Connectivity(num_elems, elem_ptr, elem_conn, 0); // 0: don't delete
  Connectivity *nodeToElem = elemToNode->reverse();
  Connectivity *nodeToNode = nodeToElem->transcon(elemToNode);
  compStruct renumber = nodeToNode->renumByComponent(1);
  this->numSub = renumber.numComp;
  this->globalNumSub = this->communicator->globalSum(this->numSub);
  oneSubPerCPU = bool(this->communicator->globalMax(this->numSub) == 1); // special case

  this->subDomain = new GenSubDomain<Scalar>*[this->numSub];
  sysMatrix = new GenSolver<Scalar>*[this->numSub];
  sparseK = new GenSparseMatrix<Scalar>*[this->numSub];
  for(int i=0; i<this->numSub; ++i) { sysMatrix[i] = 0; sparseK[i] = 0; }

  // make cpuToSub
  makeCpuToSub();

  // make glSubToLocal and localSubToGl
  this->makeSubDMaps();

  if(this->numSub > 1) {
    // split nodes into subdomains
    renumber.order = new int[num_nodes];
    for(int i = 0; i < num_nodes; ++i) renumber.order[renumber.renum[i]] = i;
    Connectivity *subToNode_local = new Connectivity(this->numSub, renumber.xcomp, renumber.order);
    // split elements into subdomains
    Connectivity *subToElem = subToNode_local->transcon(nodeToElem);
    // split dofs into subdomain
    int *ptr = new int[num_nodes+1];
    /*int index = 0;*/ ptr[0] = 0;
    for(int i = 0; i < num_nodes; ++i) {
      int flags = int(dofmap_on_node[i]);
      DofSet dofs(flags);
      // unmark constrained dofs
      for(int j=0; j<ndofs_deleted; ++j) {
        if(bc_node_ids[j] == i) { // PJSA 12-6-05
          int flag = (1 << bc_node_dof[j]);
          dofs.unmark(flag);
        }
      }
      ptr[i+1] = ptr[i]+dofs.count();
    }
    int *tgt = new int[ptr[num_nodes]];
    for(int i = 0; i < ptr[num_nodes]; ++i) tgt[i] = i;

    Connectivity *nodeToDof = new Connectivity(num_nodes, ptr, tgt);
    Connectivity *subToDof = subToNode_local->transcon(nodeToDof);

    for(int i = 0; i < this->numSub; ++i) {  // this could be done in parallel
      // extract subdomain node data
      int sub_num_nodes = subToNode_local->num(i);
      double *sub_x = new double[sub_num_nodes];
      double *sub_y = new double[sub_num_nodes]; 
      double *sub_z = new double[sub_num_nodes];
      int *sub_global_node_nums = new int[sub_num_nodes];
      unsigned short int *sub_dofmap_on_node = new unsigned short int[sub_num_nodes];  // CURRENTLY NOT DELETED ANYWHERE !!
      int *local_sub_node_nums = new int[num_nodes];
      for(int j = 0; j < sub_num_nodes; ++j) {
        int k = (*subToNode_local)[i][j];
        sub_x[j] = x[k];
        sub_y[j] = y[k];
        sub_z[j] = z[k];
        sub_global_node_nums[j] = global_node_nums[k];
        sub_dofmap_on_node[j] = int(dofmap_on_node[k]);
        local_sub_node_nums[k] = j;
      }
      // extract subdomain element data
      int sub_num_elems = subToElem->num(i);
      int *sub_elem_ptr = new int[sub_num_elems+1];
      int numtarget = 0; sub_elem_ptr[0] = 0;
      for(int j = 0; j < sub_num_elems; ++j) {
        int k = (*subToElem)[i][j];
        numtarget += elem_ptr[k+1] - elem_ptr[k];
        sub_elem_ptr[j+1] = numtarget;
      }
      int *sub_elem_conn = new int[numtarget];
      numtarget = 0;
      for(int j = 0; j < sub_num_elems; ++j) {
        int k = (*subToElem)[i][j];
        for(int l=elem_ptr[k]; l<elem_ptr[k+1]; ++l) sub_elem_conn[numtarget++] = local_sub_node_nums[elem_conn[l]];
      }
      int *sub_map = new int[subToDof->num(i)];
      for(int j=0; j<subToDof->num(i); ++j) sub_map[j] = (*subToDof)[i][j];    

      // construct the subdomain
      int subNumber = (*this->cpuToSub)[this->myCPU][i];  // global subdomain id
      GenSandiaSubD<Scalar> *subd = new GenSandiaSubD<Scalar>(subNumber, i, sub_num_nodes, sub_x, sub_y, sub_z, 
                                                              sub_global_node_nums, sub_dofmap_on_node,
                                                              sub_num_elems, sub_elem_ptr, sub_elem_conn,
                                                              map, sub_map);
      this->subDomain[i] = subd;
     
      // delete
      delete [] sub_x; delete [] sub_y; delete [] sub_z;
      delete [] sub_elem_ptr; delete [] sub_elem_conn;
      delete [] local_sub_node_nums;
    }
    delete subToNode_local; 
    delete subToElem; delete nodeToDof; delete subToDof; 
  }
  else {
    delete [] renumber.xcomp;
    int subNumber = (*this->cpuToSub)[this->myCPU][0];  // global subdomain id
    GenSandiaSubD<Scalar> *subd = new GenSandiaSubD<Scalar>(subNumber, 0, num_nodes, x, y, z,
                                                            global_node_nums, dofmap_on_node, 
                                                            num_elems, elem_ptr, elem_conn, map, 0);
    this->subDomain[0] = subd;
  }
  paralApply(this->numSub, this->subDomain, &BaseSub::makeDSA);
  
  delete [] renumber.renum;
  delete elemToNode; delete nodeToElem; delete nodeToNode;
  findNumGlobNodes(); 

  makeSubToSubEtc(); 
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::makeCpuToSub()
{
  int *ptr = new int[this->numCPU+1];
  int *tgt = new int[this->globalNumSub];
  if(oneSubPerCPU) {  // 1 subdomain per CPU
    for(int i = 0; i < this->numCPU; ++i)
      ptr[i] = tgt[i] = i;
    ptr[this->numCPU] = this->numCPU;
  }
  else {  // multiple subdomains per CPU
    int *numSubsOnCpu = new int[this->numCPU];
    for(int i=0; i<this->numCPU; ++i) numSubsOnCpu[i] = 0;
    numSubsOnCpu[this->myCPU] = this->numSub;
    this->communicator->globalSum(this->numCPU, numSubsOnCpu);
    ptr[0] = 0;
    for(int i=0; i<this->numCPU; ++i)
      ptr[i+1] =  ptr[i] + numSubsOnCpu[i];
    for(int i=0; i<this->globalNumSub; ++i) tgt[i] = i;
    delete [] numSubsOnCpu;
  }
  this->cpuToSub = new Connectivity(this->numCPU, ptr, tgt);
  geoSource->setCpuToSub(this->cpuToSub);  
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::makeSComm()
{
  if(oneSubPerCPU) { 
    Connectivity *interfNodes = new Connectivity(numads, com_ptr, com_lis,0);
    GenSubDomain<Scalar> **subds = new GenSubDomain<Scalar> * [numads];
    int *remoteId = new int[numads];
    int *subNums = new int[numads];
    for(int i=0; i<numads; ++i) {
      // all neighbors are on remote cpus
      subds[i] = NULL;
      remoteId[i] = -1;
      subNums[i] = procad[i];
    }
    SComm *sc = new SComm(numads, subNums, subds, remoteId, interfNodes->copy());
    delete [] subds;
    sc->locSubNum = 0;
    sc->glSubToLocal = this->glSubToLocal;
    this->subDomain[0]->setSComm(sc);
    this->subDomain[0]->renumberSharedNodes();
  }
  else {
    this->getSharedNodes();
  }
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::makeSubToSubEtc()
{
  // make subtoNode, nodeToSub and subToSub connectivities
  if(oneSubPerCPU) {
    this->subToNode = this->nodeToSub = 0;
    int *ptr = new int[this->globalNumSub+1];
    for(int i =0; i < this->globalNumSub; ++i) ptr[i] = 0;
    ptr[this->myCPU] = numads+1;
    this->communicator->globalSum(this->globalNumSub, ptr);
    int total = 0;
    for(int i = 0; i < this->globalNumSub; ++i) {
      int tmp = ptr[i];
      ptr[i] = total;
      total += tmp;
    }
    ptr[this->globalNumSub] = total;
    int *tg = new int[total];
    for(int i = 0; i < total; ++i) tg[i] = 0;
    tg[ptr[this->myCPU]] = this->myCPU;
    for(int i = 0; i < numads; ++i) tg[i+1+ptr[this->myCPU]] = procad[i];
    this->communicator->globalSum(total, tg);
    this->subToSub = new Connectivity(this->globalNumSub, ptr, tg);
  }
  else {
    int *nn_perSub = new int[this->globalNumSub];
    for(int i=0; i<this->globalNumSub; ++i) nn_perSub[i] = 0;
    for(int i=0; i<this->numSub; ++i) nn_perSub[this->subDomain[i]->subNum()] = this->subDomain[i]->numNodes();
    this->communicator->globalSum(this->globalNumSub, nn_perSub);
    int *ptr =  new int[this->globalNumSub+1];
    ptr[0] = 0;
    for(int i=0; i<this->globalNumSub; ++i) ptr[i+1] = ptr[i] + nn_perSub[i];
    delete [] nn_perSub;
    int *tgt = new int[ptr[this->globalNumSub]];
    for(int i=0; i<ptr[this->globalNumSub]; ++i) tgt[i] = 0;
    for(int i=0; i<this->numSub; ++i) {
      for(int j=0; j<this->subDomain[i]->numNodes(); ++j) {
        tgt[ptr[this->subDomain[i]->subNum()]+j] = this->subDomain[i]->localToGlobal(j);
      }
    }
    this->communicator->globalSum(ptr[this->globalNumSub],tgt);
    this->subToNode = new Connectivity(this->globalNumSub, ptr, tgt);
    this->nodeToSub = this->subToNode->reverse();
    this->subToSub = this->subToNode->transcon(this->nodeToSub);
  }
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::setMPCs(int _numberMpc, MpcLocal *_mpclocal)
{
  numberMpc = _numberMpc; // this is the global number of MPCs
  mpclocal = _mpclocal;
  //if(this->myCPU == 0) cerr << "GenSandiaDomain<Scalar>::setMPCs, numberMpc = " << numberMpc << endl;

  // convert from Sandia MpcLocal to FEM LMPCons and add to domain
  for(int i=0; i<numberMpc; ++i) {
    LMPCons *lmpc = new LMPCons(i,mpclocal[i].G());
    for(int j = 0; j < mpclocal[i].NumEntries(); ++j) {
      LMPCTerm *term = new LMPCTerm(global_node_nums[mpclocal[i].LocalId(j)],
                                    mpclocal[i].NodalDof(j), mpclocal[i].Coef(j));
      lmpc->addterm(term);
    }
    //lmpc->type = 1; // XXXX DEBUG CONTACT make all mpcs inequalities
    domain->addLMPC(lmpc);
  }
}

template<class Scalar>
int GenSandiaDomain<Scalar>::updateLocalMpcs(MpcLocal **_mpclocal)
{
  // initialize coefs to zero
  for(int i = 0; i<numberMpc; ++i) {
    for(int k=0; k<mpclocal[i].NumEntries(); k++) {
      mpclocal[i].Coef(k,0.0);
    }
  }
  // loop sequentially over subdomains and add split coefs (2 subdomains may share same term)
  for(int i = 0; i<this->numSub; ++i) {
    for(int j = 0; j<this->subDomain[i]->numMPC; ++j) {
      int gj = this->subDomain[i]->localToGlobalMPC[j]; 
      mpclocal[gj].G(ScalarTypes::Real(this->subDomain[i]->mpc[j]->rhs));
      for(int k=0; k<mpclocal[gj].NumEntries(); k++) {
        double coef =  mpclocal[gj].Coef(k) + ScalarTypes::Real(this->subDomain[i]->mpc[j]->terms[k].coef);
        mpclocal[gj].Coef(k,coef);
      }
    }
    for(int j = 0; j<this->subDomain[i]->numMPC_primal; ++j) {
      int gj = this->subDomain[i]->localToGlobalMPC_primal[j];
      mpclocal[gj].G(ScalarTypes::Real(this->subDomain[i]->mpc_primal[j]->rhs));
      for(int k=0; k<mpclocal[gj].NumEntries(); k++) {
        double coef =  mpclocal[gj].Coef(k) + ScalarTypes::Real(this->subDomain[i]->mpc_primal[j]->terms[k].coef);
        mpclocal[gj].Coef(k,coef);
      }
    }
  }

  *_mpclocal = mpclocal;
  return numberMpc;
}  

template<class Scalar>
void
GenSandiaDomain<Scalar>::setWaveNumbers(FetiParams *params)
{
   paralApply(this->numSub, this->subDomain, &BaseSub::setWaveNumbers, params->WaveNumbers());
 
  // need to send & collect wave numbers for dph before calling SubDomain::makeQ
  if(finfo->waveMethod == FetiInfo::averageK) {
    FSCommPattern<double> *kPat = new FSCommPattern<double>(this->communicator, this->cpuToSub, this->myCPU, FSCommPattern<double>::CopyOnSend);
    for(int i=0; i<this->numSub; ++i) this->subDomain[i]->setCommSize(kPat, 3);
    kPat->finalize();
    paralApply(this->numSub, this->subDomain, &BaseSub::sendWaveNumbers, kPat);
    kPat->exchange();
    paralApply(this->numSub, this->subDomain, &BaseSub::collectWaveNumbers,kPat);
    delete kPat;
  } 
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::constrain()
{
  BCond *dbc = new BCond[ndofs_deleted];
  for(int i = 0; i < ndofs_deleted; ++i) {
    dbc[i].nnum = global_node_nums[bc_node_ids[i]];
    dbc[i].dofnum = bc_node_dof[i];
    dbc[i].val = 0.0;
  }
  //if(this->numSub == 1) { // don't have nodeToSub 
  if(oneSubPerCPU) { // XXXX
    this->subDomain[0]->setDirichlet(ndofs_deleted, dbc);
    this->subDomain[0]->renumberDirichlet();
    this->subDomain[0]->makeCDSA();
  }
  else {
    geoSource->setDirichlet(ndofs_deleted, dbc);
    this->preProcessBCsEtc();
  }
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::preProcess()
{
  this->makeCorners();
  this->getSharedDOFs();
  this->preProcessMPCs();
  this->getSharedMPCs();
  paralApply(this->numSub, this->subDomain, &BaseSub::mergeInterfaces);
  paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::splitData, this->cornerWeight);
  if(this->cornerWeight) { delete [] this->cornerWeight; this->cornerWeight = 0; }
  paralApply(this->numSub, this->subDomain, &GenSubDomain<Scalar>::initSrc);
  if(finfo->mpc_scaling == FetiInfo::kscaling)
    paralApply(this->numSub, this->subDomain, &GenSandiaSubD<Scalar>::setMpcGlobalTermIds, mpclocal, global_node_nums, oneSubPerCPU);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::initSysMatrix()
{
  if(fetiSolver) fetiSolver->resetOrthoSet();  // PJSA 9-18-2007
  execParal(this->numSub, this, &GenSandiaDomain<Scalar>::initSubSysMatrix);
}
 
template<class Scalar>
void
GenSandiaDomain<Scalar>::initSubSysMatrix(int iSub)
{
  SubDomain(iSub)->initSysMatrix(neq, sysMatrix[iSub], sparseK[iSub]);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::zeroSysMatrix()
{
  if(fetiSolver) fetiSolver->resetOrthoSet();  // PJSA 9-18-2007
  execParal(this->numSub, this, &GenSandiaDomain<Scalar>::zeroSubSysMatrix);     
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::zeroSubSysMatrix(int iSub)
{
  SubDomain(iSub)->zeroSysMatrix(sysMatrix[iSub], sparseK[iSub]);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::addSysMatrix(const double *Kaa, const int *Kaj,
                                      const int *Kai, Scalar multiplier, bool isK)
{
  execParal(this->numSub, this, &GenSandiaDomain<Scalar>::addSubSysMatrix, Kaa, Kaj, Kai, multiplier, isK);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::addSubSysMatrix(int iSub, const double *Kaa, const int *Kaj,
                                         const int *Kai, Scalar multiplier, bool isK)
{
  SubDomain(iSub)->addSysMatrix(neq, Kaa, Kaj, Kai, multiplier, sysMatrix[iSub], sparseK[iSub], isK);
}

template<class Scalar>
GenFetiDPSolver<Scalar> *
GenSandiaDomain<Scalar>::getSolver()
{
  Rbm **rbmPointer = 0;
  computeRbms = (computeRbms || (problem_type == 0 && this->numDualMpc > 0)); // PJSA 9-20-07

  if(have_contact) { finfo->useMRHS = false; finfo->nlPrecFlg = 0; } // disable multiple LHS and multiple RHS acceleration

  fetiSolver = new GenFetiDPSolver<Scalar>(this->numSub, this->subDomain, this->subToSub, finfo, this->communicator,
    this->glSubToLocal, this->mpcToSub_dual, this->mpcToSub_primal, this->mpcToMpc, this->mpcToCpu, this->cpuToSub, this->grToSub, 
    sysMatrix, sparseK, rbmPointer, problem_type+1, computeRbms); 

  // compute the geometric gap
  if(have_contact && finfo->geometric_gap) {
    if(this->myCPU == 0) cerr << "computing geometric gap\n";
    if(!initial_disp) { 
      initial_disp = new double[neq];
      for(int i=0; i<neq; ++i) initial_disp[i] = 0.0;
    }
    int index = 0;
    for(int i = 0; i < num_nodes; ++i) {
      int flags = int(dofmap_on_node[i]);
      DofSet dofs(flags);
      // unmark constrained dofs
      for(int j=0; j<ndofs_deleted; ++j) {
        if(bc_node_ids[j] == i) { 
          int flag = (1 << bc_node_dof[j]);
          dofs.unmark(flag);
        }
      }
      if(dofs.contains(DofSet::Xdisp)) initial_disp[index++] += x[i]; 
      if(dofs.contains(DofSet::Ydisp)) initial_disp[index++] += y[i];
      if(dofs.contains(DofSet::Zdisp)) initial_disp[index++] += z[i];
    }
  }
  if(initial_disp) { setContactGap(initial_disp); delete [] initial_disp; initial_disp = 0; }

  return (GenFetiDPSolver<Scalar> *) fetiSolver;
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::updateSysMatrixInSolver()
{
  fetiSolver->setSysMatrices(sysMatrix, sparseK);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::getRBMs(Scalar *salinasRBMs)
{
  int nRBM = fetiSolver->numRBM();

  int size  = 0;
  for(int i=0; i<this->numSub; ++i) size += this->subDomain[i]->numUncon();
  Scalar *fetiRBMs = new Scalar[nRBM*size];

  fetiSolver->getRBMs(fetiRBMs);
  paralApply(this->numSub, this->subDomain, &GenSandiaSubD<Scalar>::getRBMs, size, nRBM, salinasRBMs, fetiRBMs);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::makeF(Scalar *f, GenDistrVector<Scalar> *df)
{
  execParal(this->numSub, this, &GenSandiaDomain<Scalar>::makeSubF, f, df);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::makeSubF(int iSub, Scalar *f, GenDistrVector<Scalar> *df)
{
  SubDomain(iSub)->makeF(neq, f, df->subData(iSub));
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::getD(Scalar *dis, GenDistrVector<Scalar> *ddis, bool assemble)
{
  if(assemble) for(int i=0; i<neq; ++i) dis[i] = 0.0;
  execParal(this->numSub, this, &GenSandiaDomain<Scalar>::getSubD, dis, ddis, assemble);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::getSubD(int iSub, Scalar *dis, GenDistrVector<Scalar> *ddis, bool assemble)
{
  SubDomain(iSub)->getD(neq, dis, ddis->subData(iSub), assemble);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::constraintProduct(int num_vect, const double* R[],
                                           Scalar** V, int trans)
{
  execParal(this->numSub, this, &GenSandiaDomain<Scalar>::subConstraintProduct, 
            num_vect, R, V, trans);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::subConstraintProduct(int iSub, int num_vect, const double* R[],
                                              Scalar** V, int trans)
{
  cerr << "GenSandiaDomain<Scalar>::subConstraintProduct(...) is not implemented \n";
  exit(-1);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::mpcForces(double *lambdas)
{
  execParal(this->numSub, this, &GenSandiaDomain<Scalar>::subMpcForces, lambdas);
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::subMpcForces(int iSub, double *lambdas)
{
  double *local_lambdas = (this->numDualMpc) ? new double[this->subDomain[iSub]->numMPCs()] : new double[this->subDomain[iSub]->numMPCs_primal()];
  fetiSolver->getLocalMpcForces(iSub, local_lambdas);
  int m = 0;
  for(int i=0; i<numberMpc; i++) {
    if(mpclocal[i].NumEntries() > 0) {
      int locMpcNb = (this->numDualMpc) ? this->subDomain[iSub]->globalToLocalMPC[i] : this->subDomain[iSub]->globalToLocalMPC_primal[i];
      if(locMpcNb > -1) {
        lambdas[m] = local_lambdas[locMpcNb];
      }
      m++;
    }
  }
  if(local_lambdas) delete [] local_lambdas;
}

//#define MORTAR_DEBUG

template<class Scalar>
void
GenSandiaDomain<Scalar>::setSurfaceInteractions(int itype, int numNeighb, int *neighb, int *neighbPtr,
                                                enum face_type *faceEl_type, int *faceEl_ptr, int *faceEl_conn,
                                                double (*faceEl_normal)[3], double (*gap_vector)[3], int *glbl_surface_id,
                                                double *normal_tol, double *tangent_tol)
{
#ifdef MORTAR_DEBUG
  for(int icpu=0; icpu<this->numCPU; icpu++) {
    this->communicator->sync();
    if(icpu==this->myCPU) {
      cerr << " CPU " <<this->myCPU<<", numNeighb = "<<numNeighb<<endl; cerr.flush();
      for(int i=0; i<numNeighb; ++i) {
        cerr<< " -- i = " << i << ", neighb["<<i<<"] = "<<neighb[i];
        cerr<< ", glbl_surface_id["<<i<<"] = "<<glbl_surface_id[i];
        cerr<< ", normal_tol["<<i<<"] = "<<normal_tol[i]<<", tangent_tol["<<i<<"] = "<<tangent_tol[i]<<endl; cerr.flush();
        for(int j=neighbPtr[i]; j<neighbPtr[i+1]; ++j) {
          cerr<<"    -> j = "<<j<<", faceEl_type["<<j<<"] = "<<faceEl_type[j]<<endl;
          cerr<<"       normal = ("<<faceEl_normal[j][0]<<","<<faceEl_normal[j][1]<<","<<faceEl_normal[j][2]<<")\n";
          if(faceEl_type[j] != POINT_FACE)
            cerr<<"       gap = ("<<gap_vector[j][0]<<","<<gap_vector[j][1]<<","<<gap_vector[j][2]<<")\n";
          for(int k=faceEl_ptr[j]; k<faceEl_ptr[j+1]; ++k) {
            cerr<<"      -> k = "<<j<<", faceEl_conn["<<k<<"] = "<<faceEl_conn[k]<<" -> global Id = "
                <<global_node_nums[faceEl_conn[k]]<<endl;
            cerr.flush();
          }
        }
      }
    }
    this->communicator->sync();
  }
#endif
  if(itype == 1) have_contact = true;
  // currently mixing nodeToNode and mortar data is not supported
  // first check if there are any POINT_FACEs in faceEl_type, if so then assume all interactions are node-to-node
  bool nn = setNodeToNodeInteractions(itype, numNeighb, neighb, neighbPtr, faceEl_type, faceEl_ptr, faceEl_conn, faceEl_normal, gap_vector);
  if(!nn) setMortarInteractions(itype, numNeighb, neighb, neighbPtr, faceEl_type, faceEl_ptr, faceEl_conn, glbl_surface_id, normal_tol, tangent_tol);
}

template<class Scalar>
bool 
GenSandiaDomain<Scalar>::setNodeToNodeInteractions(int itype, int numNeighb, int *neighb, int *neighbPtr,
                                                   enum face_type *faceEl_type, int *faceEl_ptr, int *faceEl_conn,
                                                   double (*faceEl_normal)[3], double (*gap_vector)[3])
{
  // note: currently salinas is not passing the gap for node-to-node. It is passed later on in setContactGap

  // make cpuToDec connectivity XXXX this should only have the cpus with nodeToNode interations
  int *ptr = new int[this->numCPU+1];
  int *tgt = new int[this->numCPU];
  for(int i = 0; i < this->numCPU; ++i) ptr[i] = tgt[i] = i;
  ptr[this->numCPU] = this->numCPU;
  Connectivity *cpuToDec = new Connectivity(this->numCPU, ptr, tgt);

  // need to have a unique global id for each mpc. first step is to count the number of master nodes on each cpu and send to all
  int *master_count = (int*) dbg_alloca((this->numCPU)*sizeof(int));
  for(int i=0; i<this->numCPU; ++i) master_count[i] = 0;
  for(int i=0; i<numNeighb; ++i) {
    if(this->myCPU > neighb[i]) {
      for(int j=neighbPtr[i]; j<neighbPtr[i+1]; ++j)
        if(faceEl_type[j] == POINT_FACE) master_count[this->myCPU]++;
    }
  }
  this->communicator->globalSum(this->numCPU,master_count);
  int pair_offset = 0; for(int i=0; i<this->myCPU; ++i) pair_offset += master_count[i];
  int pair_total = 0; for(int i=0; i<this->numCPU; ++i) pair_total += master_count[i];

  if(pair_total == 0) { delete cpuToDec; return false; } // no node-to-node interations

  FSCommPattern<int> *pat = new FSCommPattern<int>(this->communicator, cpuToDec, this->myCPU, FSCommPattern<int>::CopyOnSend);
  FSCommPattern<double> *pat2 = new FSCommPattern<double>(this->communicator, cpuToDec, this->myCPU, FSCommPattern<double>::CopyOnSend);
  int *len = new int[this->numCPU];
  for(int i=0; i<this->numCPU; ++i) len[i] = 0;
  for(int i=0; i<numNeighb; ++i) 
    for(int j=neighbPtr[i]; j<neighbPtr[i+1]; ++j)
      if(faceEl_type[j] == POINT_FACE) {
        if(this->myCPU == neighb[i]) {
          cerr << " *** ERROR on CPU " << this->myCPU << ": node-to-node contact within a subdomain is NOT supported\n";
          exit(-1);
        }
        if((faceEl_ptr[j+1]-faceEl_ptr[j]) > 1)  { 
          cerr << " *** ERROR on CPU " << this->myCPU << ": more than one node per faceEl ("  << j << "," << faceEl_ptr[j+1]-faceEl_ptr[j] 
               << "). Only face type POINT_FACE is supported\n";
          exit(-1);
        }
        len[neighb[i]]++;
      } 
  for(int i=0; i<this->numCPU; ++i)
    if(this->myCPU != i) { pat->setLen(this->myCPU, i, len[i]); pat2->setLen(this->myCPU, i, 3*len[i]); }
  delete [] len;
  pat->finalize(); pat2->finalize();

  // send the pair_id from master to slave and the normals
  int pair_count = 0;
  for(int i=0; i<numNeighb; ++i) {
    if(this->myCPU != neighb[i]) {
      FSSubRecInfo<int> sInfo = pat->getSendBuffer(this->myCPU, neighb[i]);
      FSSubRecInfo<double> sInfo2 = pat2->getSendBuffer(this->myCPU, neighb[i]);
      int index = 0;
      for(int j=neighbPtr[i]; j<neighbPtr[i+1]; ++j)
        if(faceEl_type[j] == POINT_FACE) {
          sInfo.data[index] = (this->myCPU > neighb[i]) ? pair_offset+pair_count++ : -1;
          for(int k=0; k<3; ++k) sInfo2.data[3*index+k] = faceEl_normal[j][k];
          index++;
        }
    }
  }    
  pat->exchange(); pat2->exchange();
  
  // add new LMPCons objects to Domain::lmpc array
  int nNodeToNodeLMPCs = (itype == 1) ? pair_total : finfo->spaceDimension*pair_total;
  //if(this->myCPU == 0) cerr << "generating " <<  nNodeToNodeLMPCs << " LMPCs for node-to-node interaction\n";
  for(int i=numberMpc; i<numberMpc+nNodeToNodeLMPCs; ++i) {
    LMPCons *lmpc = new LMPCons(i,0.0);
    lmpc->type = itype;
    domain->addLMPC(lmpc);
  }
  pair_count = 0;
  for(int i=0; i<numNeighb; ++i) {
    if(this->myCPU != neighb[i]) {
      FSSubRecInfo<int> rInfo = pat->recData(neighb[i], this->myCPU);
      FSSubRecInfo<double> rInfo2 = pat2->recData(neighb[i], this->myCPU);
      int index = 0;
      for(int j=neighbPtr[i]; j<neighbPtr[i+1]; ++j) {
        if(faceEl_type[j] == POINT_FACE) {
          bool master = (this->myCPU > neighb[i]);
          int pair_id = (master) ? pair_offset+pair_count++ : rInfo.data[index];
          int lmpcnum = (itype == 1) ? numberMpc + pair_id : numberMpc + finfo->spaceDimension*pair_id;
          int n1 = (master) ? global_node_nums[faceEl_conn[faceEl_ptr[j]]] : -1;
          int n2 = (master) ? -1 : global_node_nums[faceEl_conn[faceEl_ptr[j]]];
          double gap[3] = { 0.0, 0.0, 0.0 }; // gaps are not passed at this time
          double sign = (master) ? 1.0 : -1.0;
          double normal[3]; for(int k=0; k<3; ++k) normal[k] = sign*(faceEl_normal[j][k] - rInfo2.data[3*index+k])/2.0; // average of master and slave
          domain->addNodeToNodeLMPCs(lmpcnum, n1, n2, normal, gap, itype);
          index++;
        }
      }
    }
  }
  delete pat; delete pat2; delete cpuToDec;

#ifdef MORTAR_DEBUG
  for(int icpu=0; icpu<this->numCPU; icpu++) {
    this->communicator->sync();
    if(icpu==this->myCPU) {
      cerr << " CPU " <<this->myCPU<<", numberMpc = "<<numberMpc<<", pair_offset = "<<pair_offset<<", pair_total = "<<pair_total<<", pair_count = "<<pair_count<<endl; cerr.flush();
      cerr << "master_count = "; for(int i=0; i<this->numCPU; ++i) cerr << master_count[i] << " "; cerr << endl; cerr.flush();
      cerr << "oneSubPerCPU = " << oneSubPerCPU << ", numSub = " << this->numSub << endl;
      domain->printLMPC2();
    }
    this->communicator->sync();
  }
#endif

  // update mpclocal and numberMpc
  augmentMpclocal();

  return true;
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::setMortarInteractions(int itype, int numNeighb, int *neighb, int *neighbPtr,
                                               enum face_type *faceEl_type, int *faceEl_ptr, int *faceEl_conn,
                                               int *glbl_surface_id, double *normal_tol, double *tangent_tol)
{
/*
  if(itype == 1) { // save pointers to this data so I can rebuild the mortar lmpcs later on using updated geometry
    _numNeighb = numNeighb; _neighb = neighb; _neighbPtr = neighbPtr;
    _faceEl_type = faceEl_type; _faceEl_ptr = faceEl_ptr; _faceEl_conn = faceEl_conn;
    _glbl_surface_id = glbl_surface_id; _normal_tol = normal_tol; _tangent_tol = tangent_tol;
  }
*/

  // make cpuToDec connectivity
  int *ptr = new int[this->numCPU+1];
  int *tgt = new int[this->numCPU];
  for(int i = 0; i < this->numCPU; ++i) ptr[i] = tgt[i] = i;
  ptr[this->numCPU] = this->numCPU;
  Connectivity *cpuToDec = new Connectivity(this->numCPU, ptr, tgt);

  // count number of surfaces to shared with each cpu
  int *surfCount = (int*) dbg_alloca((this->numCPU)*sizeof(int));
  for(int i=0; i<this->numCPU; ++i) surfCount[i] = 0;
  for(int i=0; i<numNeighb; ++i) surfCount[neighb[i]]++;

  // Step 1. mpi exchange & broadcast to find number of mortars (nMortarCond) 
  // and glbl_surface_id pairs for each mortar
  int i, /*j,*/ k;
  FSCommPattern<int> *sPat = new FSCommPattern<int>(this->communicator, cpuToDec, this->myCPU, FSCommPattern<int>::CopyOnSend);
  for(i=0; i<this->numCPU; ++i) {
    if(surfCount[i] == 0) continue;
    if(this->myCPU != i) { // PJSA 
      sPat->setLen(this->myCPU, i, surfCount[i]);
    }
  }
  sPat->finalize();
  int max_glbl_surface_id = 0;
  for(int i=0; i<this->numCPU; ++i) surfCount[i] = 0; // reset counter
  for(i=0; i<numNeighb; ++i) {
    if(this->myCPU != neighb[i]) { // PJSA 
      FSSubRecInfo<int> sInfo = sPat->getSendBuffer(this->myCPU, neighb[i]);
      sInfo.data[surfCount[neighb[i]]] = glbl_surface_id[i];
      surfCount[neighb[i]]++;
    }
    if(glbl_surface_id[i] > max_glbl_surface_id) max_glbl_surface_id = glbl_surface_id[i];
  }
  max_glbl_surface_id = this->communicator->globalMax(max_glbl_surface_id);
  int *neighb_glbl_surface_id = (int*) dbg_alloca(numNeighb*sizeof(int));
  sPat->exchange();
  for(int i=0; i<this->numCPU; ++i) surfCount[i] = 0; // reset counter
  bool first=true;
  for(i=0; i<numNeighb; ++i) {
    if(this->myCPU != neighb[i]) { // PJSA 
      FSSubRecInfo<int> rInfo = sPat->recData(neighb[i], this->myCPU);
      neighb_glbl_surface_id[i] = rInfo.data[surfCount[neighb[i]]];
      surfCount[neighb[i]]++;
    }
    else { 
      if(first) { 
        neighb_glbl_surface_id[i] = glbl_surface_id[i+1]; 
        first = false;
      }
      else { 
        neighb_glbl_surface_id[i] = glbl_surface_id[i-1];
        first = true;
      }
    }
  }
  delete sPat;
  int *surface_pair = (int*) dbg_alloca((max_glbl_surface_id+1)*sizeof(int));
  int *surface_pair_count = (int*) dbg_alloca((max_glbl_surface_id+1)*sizeof(int));
  for(i=0; i<=max_glbl_surface_id; ++i) {
    surface_pair[i] = 0;
    surface_pair_count[i] = 0;
  }

  for(i=0; i<numNeighb; ++i) {
    surface_pair[glbl_surface_id[i]] += neighb_glbl_surface_id[i];
    surface_pair_count[glbl_surface_id[i]]++;
  }
  this->communicator->globalSum((max_glbl_surface_id+1), surface_pair_count);
  this->communicator->globalSum((max_glbl_surface_id+1), surface_pair);
  int zeroCount = 0;
  for(i=0; i<=max_glbl_surface_id; ++i) {
    if(surface_pair_count[i] == 0)
      zeroCount++;
    else surface_pair[i] /= surface_pair_count[i];
  }

  // 2. mpi global broadcast & sum to build SurfEntity object for all glbl_surface_ids on every processor
  // (global node numbers and element types and node coordinates)
  SurfaceEntity **surfEntities = new SurfaceEntity * [max_glbl_surface_id+1];
  for(i=0; i<=max_glbl_surface_id; ++i) surfEntities[i] = new SurfaceEntity(i);
  bool *surfDone = (bool*) dbg_alloca((max_glbl_surface_id+1)*sizeof(bool)); 
  int numNodes[7] = {4, 8, 3, 6, 4, 2, 1};
  for(i=0; i<=max_glbl_surface_id; ++i) {
    if(surface_pair_count[i] == 0) continue;
    int *node_stride = (int*) dbg_alloca(this->numCPU*sizeof(int));
    int node_send_data[1] = { 0 };
    int *elem_stride = (int*) dbg_alloca(this->numCPU*sizeof(int));
    int elem_send_data[1] = { 0 };
    for(int j=0; j<=max_glbl_surface_id; ++j) surfDone[j] = false; 
    for(int j=0; j<numNeighb; ++j)
      if((glbl_surface_id[j] == i) && !surfDone[glbl_surface_id[j]]) {  
        node_send_data[0] += (faceEl_ptr[neighbPtr[j+1]] - faceEl_ptr[neighbPtr[j]]);
        elem_send_data[0] += (neighbPtr[j+1] - neighbPtr[j]);
        surfDone[glbl_surface_id[j]] = true; 
      }
    this->communicator->allGather(node_send_data,1,node_stride,1);
    this->communicator->allGather(elem_send_data,1,elem_stride,1);
    int surface_numElem = this->communicator->globalSum(elem_send_data[0]);
    int *node_displs = (int*) dbg_alloca(this->numCPU*sizeof(int));
    int *elem_displs = (int*) dbg_alloca(this->numCPU*sizeof(int));
    node_displs[0] = elem_displs[0] = 0;
    for(int j=1; j<this->numCPU; ++j) { 
      node_displs[j] = node_displs[j-1] + node_stride[j-1];
      elem_displs[j] = elem_displs[j-1] + elem_stride[j-1];
    }
    int node_bufsize = node_displs[this->numCPU-1] + node_stride[this->numCPU-1];
    int elem_bufsize = elem_displs[this->numCPU-1] + elem_stride[this->numCPU-1];
    //int *surfaceNodes = (int*) dbg_alloca(node_bufsize*sizeof(int));
    int *surfaceNodes = new int[node_bufsize];
    //double *xCoords = (double*) dbg_alloca(node_bufsize*sizeof(double));
    double *xCoords = new double[node_bufsize];
    //double *yCoords = (double*) dbg_alloca(node_bufsize*sizeof(double));
    double *yCoords = new double[node_bufsize];
    //double *zCoords = (double*) dbg_alloca(node_bufsize*sizeof(double));
    double *zCoords = new double[node_bufsize];
    //int *surfaceElemTypes = (int*) dbg_alloca(elem_bufsize*sizeof(int));
    int *surfaceElemTypes = new int[elem_bufsize];
    //int *surfaceNodesTmp = (int*) dbg_alloca(node_stride[this->myCPU]*sizeof(int));
    int *surfaceNodesTmp = new int[node_stride[this->myCPU]];
    //double *xCoordsTmp = (double*) dbg_alloca(node_stride[this->myCPU]*sizeof(double));
    double *xCoordsTmp = new double[node_stride[this->myCPU]];
    //double *yCoordsTmp = (double*) dbg_alloca(node_stride[this->myCPU]*sizeof(double));
    double *yCoordsTmp = new double[node_stride[this->myCPU]];
    //double *zCoordsTmp = (double*) dbg_alloca(node_stride[this->myCPU]*sizeof(double));
    double *zCoordsTmp = new double[node_stride[this->myCPU]];
    //int *surfaceElemTypesTmp = (int*) dbg_alloca(elem_stride[this->myCPU]*sizeof(int));
    int *surfaceElemTypesTmp = new int[elem_stride[this->myCPU]];
    int node_index = 0;
    int elem_index = 0;
    for(int j=0; j<=max_glbl_surface_id; ++j) surfDone[j] = false; 
    for(int j=0; j<numNeighb; ++j) {
      if((glbl_surface_id[j] == i) && !surfDone[glbl_surface_id[j]]) { 
        for(int k=faceEl_ptr[neighbPtr[j]]; k<faceEl_ptr[neighbPtr[j+1]]; ++k) {
          int locNode = faceEl_conn[k];
          surfaceNodesTmp[node_index] = global_node_nums[locNode];
          xCoordsTmp[node_index] = x[locNode];
          yCoordsTmp[node_index] = y[locNode];
          zCoordsTmp[node_index] = z[locNode];
          node_index++;
        }
        for(int k=neighbPtr[j]; k<neighbPtr[j+1]; ++k) {
          surfaceElemTypesTmp[elem_index] = faceEl_type[k];
          elem_index++;
        }
        surfDone[glbl_surface_id[j]] = true;
      }
    }
    this->communicator->allGatherv(surfaceNodesTmp, node_stride[this->myCPU], surfaceNodes, node_stride, node_displs);
    this->communicator->allGatherv(xCoordsTmp, node_stride[this->myCPU], xCoords, node_stride, node_displs);
    this->communicator->allGatherv(yCoordsTmp, node_stride[this->myCPU], yCoords, node_stride, node_displs);
    this->communicator->allGatherv(zCoordsTmp, node_stride[this->myCPU], zCoords, node_stride, node_displs);
    this->communicator->allGatherv(surfaceElemTypesTmp, elem_stride[this->myCPU], surfaceElemTypes, elem_stride, elem_displs);
 
    std::map<int,Node> *nodeCoordMap = new std::map<int,Node>();
    node_index = 0;
    for(int j=0; j<surface_numElem; ++j) {
      int nnodes = numNodes[surfaceElemTypes[j]];
      int *nodes = (int*) dbg_alloca(nnodes*sizeof(int));
      for(k=0; k<nnodes; ++k) {
        nodes[k] = surfaceNodes[node_index];
        double xyz[3] = { xCoords[node_index], yCoords[node_index], zCoords[node_index] };
        nodeCoordMap->insert(std::map<int, Node>::value_type(nodes[k], Node(xyz)));
        node_index++;
      }
      surfEntities[i]->AddFaceElement(j, surfaceElemTypes[j]+1, nnodes, nodes);
    }
    delete [] surfaceNodes;
    delete [] xCoords;
    delete [] yCoords;
    delete [] zCoords;
    delete [] surfaceElemTypes;
    delete [] surfaceNodesTmp;
    delete [] xCoordsTmp;
    delete [] yCoordsTmp;
    delete [] zCoordsTmp;
    delete [] surfaceElemTypesTmp;
#ifdef MORTAR_DEBUG
    if(this->myCPU==0) {
      cerr << " TOPOLOGY OF SURF "<<i<<" BEFORE RENUMBERING\n"; cerr.flush();
      surfEntities[i]->PrintFaceElemSet();
    }
    this->communicator->sync();
#endif
    surfEntities[i]->SetUpData();
    surfEntities[i]->MakeNodeSet((*nodeCoordMap));
    surfEntities[i]->Renumber();
#ifdef MORTAR_DEBUG
    if(this->myCPU==0) {
      cerr << " DATA OF SURF "<<i<<" AFTER RENUMBERING\n"; cerr.flush();
      surfEntities[i]->Print();
      filePrint(stderr," ------------------------------------------------------------------------\n");
      filePrint(stderr,"  average normal of face element of surface %2d AFTER local renumbering\n",surfEntities[i]->ID());
      filePrint(stderr," ------------------------------------------------------------------------\n");
      surfEntities[i]->PrintFaceNormal(surfEntities[i]->GetNodeSet());
    }
    this->communicator->sync();
#endif
    delete nodeCoordMap;
  }

  // broadcast the search tolerances
  double *all_normal_tols = new double[max_glbl_surface_id+1];
  double *all_tangent_tols = new double[max_glbl_surface_id+1];
  for(i=0; i<=max_glbl_surface_id; ++i) all_normal_tols[i] = all_tangent_tols[i] = 0.0;
  for(i=0; i<numNeighb; ++i) {
    all_normal_tols[glbl_surface_id[i]] = normal_tol[i];
    all_tangent_tols[glbl_surface_id[i]] = tangent_tol[i];
  }
  this->communicator->globalMax(max_glbl_surface_id+1,all_normal_tols);
  this->communicator->globalMax(max_glbl_surface_id+1,all_tangent_tols);

  for(i=0; i<=max_glbl_surface_id; ++i) {
    // note: surface with lower numbered global id is the master, higher number is slave
    if((surface_pair_count[i] != 0) && (i < surface_pair[i])) {
      int masterId = i, slaveId = surface_pair[i];
      MortarHandler *mortarCond = new MortarHandler(masterId, slaveId);
      mortarCond->SetInteractionType(itype); // 0 for tied, 1 for contact, 2 for fluid-structure
      mortarCond->SetSearchTol(all_normal_tols[i],all_tangent_tols[i]); // PJSA 9-7-2007
      if(mortar_type == 2) mortarCond->SetMortarType(MortarHandler::DUAL);
      mortarCond->SetPtrSurfaceEntity(surfEntities[masterId], surfEntities[slaveId]);
      domain->AddMortarCond(mortarCond);
    }
  }
  delete [] all_normal_tols; delete [] all_tangent_tols;

  if(domain->GetnMortarConds() > 0) {
    if(dofmap_on_node[0] & DofSet::Helm) { // XXXX assume mortars apply only to pressure dofs, not x,y,z displacements
      int mdofs[1] = { 6 };
      domain->ComputeMortarLMPC(1, mdofs);
    }
    else domain->ComputeMortarLMPC();
    if(domain->GetnMortarLMPCs() > 0) {
      domain->CreateMortarToMPC();
      //domain->normalizeLMPC();
      augmentMpclocal();
    }
    domain->DeleteMortarConds();
  }
  //domain->printLMPC();
  
  delete cpuToDec;
  if(surfEntities) {
    for(int i=0; i<max_glbl_surface_id+1; ++i) 
      if(surfEntities[i]) delete surfEntities[i];
    delete [] surfEntities;
  }
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::augmentMpclocal()
{
  // add new LMPCons from Domain::lmpc to salinas mpclocal array
  int numberMpc_copy = numberMpc;
  numberMpc = domain->getNumLMPC();
  ResizeArray<LMPCons *> *lmpc = domain->getLMPC();

  // make global to local node map
  if(!local_node_nums) { // make global to local node map
    global_node_max = 0;
    for(int i = 0; i < num_nodes; ++i)
      if(global_node_nums[i] > global_node_max) global_node_max = global_node_nums[i];
    local_node_nums = new int[global_node_max+1];
    for(int i=0; i<global_node_max; ++i) local_node_nums[i] = -1;
    for(int i = 0; i < num_nodes; ++i) local_node_nums[global_node_nums[i]] = i;
  }

  MpcLocal *mpclocal_new = new MpcLocal[numberMpc];
  for(int i = 0; i < numberMpc_copy; ++i)  {
    mpclocal_new[i] = mpclocal[i]; // this should copy data
  }
  for(int i=numberMpc_copy; i < numberMpc; ++i) {
    int numLocalEntries = 0;
    for(int j = 0; j < (*lmpc)[i]->nterms; ++j) {
      int glNode = (*lmpc)[i]->terms[j].nnum;
      int locNode = (glNode > -1 && glNode <= global_node_max) ? local_node_nums[glNode] : -1;
      if(locNode != -1) numLocalEntries++;
    }
    mpclocal_new[i].NumEntries(numLocalEntries);
    mpclocal_new[i].NumEntriesGlobal((*lmpc)[i]->nterms);

    int m = 0;
    for(int j = 0; j < (*lmpc)[i]->nterms; ++j) {
      int glNode = (*lmpc)[i]->terms[j].nnum;
      int locNode = (glNode > -1 && glNode <= global_node_max) ? local_node_nums[glNode] : -1;
      if(locNode != -1) {
        mpclocal_new[i].LocalId(m, locNode);
        mpclocal_new[i].Coef(m, (*lmpc)[i]->terms[j].coef.r_value);
        mpclocal_new[i].NodalDof(m, (*lmpc)[i]->terms[j].dofnum);
        mpclocal_new[i].LocalToGlobalEntry(m,j);
        mpclocal_new[i].G((*lmpc)[i]->rhs.r_value);
        m++;
      }
    }
  }
  delete [] mpclocal;
  mpclocal = mpclocal_new;
}


template<class Scalar>
void
GenSandiaDomain<Scalar>::findNumGlobNodes()
{
  // XXXX
  // Compute maximum node number in the global model
  // it would be nice if this wasn't necessary !!! Currently used in
  // GenDistrDomain<Scalar>::finishFETIDPInterface() for corner weighting
  int cpuNMax = 0;
  for(int i=0; i<this->numSub; ++i) {
    int subNMax = this->subDomain[i]->getGlobalNMax();
    if(subNMax > cpuNMax) cpuNMax = subNMax;
  }
  int domainNMax = structCom->globalMax(cpuNMax);
  geoSource->setNumGlobNodes(domainNMax+1);
  domain->setNumNodes(domainNMax+1); // PJSA 9-25-06
}

template<class Scalar>
int
GenSandiaDomain<Scalar>::setContactNormal(const double (*normals)[3])
{
  if(this->myCPU == 0) cerr << "GenSandiaDomain<Scalar>::setContactNormal(...) is not implemented\n";
  return 0;
}

template<class Scalar>
int
//GenSandiaDomain<Scalar>::setContactGap(const double (*gapVectors)[3])
GenSandiaDomain<Scalar>::setContactGap(const double *disp)
{
  finfo->type = FetiInfo::nonlinear;  // currently this function is only called for nonlinear, but 
                                      // perhaps it should also be called for linear analysis to set the initial gap

/*
  if(!position_vec) { // first time this is called compute the initial geometry used to compute the initial gap for contact
                      // non-zero disp on first call represents an "artificial gap" or initial displacement
    position_vec = new double[neq];
    for(int i=0; i<neq; ++i) position_vec[i] = 0.0;
    dx = new double[num_nodes]; dy = new double[num_nodes]; dz = new double[num_nodes];
    for(int i=0; i<num_nodes; ++i) dx[i] = dy[i] = dz[i] = 0.0;
    int index = 0;
    for(int i = 0; i < num_nodes; ++i) {
      int flags = int(dofmap_on_node[i]);
      DofSet dofs(flags);
      // unmark constrained dofs
      for(int j=0; j<ndofs_deleted; ++j) {
        if(bc_node_ids[j] == i) { 
          int flag = (1 << bc_node_dof[j]);
          dofs.unmark(flag);
        }
      }
      if(dofs.contains(DofSet::Xdisp)) { 
        position_vec[index] = (finfo->geometric_gap) ? x[i]+disp[index] : disp[index]; 
        dx[i] = disp[index]; 
        index++; 
      }
      if(dofs.contains(DofSet::Ydisp)) { 
        position_vec[index] = (finfo->geometric_gap) ? y[i] + disp[index] : disp[index]; 
        dy[i] = disp[index]; 
        index++; 
      }
      if(dofs.contains(DofSet::Zdisp)) { 
        position_vec[index] = (finfo->geometric_gap) ? z[i] + disp[index] : disp[index]; 
        dz[i] = disp[index]; 
        index++; 
      }
    }
  }
*/
  if(fetiSolver) fetiSolver->initNewton(); // this needs to be called before the start of each newton loop in nonlinear

  if(numberMpc > 0) {
    if(fetiSolver) {
      Scalar *current_position = new Scalar[neq];
      //for(int i=0; i<neq; ++i) current_position[i] = (disp) ? -(disp[i] + position_vec[i]) : -position_vec[i];
      for(int i=0; i<neq; ++i) current_position[i] = -disp[i];

      GenDistrVector<Scalar> *cu = new GenDistrVector<Scalar>(fetiSolver->interfInfo());
      GenDistrVector<Scalar> *u = new GenDistrVector<Scalar>(fetiSolver->localInfo());
      makeF(current_position, u); // map to feti numbering
      fetiSolver->multC(*u, *cu); // interfvec = C*u
      execParal(this->numSub, this, &GenSandiaDomain<Scalar>::setMpcRhs, *cu);
      delete cu; delete u;

      delete [] current_position;
    }
    else { // save the initial disp, it may be non-zero (eg for artifical gap)
      initial_disp = new double[neq];
      for(int i=0; i<neq; ++i) initial_disp[i] = disp[i];
    }
  }
  return 0;
}

template<class Scalar>
int
GenSandiaDomain<Scalar>::getContactForces(double *ctcForces)
{
  if(fetiSolver && numberMpc) {
    GenDistrVector<Scalar> *cf = new GenDistrVector<Scalar>(fetiSolver->localInfo());
    GenDistrVector<Scalar> *lambda_total = fetiSolver->getLambdaTotal();
    fetiSolver->trMultC(*lambda_total, *cf); // cf = C^T*lambda_total
  
    Scalar *tmp = new Scalar[neq];
    getD(tmp, cf, true); // assemble and map to salinas numbering 
    for(int i=0; i<neq; ++i) ctcForces[i] = -ScalarTypes::Real(tmp[i]);
    delete [] tmp;  delete cf;
  }
  else for(int i=0; i<neq; ++i) ctcForces[i] = 0.0; // first time
  return 0;
}

template<class Scalar>
int
GenSandiaDomain<Scalar>::updateContactGap(const double *DU)
{
  if(fetiSolver && numberMpc) {
    GenDistrVector<Scalar> *cu = new GenDistrVector<Scalar>(fetiSolver->interfInfo());
    GenDistrVector<Scalar> *u = new GenDistrVector<Scalar>(fetiSolver->localInfo());
    Scalar *tmp = new Scalar[neq];
    for(int i=0; i<neq; ++i) tmp[i] = -DU[i];
    makeF(tmp, u); // map to feti numbering
    fetiSolver->multC(*u, *cu); // interfvec = C*u
    execParal(this->numSub, this, &GenSandiaDomain<Scalar>::updateMpcRhs, *cu);
    delete [] tmp; delete cu; delete u;
  }
  return 0;
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::updateMpcRhs(int iSub, GenDistrVector<Scalar> &cu)
{
  SubDomain(iSub)->updateMpcRhs(cu.subData(SubDomain(iSub)->localSubNum()));
}


template<class Scalar>
void
GenSandiaDomain<Scalar>::setMpcRhs(int iSub, GenDistrVector<Scalar> &cu)
{
  SubDomain(iSub)->setMpcRhs(cu.subData(SubDomain(iSub)->localSubNum()));
}

template<class Scalar>
void
GenSandiaDomain<Scalar>::zeroMpcRhs(int iSub)
{
  SubDomain(iSub)->zeroMpcRhs();
}

