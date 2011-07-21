#ifndef _FETI_INFO_H_
#define _FETI_INFO_H_

// maxit    = maximum number of FETI iterations
// tol      = global tolerance
// maxortho = maximum number of reorthogonalization directions

// Preconditioner types:
// precno     = 0       no preconditioner
// precno     = 1       lumped			(default)
// precno     = 2       dirichlet

// Solver types for coarse problem: (GtG)x = b
// other solvers not implemented yet.
// gtgSolver  = 0       skyline                 (default)
// gtgSolver  = 1       sparse 
// gtgSolver  = 2       blocksky

// Solver types for Subdomain matrices
// solvertype = 0       skyline			(default)
// solvertype = 1       sparse
// solvertype = 2       sgiSparse
// solvertype = 3       sgiSkyline
// solvertype = 4       pcg
// solvertype = 5       frontal
// solvertype = 6       blocksky
// solvertype = 7       parallelSparse
// solvertype = 8       spooles
// solvertype = 9       mumps
// solvertype = 10      diagonal

// Projector types:
// nonLocalQ = 0	basic projector		(default)
// nonLocalQ = 1	Q(preconditioner type)  projector

// Scaling types:
// scaling = 1		kscaling (stiffness)	(default)
// scaling = 2		tscaling (topology)

// FETI Version
// version = 0		FETI 1			(default)
// version = 1		FETI 2

// Whether to solve coarse problem or not
// only valid for dynamics problems
// noCoarse = 0		coarse problem used	(default)
// noCoarse = 1		coarse problem not used

// Global rigid body mode relative tolerance
// grbm_tol = 1.0E-6 (default)

// Number of iterations to print the error
// printNumber = -1     no printing
// printNumber =  1	every iteration
// printNumber =  5     every 5 iterations
// printNumber =  n     every n iterations

// contactPrintFlag = 2         print status change and also other info
// contactPrintFlag = 1   	print status change info for FETI-DPC (x+- etc)
// contactPrintFlag = 0		don't print status change info (default)

// Nonlinear FETI information
// nTang = rebuild Tangent Stiffness matrix every N iterations
//         
// nPrec = rebuild FETI preconditioner every N iterations 
//         NOTE: if the user specifies to rebuild FETI every
//               n iterations and the preconditioner every N
//               iterations, then only rebuild preconditioner
//               when FETI is being rebuilt.

// Nonlinear Krylov information
// nlPrecFlg  = 0       do not use Krylov preconditioning
// nlPrecFlg  = 1       use Krylov preconditioning
// nlPrecFlg  = 2       use Krylov preconditioning per load step

// We need to keep track of the number of load steps while using FETI
//
// numLoadSteps = Number of Load Steps
//              = 0 to start and is incremented at new Load step

// Which corner degrees of freedom to clamp.
// corners = allCorners6    clamp all corner dofs (default)
//         = allCorners3
//         = noEndCorners3
//         = noEndCorners6

// FSI corner types: 
// fsi_corner = 0;    no corners on fluid/structure wet interface
// fsi_corner = 1;    all fluid wet subdomain interface nodes are corners 
// fsi_corner = 2;    all fluid and structure wet subdomain interface nodes are corners 
// fsi_corner = 3;    all fluid and a few structure wet subdomain interface nodes are corners

// f_projector = 0	don't project the force (default)
//	       = 1	use averaging projector for force
// 	       = 2	use trimming projector for force
// e_projector = 0	don't project estar (default)
//	       = 1	use averaging projector for estar
// u_projector = 0	don't project u
//	       = 1	project u (default)

// MPC Preconditioner types: (for Rixen method)
// mpc_precno     = 0       no preconditioner
// mpc_precno     = 1       diagonal                  
// mpc_precno     = 2       global 
// mpc_precno     = 3       topo block diagonal 
// mpc_precno     = 4       sub block diagonal
// mpc_precno     = 5       mortar edge block diagonal
// mpc_precno     = 6       auto select (default)

// MPC type
// mpcflag        = 0	    ignore mpcs
// mpcflag	  = 1  	    "dual" Rixen method (mpc lagrange multipliers) (default)
// mpcflag	  = 2	    "primal" include mpcs in coarse problem
// mpcflag        = 3       mixed dual and primal

// MPC Scaling types:
// mpc_scaling = 1          kscaling (stiffness) 
// mpc_scaling = 2          tscaling (topology) (default)

// Solver types for CC^t matrices
// cctSolver = 0       skyline                 (default)
// cctSolver = 1       sparse

//HB: Scaling for CC^t solver (currently ONLY supported for skyline solver)
// cctScaled = false    scaling OFF (default)
// cctScaled = true     scaling ON

//HB: overlap level for mortar block CC^t approximate solve 
// mpcBlkOverlap = 0 (default)

// For FETI-H,
// - Construction of the mass interface matrix : lumpedinterface
// - Number of directions for the coarse grid : numcgm

// For FETI-DPH
// outerloop = 0                    use CG solver
// outerloop = 1                    use GMRES solver
// outerloop = 2                    use GCR solver
// outerloop = 3                    use CGAL solver

// numdir    = 3 (default)          number of wave directions added in Q matrix
// orthotol  = 1.0E-02 (default)    relative tolerance value in orthogonalizing Q matrix
// orthotol2 = 1.0E-02 (default)    absolute tolerance value in orthogonalizing Q matrix

class FetiInfo {

  public:

    // Constructor
    FetiInfo();

    // Data members
    int    maxit;
    double tol;
    double absolute_tol;
    double stagnation_tol;
    double absolute_stagnation_tol;
    double grbm_tol;
    double crbm_tol;
    double cct_tol; // used to factorize CC^t in rixen mpc method
    int    uproj;
    int    maxortho;
    int    noCoarse;
    int    nonLocalQ;
    int    nQ;
    int    primalFlag; // whether to output primal residual
    int    printMatLab;

    int    printNumber;

    // Nonlinear Data members
    int    nPrec;
    int    nTang;
    int    nlPrecFlg;
    int    numLoadSteps;

    enum Preconditioner { noPrec, lumped, dirichlet, identity } precno;
    enum PreconditionerType { nonshifted, shifted } prectype;
    enum MpcPreconditioner { noMpcPrec=0, diagCCt, globalCCt, blockDiagCCt, subBlockDiagCCt, superBlockDiagCCt, autoSelectCCt } mpc_precno;
    enum MpcBlock { subBlock, topoBlock, mortarBlock } mpc_block;
    int mpcflag;
    enum Solvertype { skyline, sparse, sgisparse, sgisky, pcg, frontal,
         blocksky, parallel_sparse, spooles, mumps, diagonal } solvertype, gtgSolver, auxCoarseSolver, cctSolver;
    enum Scaling { noscaling=0, kscaling=1, tscaling=2 } scaling, mpc_scaling, fsi_scaling;
    enum Version { feti1, feti2, feti3, fetidp } version;
    bool rescalef; // if this is true then reassemble and apply scaling to f for every system, not just the first
                   // effects the convergence criteria for nonlinear and dynamics since relative primal error will be 
                   // defined using norm of rescaled f which is typically much smaller than the norm of the unscaled f
    enum Feti2Version { fullCoarse, sparseCoarse } feti2version;
    enum Type { linear, nonlinear, eigen } type;
    enum CornerType { allCorners6, allCorners3, noEndCorners6,
                      noEndCorners3, interface3, interface6, ThreeD, noCorners } corners;
    enum AugmentType { none, Gs, Edges, WeightedEdges } augment;
    int isEdgeAugmentationOn() { 
      return (((augment == Edges) || (augment == WeightedEdges)) && ((nGs > 0) || (numdir > 0))) ? 1 : 0; 
    }
    enum RbmType { translation, rotation, all,
                   averageTran, averageRot, averageAll, None,
                   pressure, temperature } rbmType;
    double nullSpaceFilterTol;

    // FETI-H
    double tolcgm;
    int numcgm; // number of coarse grid modes
    double numcgm2;
    int spaceDimension;
    int krylovtype;
    int lumpedinterface; // 0 - default (consistent) 1 - lumped
    int saveMemCoarse;

    int nGs;

    int maxiter()      { return maxit;       }
    int maxorth()      { return maxortho;    }
    int nPrecond()     { return nPrec;       } 
    double tolerance() { return tol;         }
    int numPrint()     { return printNumber; }

    enum OuterloopType { CG, GMRES, GCR, CGAL }  outerloop; 
    enum WaveType      { solid, shell, fluid, any } waveType; 
    enum WaveMethod    { averageK, averageMat, uniform } waveMethod;  // note: only use uniform if domain is homogeneous 
                                                                      // and entirely solid or shell or fluid
    int numdir;
    double orthotol, orthotol2; 
    bool dph_flag;
    int contactPrintFlag;
    bool cctScaled;
    int rebuildcct;
    int rebuildSbb;
    bool geometric_gap;
    int mpcBlkOverlap; //0=no interaction, 1=1st order interactions, 2=1st & 2nd order interactions, etc.
    double gamma;
    bool bmpc, dmpc, cmpc;
    double linesearch_tau;
    int linesearch_maxit;
    bool c_normalize;

    bool useMRHS;
    bool gmresResidual; //HB: to force computing the "primal residual" at each GMRES iteration;
    bool wetcorners;
    bool splitLocalFsi;
    int pickAnyCorner;
    bool pick_unsafe_corners; // if this is true (default) unsafe nodes will be selected as corners
                              // should only be set to false if pivoting is enabled for local solver

    bool fsi_element, mpc_element; // true means add to element set & decompose
    int fsi_corner;
    bool complex_hermitian;
    double dual_proj_tol, primal_proj_tol;
    double dual_plan_tol, primal_plan_tol;
    int dual_plan_maxit, primal_plan_maxit;
};

inline
FetiInfo::FetiInfo()
{
  maxit      = 1000;       // default max number of iterations
  tol        = 1.0E-6;     // default global tolerance
  absolute_tol = 0.0; 
  stagnation_tol = 1.0e-6;
  absolute_stagnation_tol = 0.0;
  grbm_tol = 1.0E-06;      // default global rigid body mode rel. tolerance
  crbm_tol = 1.0E-12;      // default global corner rigid body mode tolerance
  maxortho   = maxit;      // default max number of reortho vectors
                           // if it finds zero, set maxortho = maxit
  primalFlag  = 0;         // default do not output primal residual
  noCoarse    = 0;         // default use coarse problem
  precno      = dirichlet;    // default use lumped preconditioner
  prectype    = nonshifted;// default use nonshifted preconditioner (only the sitffness part) //HB
  solvertype  = sparse;   // default use skyline subdomain solver
  gtgSolver   = sparse;   // default use skyline GtG solver
  auxCoarseSolver = sparse;
  cctSolver   = sparse;   // default use sparse CCt solver
  nonLocalQ   = 0;         // default basic projector
  nQ          = 0;
  scaling     = tscaling;  // default use t scaling
  rescalef    = true;      // reassemble and apply scaling to f for every system, not just the first
  mpc_scaling = tscaling;  // default use t scaling
  fsi_scaling = tscaling;
  version     = feti1;     // default use FETI1
  feti2version= sparseCoarse;  // default use New FETI2
  printNumber = 10;        // default print error at every FETI iteration
  corners    = noEndCorners3; // default clamp all corner dofs
  augment    = Edges;       // default no Kcc augmentation
  nGs        = 6;
  rbmType    = all;
  gmresResidual = false;     // to force computing the "primal residual" at each GMRES iteration
  pickAnyCorner = 1;

  // Nonlinear information
  nTang        = 1;        // default always rebuild tangent
  nPrec        = 1;        // default always rebuild preconditioner
  nlPrecFlg    = 0;        // default Krylov preconditioner off
  type         = linear;   // default linear FETI problem
  numLoadSteps = 0;        // default 0 load step

  // defaults for FETI-H
  lumpedinterface = 0;
  numcgm       = 0;
  spaceDimension = 3;
  tolcgm = 1e-3;
  krylovtype = 8;
  saveMemCoarse = 0;

  outerloop    = CG;           // default CG solver 
  waveType     = any;      // default wave for solid/shell/fluid problem 
  waveMethod   = averageMat;  // k_p and k_s calculated using average material properties
  numdir       = 0;        // default number of wave directions added in Q matrix
  orthotol     = 1.0E-02;  // default relative tolerance value in orthogonalizing Q matrix
  orthotol2    = 0.0;      // default absolute value in orthogonalizing Q matrix
  dph_flag     = false;

  // DPC information
  uproj        = 1;   	   // default project the displacements wrt global RBMs
  contactPrintFlag = 0;
  gamma = 1.0;
  linesearch_tau = 0.6667; 
  linesearch_maxit = 100;
  cmpc = bmpc = dmpc = false;

  c_normalize = false;
  
  // MPC information
  mpc_precno   = globalCCt;  // default mpc preconditioner for Rixen method
  mpc_block    = topoBlock;
  mpcflag      = 1;
  cct_tol      = 1.0e-16;
  cctScaled    = false;      // XXXX no scaling used in the CCt solver (ONLY for skyline)
  mpcBlkOverlap= 0;         // zero/minimal overlap in mortar block CCt preconditionner
  rebuildcct   = 1; 
  rebuildSbb   = 0; 
  geometric_gap = false;

  // coupled dph 
  useMRHS = true;
  wetcorners   = false;
  fsi_scaling  = tscaling;
  splitLocalFsi= true;

  gmresResidual= false;     // to force computing the "primal residual" at each GMRES iteration
  printMatLab = 0;  

  mpc_element = fsi_element = false;
  pick_unsafe_corners = true;
  fsi_corner = 2;
  complex_hermitian = false;
  nullSpaceFilterTol = 0.0;
  dual_proj_tol = primal_proj_tol = 1.0e-16;
  dual_plan_tol = primal_plan_tol = 0.0;
  dual_plan_maxit = primal_plan_maxit = 20;
}


#endif
