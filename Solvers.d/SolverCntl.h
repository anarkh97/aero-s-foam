#ifndef _SOLVER_CNTL_
#define _SOLVER_CNTL_
#include <Feti.d/FetiInfo.h>
#include <map>

struct SolverCntl {
public:
   SolverCntl() { type = 0;
                  subtype = 1; // By default we use direct sparse solver
                  trbm = 1.0E-16;   // default zero pivot tolerance
                  trbm2 = 1.0E-16;  // default zero pivot tolerance
                  sparse_renum = 0;
                  sparse_maxsup = 100;
                  sparse_defblk = 30;
                  pivot = false;
                  unsymmetric = false;
                  scaled = false;
                  spooles_scale = 0;
                  spooles_tau = 100.;
                  spooles_seed = 532196;
                  spooles_maxsize = 64;
                  spooles_maxdomainsize = 24;
                  spooles_maxzeros = 0.04;
                  spooles_msglvl = 0;
                  spooles_renum = 0;
                  mumps_icntl[3] = 0; // supress diagnostic output
                  goldfarb_tol = 1.0;
                  goldfarb_check = false;
                  precond = 0;
                  tol = 1.0e-8;
                  maxit = 1000;
                  iterType = 0;
                  iterSubtype = 3;
                  maxvecsize = 0;
                  printMatLab        = false;
                  printMatLabFile    = "";
                  verbose = 1;
                  ilu_droptol = 1e-11;
   }
   int type;     // 0 = direct, 1 = iterative, 2 = FETI, 3 = Block Diag
   int subtype;  // subtype ... 9 is mumps  10 is diag
   double trbm;         // algebraic rbm tolerance
   double trbm2;        // algebraic rbm tolerance used for sparse/skyline when GRBM is activated
   int sparse_renum;  // renumbering scheme for BLKSparseMatrix: 0 = esmond MMD (default), 1 = metis ND
   int sparse_maxsup, sparse_defblk;
   bool pivot;  // true if pivoting is to be used in spooles/mumps solvers
   bool unsymmetric;
   bool scaled; // true if scaling is to be used in skyline solver
   int spooles_scale; // true if scaling is to be used in spooles solver
   double spooles_tau;  // used when pivoting is enabled, all entries in L and U have magnitude
                        // less than or equal to tau, default is 100.
   double spooles_maxzeros; // see Solvers.d/Spooles.C for description
   int spooles_maxsize, spooles_maxdomainsize, spooles_seed, spooles_msglvl; // see Solvers.d/Spooles.C for description
   int spooles_renum; // renumbering scheme for spooles: 0 = best of ND and MS, 1 = MMD, 2 = MS, 3 = ND
   std::map<int, int> mumps_icntl;
   std::map<int, double> mumps_cntl;
   double goldfarb_tol;
   bool goldfarb_check;
   int iterType; // 0 = CG, 1 = GMRES, 2 = GCR, 4 = BCG, 5 = CR
   int iterSubtype; // matrix storage
   int precond;  // preconditioner 0 = none, 1 = jacobi
   int maxit;    // maximum number of iterations
   double tol;   // tolerance for convergence
   int maxvecsize;  // for pcg # of krylov vectors to store default = 0
   FetiInfo fetiInfo;
   bool printMatLab;
   const char * printMatLabFile;
   int verbose;
   double ilu_droptol;
};

#endif
