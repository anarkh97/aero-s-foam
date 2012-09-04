#ifndef SALINAS_CU_Feti_PARAMS_H
#define SALINAS_CU_Feti_PARAMS_H

#include <iostream>

#define PRECOND_LUMPED     1
#define PRECOND_DIRICHLET  2
#define RBM_GEOMETRIC      1
#define RBM_ALGEBRAIC      2

enum CORNER_AUGMENTATION { no_augmentation, subdomain_augmentation, 
                           edge_augmentation, weighted_edge_augmentation };

enum OUTERLOOP_SOLVER { cg, gmres, gcr, cgal }; // PJSA 6-15-07

enum CORNER_AUG_RBM_TYPE { translation, all, none };

enum SOLVER_TYPE { Skyline, Sparse, Spooles, Auto };

enum WEIGHTING_TYPE { Topological, Stiffness };

// MPC Method options
struct MPC_Method {
  enum Method { Dual, Primal } method; // : 8;
  enum SubMethod { None, Full, PerFace, PerSub, Diag, BlockDiag, Auto } submethod; // : 8;
  enum SOLVER_TYPE solver; // : 8;
  enum WEIGHTING_TYPE weighting; // : 8;
  double cct_tol;  // tolerance used to factor CC^T matrix for dual MPCs
};

struct ContactParams {
  double eq_tol, ieq_tol; // tolerances for equality and inequality constraints
  double gamma; // proportioning tolerance
  double rho;   // initial penalty parameter (for augmented lagrangian)
  double beta;  // penalty update factor (for augmented lagrangian)
  double M, eta; // precision control parameters (for augmented lagrangian)
  int cgal_prec; // 0 = augmented lagrangian precondititoning off, 1 = on
  double alphabar; // used for 2-step expansion, currently used for cgal
  bool rebuild_cct; // PJSA 9-19-07
  bool geometric_gap; // PJSA 9-20-07
};
  
class FetiParams {
private:  // declare, don't define
    FetiParams( const FetiParams&);
    FetiParams& operator=( const FetiParams&);

private:
  int max_iterations;
  double solver_tol;      // controls the outer iteration loop
  double stagnation_tol;  // controls stagnation 
  int max_N_orthog_vecs;  // this is for PCCG. set to 1 for no orthog 
  int precondition;       // precondition method. bitwise OR of 
                          // lumped | dirichlet
  int rbm_method;         // method for finding rigid body modes. geometric or
                          // algebraic. Note that if boundary conditions are applied prior to
                          // calling Feti, then only algebraic methods can be used.
  double rbm_tol_svd;     // RBM tolerance within the SVD. geometric method only 
  double rbm_tol_mech;    // RBM mechanical tolerance. either method
  double wave_ortho_tol;  // tolerance for Gram-Schmidt orthog of the wave 
                          // Q Vector (default 1e-2)
  int feti_level;         // 1 = FETI-1
                          // 2 = FETI-DP
                          // 3 = FETI-DPC
                          // 4 = FETI-DPH
  int projector;          // 0 std, 1 = Q projector
  int solution_type;      // 0 = statics
                          // 1 = eigen
                          // 2 = implicit transient
                          // 3 = direct frequency response
  int bailOut;            // controls what happens when stagnation or max iterations occurs
                          // during a feti solve.  default is to exit if either occurs before
                          // convergence, and to continue if bailout option is set.
  unsigned int verbose_flag;
                          // signals level of printed output.
                          // 1st bit: 1=do, (0=do not) print summary timer info
                          // 2nd bit: 1=do print # rbm in each subdomain
                          // 3rd and 4th bits: debug info
                          // other bits currently undefined
  enum SOLVER_TYPE local_solver;         // indicates solver to use for local ludcmp
  enum SOLVER_TYPE precondition_solver;  // indicates solver to use for preconditioner
  enum SOLVER_TYPE coarse_solver;        // indicates solver to use for primary coarse problem (Kcc* in DP)
  enum SOLVER_TYPE aux_coarse_solver;    // indicates solver to use for auxillary coarse problem (GtG in DP)
  int pivoting;           // 0 = no pivoting, 1 = pivoting 
                          // (only applies to solvers that support pivoting, eg. spooles
  enum OUTERLOOP_SOLVER outerloop_solver;
  bool multiple_rhs;
  unsigned int coarse_skip;
                          // 1 indicates skip computation of coarse problem
  double grbm_tol;        // the relative tolerance used to determine which
                          // small pivots in the LDL' factorization of G'FG
                          // correspond to global rigid body modes.
  double crbm_tol;        // relative tolerance Kcc factorization
  int corner_dimensionality;  // 3 or 6 dofs/corner
  int corner_algorithm;   // 1=std_old, 3=3 per neighbor, 5, 6, 7, 8 use
                          // nQ = 1, 2, 3 & 4 respectively
                          // where nQ+2 (here) is # of touching subdomains

  enum CORNER_AUGMENTATION corner_augmentation;
  enum CORNER_AUG_RBM_TYPE corner_aug_rbm_type;
  enum WEIGHTING_TYPE weighting;

  struct MPC_Method mpc_method;    
  int numWaveDirections;  // number of wave directions used to constuct Q matrix in FETI-DPH
                          // range 0-13, default = 3
  double waveNumbers[3];  // waveNumbers[0] = pressure wave number for this subdomain in FETI-DPH
                          // waveNumbers[1] = 1st shear wave number for this subdomain in FETI-DPH
                          // waveNumbers[2] = 2nd shear wave number (same as waveNumbers[1] for isotropic)
  int nlPrecFlg;          // nlPrecFlg  = 0  do not use Krylov preconditioning (multiple LHS acceleration)
                          // nlPrecFlg  = 1  use Krylov preconditioning
                          // nlPrecFlg  = 2  use Krylov preconditioning per load step
  int numLocalSubdomain;  // currently not used
  int helmholtz_type;	  // 0 for structural, 1 for acoustic, 2 for coupled
  int space_dimension;	  // 1 = 1-D, 2 = 2-D, 3 = 3-D (used for DPH augmentation only)
  int nrbms;
  double shift;
  struct ContactParams contact_params; // PJSA 6-15-07
  int initial_lambda; // PJSA 9-11-07
  int mortar_type; // PJSA 9-11-07: 1 = standard, 2 = dual

public:
    friend void readFetiparams();

// constructor
     FetiParams() {
       max_iterations        = 1000; 
       solver_tol            = 1.0e-08; 
       stagnation_tol        = 1.0e-06;   
       max_N_orthog_vecs     = 1000; 
       precondition          = PRECOND_DIRICHLET;
       rbm_method            = RBM_GEOMETRIC;
       rbm_tol_svd           = 1.0e-10;
       rbm_tol_mech          = 1.0e-8;
       wave_ortho_tol        = 1.0e-2;
       feti_level            = 2;
       projector             = 0;
       solution_type         = 0;
       verbose_flag          = 1;
       bailOut               = 0;
       local_solver          = Sparse;
       precondition_solver   = Sparse;
       coarse_solver         = Sparse;
       aux_coarse_solver     = Skyline;
       outerloop_solver      = cg;
       multiple_rhs          = true;
       coarse_skip           = 0;
       grbm_tol              = 1e-6;
       crbm_tol              = 1e-6;
       corner_dimensionality = 6;
       corner_algorithm      = 1;
       corner_augmentation   = edge_augmentation;
       corner_aug_rbm_type   = translation;
       mpc_method.method     = MPC_Method::Dual;
       mpc_method.submethod  = MPC_Method::Auto; // PJSA 9-21-07
       mpc_method.solver     = Skyline;
       mpc_method.weighting  = Topological;
       mpc_method.cct_tol    = 1.0e-8; // PJSA 6-15-07
       numWaveDirections     = 3;
       waveNumbers[0]        = 0.0;
       waveNumbers[1]        = 0.0;
       waveNumbers[2]        = 0.0;
       numLocalSubdomain     = 0;
       nlPrecFlg             = 0;
       weighting             = Topological; // PJSA 9-21-07
       pivoting              = 1;
       space_dimension       = 3;
       helmholtz_type        = 0;
       nrbms                 = 0;
       shift                 = 0.0;
       contact_params.eq_tol = 1.0e-6;
       contact_params.ieq_tol = 1.0e-6;
       contact_params.gamma  = 1.0;
       contact_params.rho    = 1.0;
       contact_params.beta   = 10.0;
       contact_params.M      = 1.0;
       contact_params.eta    = 1.0e-4;
       contact_params.cgal_prec = 0;
       contact_params.alphabar = 2.0;
       contact_params.rebuild_cct = true; // PJSA 9-19-07
       contact_params.geometric_gap = false; // PJSA 9-19-07
       initial_lambda         = 0; // PJSA 9-11-07
       mortar_type            = 1; // PJSA 9-11-07
     }

    // data access routines .................................
    int Max_iterations() const { return max_iterations; }
    double Solver_tol() const  { return solver_tol; }
    double Stagnation_tol() const { return stagnation_tol; }
    int Max_N_orthog_vecs() const { return max_N_orthog_vecs; }
    int Precondition() const { return precondition; }
    int Rbm_method() const { return rbm_method; }
    double Rbm_tol_svd() const { return rbm_tol_svd; }   
    double Rbm_tol_mech() const { return rbm_tol_mech; }  
    int Feti_level() const { return feti_level; }       
    int Projector() const { return projector; }        
    int Solution_type() const { return solution_type; }
    int BailOut() const { return bailOut; }          
    unsigned int Verbose_flag() const { return verbose_flag; }
    enum SOLVER_TYPE Local_solver() const { return local_solver; }
    enum SOLVER_TYPE Precondition_solver() const { return precondition_solver; }
    enum SOLVER_TYPE Coarse_solver() const { return coarse_solver; }
    enum SOLVER_TYPE Aux_Coarse_solver() const { return aux_coarse_solver; }
    int Pivoting() const { return pivoting; }

    unsigned int Coarse_skip() const { return coarse_skip; }
    double Grbm_tol() const { return grbm_tol; }
    double Crbm_tol() const { return crbm_tol; }
    int Corner_dimensionality() const { return corner_dimensionality; }  
    int Corner_algorithm() const { return corner_algorithm; }  
    enum CORNER_AUGMENTATION Corner_augmentation() const { return corner_augmentation; }
    enum CORNER_AUG_RBM_TYPE Corner_aug_rbm_type() const { return corner_aug_rbm_type; }
    MPC_Method Mpc_method() const { return mpc_method; }        
    int NumLocalSubdomain() const { return numLocalSubdomain; } 
    int NlPrecFlg() const { return nlPrecFlg; }
    enum WEIGHTING_TYPE Weighting() const { return weighting; }

    enum OUTERLOOP_SOLVER Outerloop_solver() const { return outerloop_solver; }
    bool Multiple_rhs() const { return multiple_rhs; }
    int NumWaveDirections() const { return numWaveDirections; }    
    double Wave_ortho_tol() const { return wave_ortho_tol; }
    double* WaveNumbers() { return waveNumbers; }
    int Helmholtz_type() { return helmholtz_type; }
    int Space_dimension() { return space_dimension; }
    int Nrbms() { return nrbms; } 
    double Shift() { return shift; }
    ContactParams Contact_params() const { return contact_params; } // PJSA 6-15-07
    int Initial_lambda() const { return initial_lambda; } // PJSA 9-11-07
    int Mortar_type() const { return mortar_type; } // PJSA 9-11-07

    // routines to set variables
    void Max_iterations(int m)            { max_iterations=m; }
    void Solver_tol(double t)             { solver_tol=t; }
    void Stagnation_tol(double t)         { stagnation_tol=t; }
    void Max_N_orthog_vecs(int m)         { max_N_orthog_vecs=m; }
    void Precondition(int m)              { precondition=m; }
    void Rbm_method(int m)                { rbm_method=m; }
    void Rbm_tol_svd(double t)            { rbm_tol_svd=t; }   
    void Rbm_tol_mech(double t )          { rbm_tol_mech=t; }  
    void Feti_level(int m)                { feti_level=m; }
    void Projector(int m)                 { projector=m; }        
    void Solution_type(int m)             { solution_type=m; }
    void BailOut(int m)                   { bailOut=m; }          
    void Verbose_flag(unsigned int n)     { verbose_flag=n; }
    void Local_solver(enum SOLVER_TYPE m) { local_solver=m; }
    void Precondition_solver(enum SOLVER_TYPE m) { precondition_solver=m; }
    void Coarse_solver(enum SOLVER_TYPE m) { coarse_solver=m; }
    void Aux_coarse_solver(enum SOLVER_TYPE m) { aux_coarse_solver=m; }
    void Pivoting(int m)                  { pivoting = m; }

    void Coarse_skip(unsigned int n)      { coarse_skip=n; }
    void Grbm_tol(double t)               { grbm_tol=t; }
    void Crbm_tol(double t)               { crbm_tol=t; }
    void Corner_dimensionality(int m)     { corner_dimensionality=m; }  
    void Corner_algorithm(int m)          { corner_algorithm=m; }  
    void Corner_augmentation(enum CORNER_AUGMENTATION a) 
	{ corner_augmentation=a; }
    void Corner_aug_rbm_type(enum CORNER_AUG_RBM_TYPE a) 
	{ corner_aug_rbm_type=a; }
    void Mpc_method(MPC_Method meth)      { mpc_method=meth; }        
    void NumLocalSubdomain(int m)         { numLocalSubdomain=m; } 
    void Weighting(enum WEIGHTING_TYPE m)     { weighting = m; }
    void NlPrecFlg(int m)                 { nlPrecFlg = m; }

    void Outerloop_solver(enum OUTERLOOP_SOLVER a) { outerloop_solver=a; }
    void Multiple_rhs(bool m) { multiple_rhs=m; }
    void NumWaveDirections(int m)         { numWaveDirections=m;   }    
    void Wave_ortho_tol(double t)         { wave_ortho_tol=t; }
    void WaveNumbers(const double* wn) {
         waveNumbers[0] = wn[0];
         waveNumbers[1] = wn[1];
         waveNumbers[2] = wn[2];
     }
    void Helmholtz_type(int i)		   { helmholtz_type = i; }
    void Space_dimension(int i)		   { space_dimension = i; }
    void Nrbms(int n)                      { nrbms = n; }
    void Shift(double s)                   { shift = s; }
    void Contact_params(ContactParams cp)  { contact_params = cp; } // PJSA 6-15-07
    void Initial_lambda(int i)             { initial_lambda = i; } // PJSA 9-11-07
    void Mortar_type(int m)                { mortar_type = m; } // PJSA 9-11-07

    void print() {
      cerr << "***** Feti Parameters *****\n";
      cerr << "max_iterations        = " << max_iterations << endl;
      cerr << "solver_tol            = " << solver_tol << endl;
      cerr << "stagnation_tol        = " << stagnation_tol << endl;
      cerr << "max_N_orthog_vecs     = " << max_N_orthog_vecs << endl;
      cerr << "precondition          = " << precondition << endl;
      cerr << "rbm_method            = " << rbm_method << endl;
      cerr << "rbm_tol_svd           = " << rbm_tol_svd << endl;
      cerr << "rbm_tol_mech          = " << rbm_tol_mech << endl;
      cerr << "wave_ortho_tol        = " << wave_ortho_tol << endl;
      cerr << "feti_level            = " << feti_level << endl;
      cerr << "projector             = " << projector << endl;
      cerr << "solution_type         = " << solution_type << endl;
      cerr << "verbose_flag          = " << verbose_flag << endl;
      cerr << "bailOut               = " << bailOut << endl;
      cerr << "local_solver          = " << local_solver << endl;
      cerr << "precondition_solver   = " << precondition_solver << endl;
      cerr << "coarse_solver         = " << coarse_solver << endl;
      cerr << "aux_coarse_solver     = " << aux_coarse_solver << endl;
      cerr << "pivoting              = " << pivoting << endl;
      cerr << "outerloop_solver      = " << outerloop_solver << endl;
      cerr << "multiple_rhs          = " << multiple_rhs << endl;
      cerr << "coarse_skip           = " << coarse_skip << endl;
      cerr << "grbm_tol              = " << grbm_tol << endl;
      cerr << "crbm_tol              = " << crbm_tol << endl;
      cerr << "corner_dimensionality = " << corner_dimensionality << endl;
      cerr << "corner_algorithm      = " << corner_algorithm << endl;
      cerr << "corner_augmentation   = " << corner_augmentation << endl;
      cerr << "corner_aug_rbm_type   = " << corner_aug_rbm_type << endl;
      cerr << "weighting             = " << weighting << endl;
      cerr << "mpc_method.method     = " << mpc_method.method << endl;
      cerr << "mpc_method.submethod  = " << mpc_method.submethod << endl;
      cerr << "mpc_method.solver     = " << mpc_method.solver << endl;
      cerr << "mpc_method.weighting  = " << mpc_method.weighting << endl;
      cerr << "mpc_method.cct_tol    = " << mpc_method.cct_tol << endl; // PJSA 9-19-07
      cerr << "numWaveDirections     = " << numWaveDirections << endl;
      cerr << "waveNumbers           = " << waveNumbers[0] << " " << waveNumbers[1]
                                         << " " << waveNumbers[2] << endl;
      cerr << "numLocalSubdomain     = " << numLocalSubdomain << endl;
      cerr << "nlPrecFlg             = " << nlPrecFlg << endl;
      cerr << "helmholtz_type        = " << helmholtz_type << endl;
      cerr << "space_dimension       = " << space_dimension << endl;
      cerr << "nrbms                 = " << nrbms << endl;
      cerr << "shift                 = " << shift << endl;
      cerr << "contact_params.eq_tol = " << contact_params.eq_tol << endl;
      cerr << "contact_params.ieq_tol = " << contact_params.ieq_tol << endl;
      cerr << "contact_params.gamma  = " << contact_params.gamma << endl;
      cerr << "contact_params.rho    = " << contact_params.rho << endl;
      cerr << "contact_params.beta   = " << contact_params.beta << endl;
      cerr << "contact_params.M      = " << contact_params.M << endl;
      cerr << "contact_params.eta    = " << contact_params.eta << endl;
      cerr << "contact_params.cgal_prec = " << contact_params.cgal_prec << endl;
      cerr << "contact_params.alphabar = " << contact_params.alphabar << endl;
      cerr << "contact_params.rebuild_cct = " << contact_params.rebuild_cct << endl;
      cerr << "contact_params.geometric_gap = " << contact_params.geometric_gap << endl;
      cerr << "initial_lambda        = " << initial_lambda << endl;
      cerr << "mortar_type           = " << mortar_type << endl;
      cerr << "***************************\n";
   }

}; 

#endif
