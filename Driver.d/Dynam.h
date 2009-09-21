#ifndef _DYNAM_H_
#define _DYNAM_H_

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenFetiSolver;
typedef GenFetiSolver<double> FetiSolver;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
class Rbm;

template <class Scalar>
class GenDynamMat {
 public:
   GenSolver<Scalar>    *dynMat;        // use to solve (coeM*M +coeC*C + coeK*K)x = b 
   GenSolver<Scalar>    *Msolver;       // use to solve Mx=b
   GenSparseMatrix<Scalar> *K;		// Stiffness Matrix
   GenSparseMatrix<Scalar> *refK;	// Stiffness Matrix for eigensolves
   GenSparseMatrix<Scalar> *C;		// Damping Matrix
   GenSparseMatrix<Scalar> *Cuc;	// constrained to unconstrained Damping Matrix
   GenSparseMatrix<Scalar> *M;		// Mass Matrix
   GenSparseMatrix<Scalar> *Muc;        // constrained to unconstrained Mass Matrix
   GenSparseMatrix<Scalar> *Mcc;        // constrained to constrained Mass Matrix
   GenSparseMatrix<Scalar> *kuc;
   int           numdofs;	// number of dof
   Rbm* rigidBodyModes;

   // Constructor
   GenDynamMat() { dynMat = 0; Msolver = 0; C = 0; M = 0; Cuc = 0; Muc = 0; Mcc = 0; refK = 0; kuc = 0; rigidBodyModes = 0; }
   GenDynamMat(GenDynamMat *d) { 
       dynMat = (*d).dynMat; Msolver = (*d).Msolver;
       K = (*d).K; refK = (*d).refK; 
       C = (*d).C; Cuc = (*d).Cuc; 
       M = (*d).M; Muc = (*d).Muc; Mcc = (*d).Mcc;
       numdofs = (*d).numdofs;
       rigidBodyModes = (*d).rigidBodyModes;
   }
};
typedef GenDynamMat<double> DynamMat;

class PitaDynamMat: public DynamMat{

 //CD: for Pita, solver on the coarse time grid
 public:
  Solver  *dynMat_Dt;

  PitaDynamMat() : DynamMat() { dynMat_Dt = 0; }
  PitaDynamMat(DynamMat *d) : DynamMat(d) { dynMat_Dt = 0;}

};


class EigenMat {
 public:
   FetiSolver *dynMat;
   SparseMatrix *M;             // Mass Matrix
   int           numdofs;       // number of dof
   Vector       *rbm;           // rigid Body Modes

   // Constructor
   EigenMat() { dynMat = 0; M = 0; rbm = 0; };

   // Member functions
   void getJacobi(double *kappa,double * mu, FullSquareMatrix &xx,
               double *eigVal,int nsmax,int subSpaceSize,double tolJac);
   void ortho(Vector *v1, Vector *vr, int nsub, int nrmod);
   void ortho(VectorSet& v1, VectorSet& vr, int nsub, int nrmod);

};

#endif
