#ifndef DOMAIN_GROUP_TASK_H_
#define DOMAIN_GROUP_TASK_H_

#include <Threads.d/Paral.h>

template <class Scalar> class GenSubDomain;
template <class Scalar> class GenSolver;
template <class Scalar> class GenSparseMatrix;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class Rbm;
class FSCommunicator;

template<class Scalar>
class GenDomainGroupTask : public TaskDescr {
 public:
   int nsub;
   GenSubDomain<Scalar> **sd;
   GenSolver<Scalar> **dynMats;
   GenSparseMatrix<Scalar> **spMats;
   GenSparseMatrix<Scalar> **M;
   GenSparseMatrix<Scalar> **Muc;
   GenSparseMatrix<Scalar> **Mcc;
   GenSparseMatrix<Scalar> **C;
   GenSparseMatrix<Scalar> **Cuc;
   GenSparseMatrix<Scalar> **Ccc;
// RT
   GenSparseMatrix<Scalar> ***C_deriv;
   GenSparseMatrix<Scalar> ***Cuc_deriv;
// RT end
   GenSparseMatrix<Scalar> **K;
   GenSparseMatrix<Scalar> **Kuc;
   GenSparseMatrix<Scalar> **spp;
   GenSolver<Scalar> **sps;
   Rbm **rbms; // geometric based RBMs
   FullSquareMatrix **kelArray, **melArray, **celArray;
   double coeM, coeC, coeK;
   double alpha, beta;
   int numSommer;
   int solvertype;
   FSCommunicator *com;
   bool makeC, makeC_deriv;

   GenDomainGroupTask(int nsub, GenSubDomain<Scalar> **_sd, double, double, double,
                      Rbm **_rbms, FullSquareMatrix **_kelArray, double, double, 
                      int, int solvertype, FSCommunicator *, FullSquareMatrix **_melArray,
                      FullSquareMatrix **_celArray, bool elemsetHasDamping);
   virtual ~GenDomainGroupTask();
   void runFor(int isub, bool make_feti);
};

typedef GenDomainGroupTask<double> DomainGroupTask;

#endif
