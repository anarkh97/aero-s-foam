#ifndef DOMAIN_GROUP_TASK_H_
#define DOMAIN_GROUP_TASK_H_

#include <Math.d/SparseMatrix.h>
#include <Threads.d/Paral.h>

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
// RT
   GenSparseMatrix<Scalar> ***C_deriv;
   GenSparseMatrix<Scalar> ***Cuc_deriv;
// RT end
   GenSparseMatrix<Scalar> **K;
   GenSparseMatrix<Scalar> **Kuc;
   Rbm **rbms; // geometric based RBMs
   FullSquareMatrix **kelArray, **melArray;
   double coeM, coeC, coeK;
   double alpha, beta;
   int numSommer;
   int solvertype;
   FSCommunicator *com;

   GenDomainGroupTask(int nsub, GenSubDomain<Scalar> **_sd, double, double, double,
                      Rbm **_rbms, FullSquareMatrix **_kelArray, double, double, 
                      int, int solvertype, FSCommunicator *, FullSquareMatrix **_melArray);
   virtual ~GenDomainGroupTask();
   void runFor(int isub, bool make_feti);
};

typedef GenDomainGroupTask<double> DomainGroupTask;

#ifdef _TEMPLATE_FIX_
#include <Paral.d/DomainGroupTask.C>
#endif

#endif
