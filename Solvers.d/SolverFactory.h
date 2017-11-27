#ifndef __SOLVER_FACTORY_H__
#define __SOLVER_FACTORY_H__

#include <Solvers.d/Solver.h>
#include <Solvers.d/SolverCntl.h>
#include <map>
#include <memory>
#include <string>

class EqNumberer;
class DofSetArray;
class ConstrainedDSA;
template <class Scalar> class GenSparseMatrix;
class FSCommunicator;

template <typename Scalar>
struct SolverAndMatrix {
	GenSolver<Scalar> *solver;
	std::unique_ptr<GenSparseMatrix<Scalar>> sparseMatrix;
};

template<class Scalar>
class GenSolverFactory
{
 public:
  virtual ~GenSolverFactory() {}

  // used for feti coarse solvers: sequential or parallel direct solver only
  virtual GenSolver<Scalar>* createSolver(Connectivity *con, EqNumberer *eqnum, SolverCntl&, GenSparseMatrix<Scalar> *&sparse, int ngrbm = 0,
                                          FSCommunicator *com = 0, std::string name = "") const;
  virtual GenSolver<Scalar>* createDistSolver(Connectivity *con, EqNumberer *eqnum, SolverCntl&, GenSparseMatrix<Scalar> *&sparse, 
                                              FSCommunicator *com = 0, std::string name = "") const;
  // used for global solvers and feti local solvers: direct or iterative solver
  virtual GenSolver<Scalar>* createSolver(Connectivity *con, DofSetArray *dsa, ConstrainedDSA *cdsa, SolverCntl&, GenSparseMatrix<Scalar> *&sparse, Rbm *rbm,
                                          GenSparseMatrix<Scalar> *&spp, GenSolver<Scalar> *&prec, FSCommunicator *com = 0, std::string name = "") const;
	// used for feti Kii solver: sequential direct solver only
	virtual SolverAndMatrix<Scalar> createSolver(Connectivity *con, DofSetArray *dsa, int *map, SolverCntl&,
	                                             std::string name = "") const;

  static GenSolverFactory* getFactory();

};

extern std::unique_ptr<GenSolverFactory<double> >   solverFactory;
extern std::unique_ptr<GenSolverFactory<DComplex> > solverFactoryC;

#ifdef _TEMPLATE_FIX_
  #include <Solvers.d/SolverFactory.C>
#endif

#endif

