//
// Created by Michel Lesoinne on 2019-01-14.
//

#ifdef DISTRIBUTED
#include <Driver.d/Communicator.h>
#endif
#include "PardisoSolver.h"

#ifdef WITH_PARDISO
#include <complex>

#include <Math.d/DBSparseMatrix.h>
#include <mkl_pardiso.h>
#include <mkl_types.h>

namespace New {

template <typename idx_type>
using matLocPredicate = bool (*)(idx_type, idx_type);

template <typename idx_type>
struct DofForDofRange {
	idx_type *ndBegin;
	idx_type *ndEnd;

	matLocPredicate<idx_type> predicate;
};

struct CERange {

};

class SparseData {
public:
	using idx_type = int;
	SparseData(const Connectivity &cn, const EqNumberer &eqn, matLocPredicate<idx_type> row_col_select);

	idx_type numCol() const { return ptr.size()-1; };

	gsl::span<const idx_type> indices(idx_type r) const {
		return { targets.data() + ptr[r], targets.data() + ptr[r+1] };
	}

	idx_type row_offset(idx_type r) const {
		return ptr[ r ];
	}

	auto &rawPtr() const { return ptr; }
	auto &rawTargets() const { return targets; }

	template <typename I>
	std::vector<I> offsetCastPtr(I offset) const;
	template <typename I>
	std::vector<I> offsetCastTargets(I offset) const;
private:
	idx_type neq;
	std::vector<idx_type> ptr;
	std::vector<idx_type> targets;
};

SparseData::SparseData(const Connectivity &cn, const EqNumberer &eqn, matLocPredicate<idx_type> predicate)
{
	neq = static_cast<idx_type>(eqn.size());
	ptr.assign(neq+1, 0);
	for (idx_type rowNode = 0; rowNode < cn.csize(); ++rowNode)
	{
		idx_type firstRowDof = eqn.firstdof(rowNode);
		if(firstRowDof < 0)
			continue;
		idx_type lastRowDof = firstRowDof + eqn.weight(rowNode);

		for (auto colNode: cn[ rowNode ]) {
			idx_type firstColDof = eqn.firstdof(colNode);
			if(firstColDof < 0)
				continue;
			idx_type lastColDof = firstColDof + eqn.weight(colNode);
			for (auto r = firstRowDof; r < lastRowDof; ++r)
				for (auto c = firstColDof; c < lastColDof; ++c)
					if( (*predicate)(r,c))
						++ptr[r];
		}
	}
	for (idx_type i = 1; i <= neq; ++i)
		ptr[i] += ptr[i-1];
	targets.resize(ptr[neq]);
	for (idx_type rowNode = 0; rowNode < cn.csize(); ++rowNode)
	{
		idx_type firstRowDof = eqn.firstdof(rowNode);
		idx_type lastRowDof = firstRowDof + eqn.weight(rowNode);
		for (auto colNode: cn[ rowNode ]) {
			idx_type firstColDof = eqn.firstdof(colNode);
			idx_type lastColDof = firstColDof + eqn.weight(colNode);
			for (auto r = firstRowDof; r < lastRowDof; ++r)
				for (auto c = firstColDof; c < lastColDof; ++c)
					if( (*predicate)(r,c))
						targets[ --ptr[r] ] = c;
		}
	}
	for (idx_type i = 0; i < neq; ++i)
		std::sort(&targets[ptr[i]], &targets[ptr[i+1]]);
}

template<typename I>
std::vector<I> SparseData::offsetCastPtr(I offset) const
{
	std::vector<I> res;
	res.reserve(ptr.size());
	for( auto p : ptr )
		res.push_back( p + offset );
	return res;
}

template<typename I>
std::vector<I> SparseData::offsetCastTargets(I offset) const
{
	std::vector<I> res;
	res.reserve(targets.size());
	for( auto p : targets )
		res.push_back( p + offset );
	return res;
}

template <typename Scalar>
class DBSparseMatrix : public GenSparseMatrix<Scalar> {
public:
	using idx_type = int;
	DBSparseMatrix(const Connectivity &cn, const EqNumberer &eqn, matLocPredicate<idx_type> predicate);

	/// Move to the solver?
	int dim() const { return sparseData.numCol(); }
	int numCol() const { return sparseData.numCol(); }

	void add(const GenAssembledFullM<Scalar> &kel, gsl::span<const int> dofs);

	void add(const GenAssembledFullM<Scalar> &kel, const int *dofs) override {
		add(kel, gsl::span<const int>{ dofs, (int)kel.dim() });
	}

	void add(const FullSquareMatrix &, const int *dofs) override { throw "Not good"; }

	void unify(FSCommunicator *communicator);

	void zeroAll() override { a.assign(a.size(), 0); }


	Scalar diag( idx_type i ) const override { return a[ sparseData.row_offset(i)]; }

	Scalar &diag( idx_type i ) override { return a[ sparseData.row_offset(i)]; }

	std::vector<idx_type> getFPtr() const {
		return sparseData.offsetCastPtr<idx_type>(1);
	}

	std::vector<idx_type> getFTg() const {
		return sparseData.offsetCastTargets<idx_type>(1);
	}

	const std::vector<idx_type> &getCPtr() const {
		return sparseData.rawPtr();
	}

	const std::vector<idx_type> &getCTg() const {
		return sparseData.rawTargets();
	}

	auto indices(idx_type r) const { return sparseData.indices(r); }
protected:
	Scalar *getA() { return a.data(); }

	std::vector<Scalar>  &getK() { return a; }
private:
	matLocPredicate<idx_type> predicate;
	SparseData sparseData;
	std::vector<Scalar> a;
};

template<typename Scalar>
DBSparseMatrix<Scalar>::DBSparseMatrix(const Connectivity &cn, const EqNumberer &eqn, matLocPredicate<int> predicate)
	: predicate(predicate), sparseData(cn, eqn, predicate)
{
	a.assign(sparseData.rawPtr().back(), 0);
}

template<typename Scalar>
void DBSparseMatrix<Scalar>::add(const GenAssembledFullM<Scalar> &kel, gsl::span<const int> dofs)
{
	for (int i = 0; i < dofs.size() ; ++i) {                     // Loop over rows.
		auto r = dofs[i];
		if ( r < 0 ) continue;                   // Skip constrained dofs
		auto indices = sparseData.indices( r );
		Scalar *lc = a.data() + sparseData.row_offset( r );
		for ( int j = 0; j < dofs.size() ; ++j) {                       // Loop over columns.
			auto c = dofs[j];
			if( c < 0) continue;                     // Skip irrelevant dofs
			if ( !(*predicate)( r, c ) )
				continue;                 // Work with approved matrix part
			for (int m = 0; m < indices.size(); ++m) {
				if( indices[m] == c ) {
					lc[m] += kel[i][j];
					break;
				}
			}
		}
	}
}



template<typename Scalar>
void DBSparseMatrix<Scalar>::unify(FSCommunicator *communicator)
{
#ifdef DISTRIBUTED
	communicator->globalSum(a.size(), a.data());
#endif
}

}

template <typename Scalar>
class PardisoSolver :
	public New::DBSparseMatrix<Scalar>, public GenSolver<Scalar> {
public:
	PardisoSolver(const Connectivity *cn, const EqNumberer *_dsa);

	~PardisoSolver() override;

	int neqs() const override;

	void unify(FSCommunicator *communicator) override {
		New::DBSparseMatrix<Scalar>::unify(communicator);
	}

	void solve(const Scalar *rhs, Scalar *solution) override;

	void reSolve(Scalar *rhs) override;

	long size() const override;

	void factor() override;

	void parallelFactor() override;
private:
	MKL_INT n;
	std::vector<MKL_INT> ptr;
	std::vector<MKL_INT> tg;
//	/* Pardiso memory pointers, internal data. */
	void *pt[64];
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, error, msglvl;
	MKL_INT mtype;
};

template<typename Scalar>
int PardisoSolver<Scalar>::neqs() const
{
	return New::DBSparseMatrix<Scalar>::numCol();
}

template<typename Scalar>
long PardisoSolver<Scalar>::size() const
{
	// TODO get Pardiso's memory usage.
	return 0;
}
void test() {
			/* Matrix data. */
		MKL_INT n = 8;
		MKL_INT ia[9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19};
		MKL_INT ja[18] =
			{ 1,    3,       6, 7,
			  2, 3,    5,
			  3,             8,
			  4,       7,
			  5, 6, 7,
			  6,    8,
			  7,
			  8
			};
		double a[18] =
			{ 7.0,      1.0,           2.0, 7.0,
			  -4.0, 8.0,      2.0,
			  1.0,                     5.0,
			  7.0,           9.0,
			  5.0, 1.0, 5.0,
			  -1.0,      5.0,
			  11.0,
			  5.0
			};
		MKL_INT mtype = -2;       /* Real symmetric matrix */
		/* RHS and solution vectors. */
		double b[8], x[8];
		MKL_INT nrhs = 1;     /* Number of right hand sides. */
		/* Internal solver memory pointer pt, */
		/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
		/* or void *pt[64] should be OK on both architectures */
		void *pt[64];
		/* Pardiso control parameters. */
		MKL_INT iparm[64];
		MKL_INT maxfct, mnum, phase, error, msglvl;
		/* Auxiliary variables. */
		MKL_INT i;
		double ddum;          /* Double dummy */
		MKL_INT idum;         /* Integer dummy. */
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
		for ( i = 0; i < 64; i++ )
		{
			iparm[i] = 0;
		}
		iparm[0] = 1;         /* No solver default */
		iparm[1] = 2;         /* Fill-in reordering from METIS */
		iparm[3] = 0;         /* No iterative-direct algorithm */
		iparm[4] = 0;         /* No user fill-in reducing permutation */
		iparm[5] = 0;         /* Write solution into x */
		iparm[6] = 0;         /* Not in use */
		iparm[7] = 2;         /* Max numbers of iterative refinement steps */
		iparm[8] = 0;         /* Not in use */
		iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
		iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
		iparm[11] = 0;        /* Not in use */
		iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
		iparm[13] = 0;        /* Output: Number of perturbed pivots */
		iparm[14] = 0;        /* Not in use */
		iparm[15] = 0;        /* Not in use */
		iparm[16] = 0;        /* Not in use */
		iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
		iparm[18] = -1;       /* Output: Mflops for LU factorization */
		iparm[19] = 0;        /* Output: Numbers of CG Iterations */
		maxfct = 1;           /* Maximum number of numerical factorizations. */
		mnum = 1;         /* Which factorization to use. */
		msglvl = 1;           /* Print statistical information in file */
		error = 0;            /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
		for ( i = 0; i < 64; i++ )
		{
			pt[i] = 0;
		}
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
		phase = 11;
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		         &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if ( error != 0 )
		{
			printf ("\nERROR during symbolic factorization: %d", error);
			exit (1);
		}
		printf ("\nReordering completed ... ");
		printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
		printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
		phase = 22;
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		         &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if ( error != 0 )
		{
			printf ("\nERROR during numerical factorization: %d", error);
			exit (2);
		}
		printf ("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
		phase = 33;
		iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
		/* Set right hand side to one. */
		for ( i = 0; i < n; i++ )
		{
			b[i] = 1;
		}
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		         &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
		if ( error != 0 )
		{
			printf ("\nERROR during solution: %d", error);
			exit (3);
		}
		printf ("\nSolve completed ... ");
		printf ("\nThe solution of the system is: ");
		for ( i = 0; i < n; i++ )
		{
			printf ("\n x [%d] = % f", i, x[i]);
		}
		printf ("\n");
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
		phase = -1;           /* Release internal memory. */
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		         &n, &ddum, ia, ja, &idum, &nrhs,
		         iparm, &msglvl, &ddum, &ddum, &error);
}
template<typename Scalar>
PardisoSolver<Scalar>::PardisoSolver(const Connectivity *cn, const EqNumberer *_dsa)
	: New::DBSparseMatrix<Scalar>(*cn, *_dsa, [](int r, int c) { return c >= r; })
{

	n = static_cast<MKL_INT>(this->numCol());
	for (int i = 0; i < 64; ++i)
		pt[i] = nullptr;
	for(int i = 0; i < 64;++i)
		iparm[i] = 0;
	iparm[0] = 1;         /* No solver default */
	iparm[1] = 2;         /* Fill-in reordering from METIS */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* Permutation is internal. */
	iparm[5] = 0;         /* returns the solution in a separate vector. */
	iparm[6] = 0;         /* Not in use */
	iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	iparm[8] = 0;         /* Not in use */
	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 0;        /* Don't Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;        /* Not in use */
	iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[14] = 0;        /* Not in use */
	iparm[15] = 0;        /* Not in use */
	iparm[16] = 0;        /* Not in use */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	iparm[34] = 1;        /* C-Style indexing */

	iparm[26] = 1; // Check that the matrix format is OK.
	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;             /* Which factorization to use. */
	msglvl = 1;           /* Print statistical information in file */
	error = 0;            /* Initialize error flag */
	mtype = -2;       /* Real symmetric indefinite matrix */
	MKL_INT nrhs = 1;     /* Number of right hand sides. */
	double ddum;          /* Double dummy */
	MKL_INT idum;         /* Integer dummy. */

	ptr = this->getCPtr();
	tg = this->getCTg();
	for(int i = 0; i < 16 && i < ptr.size()-1; ++i) {
		std::cout << i << ": ( " << ptr[i+1]-ptr[i] <<" )";
		for(int j = ptr[i]; j < ptr[i+1]; ++j)
			std::cout << " " << tg [ j ];
		std::cout << std::endl;
//		auto indices = this->indices(i);
//		std::cout << i << ": ( " << indices.size() << " )";
//		for( auto idx: indices )
//			std::cout << " " << idx;
//		std::cout << std::endl;
	}


//	std::cout << "ptr[n-1]: " << ptr[n-1] << " ptr[n]: " << ptr[n] << std::endl;

/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
	MKL_INT phase = 11;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	         &n, this->getA(), ptr.data(), tg.data(), nullptr, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if ( error != 0 )
	{
		printf ("\nERROR during symbolic factorization: %d", error);
		exit (1);
	}
	printf ("\nReordering completed ... ");
	printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf ("\nNumber of factorization MFLOPS = %d\n", iparm[18]);
}

template<typename Scalar>
void PardisoSolver<Scalar>::solve(const Scalar *rhs, Scalar *solution)
{
	MKL_INT nrhs = 1;     /* Number of right hand sides. */
	MKL_INT phase = 33;
	iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
	iparm[6] = 0; // Stor the result

	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	         &n, this->getA(), ptr.data(), tg.data(), nullptr, &nrhs, iparm, &msglvl, const_cast<Scalar*>(rhs), &solution, &error);
	if ( error != 0 )
	{
		printf ("ERROR during  solution: %d\n", error);
		exit (2);
	}
}

template<typename Scalar>
void PardisoSolver<Scalar>::reSolve(Scalar *rhs)
{
	MKL_INT nrhs = 1;     /* Number of right hand sides. */
	MKL_INT phase = 33;
	iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
	iparm[5] = 1; // Store the result in place of the input.

	std::vector<Scalar> dx(n,-100);
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	         &n, this->getA(), ptr.data(), tg.data(), nullptr, &nrhs, iparm, &msglvl, rhs, dx.data(), &error);
	if ( error != 0 )
	{
		printf ("ERROR during  solution: %d\n", error);
		exit (2);
	}
}

template<typename Scalar>
void PardisoSolver<Scalar>::factor()
{
	MKL_INT nrhs = 1;     /* Number of right hand sides. */
	double ddum;          /* Double dummy */
	MKL_INT idum;         /* Integer dummy. */
	MKL_INT phase = 22;
	Scalar *a = this->getA();
	auto &K = this->getK();
	std::cout << "Size of K: " << K.size() <<" " << K[0] << " " << K[1] << " " << K[2] << std::endl;
	std::cout << "First 3 as: " << a[0] << " " << a[1] << " " << a[2] << std::endl;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	         &n, this->getA(), ptr.data(), tg.data(), nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error);;
	if ( error != 0 )
	{
		printf ("\nERROR during numerical factorization: %d", error);
		exit (2);
	}
	printf ("\nFactorization completed ... %d singularities %d negative EV.\n", iparm[13], iparm[23]);
	fflush(stdout);
}

template<typename Scalar>
void PardisoSolver<Scalar>::parallelFactor()
{
	factor();
}

template<typename Scalar>
PardisoSolver<Scalar>::~PardisoSolver()
{

}

template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getPardiso(const Connectivity *cn, const EqNumberer *_dsa) {
	auto solver = new PardisoSolver<Scalar>(cn, _dsa);
	return {solver, solver};
}

#else
template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getPardiso(const Connectivity *cn, const EqNumberer *_dsa) {
	return {nullptr, nullptr};
}

#endif

#include <complex>

template std::pair<GenSolver<double> *, GenSparseMatrix<double> *>
getPardiso<double>(const Connectivity *cn, const EqNumberer *_dsa);

template std::pair<GenSolver<std::complex<double>> *, GenSparseMatrix<std::complex<double>> *>
getPardiso<std::complex<double>>(const Connectivity *cn, const EqNumberer *_dsa);
