//
// Created by Michel Lesoinne on 2018-12-14.
//

#include <complex>
#include <memory>

#include "CholmodSolver.h"
#include "CholmodImp.h"
#include <Math.d/DBSparseMatrix.h>


#ifndef WITH_CHOLMOD
template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getCholmod(const Connectivity *cn, const EqNumberer *_dsa)
{
	std::cerr << "Cholmod is not available in this version of the code. Crash will ensue." << std::endl;
	return {nullptr, nullptr};
}
#else // WITH_CHOLMOD

class CholmodImp;

template <typename Scalar>
class CholmodSolver :
	public GenDBSparseMatrix<Scalar>, public GenSolver<Scalar> {
public:
	CholmodSolver(const Connectivity *cn, const EqNumberer *_dsa);

	~CholmodSolver();

	int neqs() const override;

	void solve(const Scalar *rhs, Scalar *solution) override;

	long size() const;

	void factor() override;

	void parallelFactor() override;
private:
	CholmodImp &getImplementation();

	std::unique_ptr<CholmodImp> impl;
};
template <typename Scalar>
CholmodSolver<Scalar>::CholmodSolver(const Connectivity *cn, const EqNumberer *_dsa)  :
	GenDBSparseMatrix<Scalar>(cn, _dsa)
{

}

template <typename Scalar>
CholmodSolver<Scalar>::~CholmodSolver()
{

}

template<typename Scalar>
int CholmodSolver<Scalar>::neqs() const
{
	return GenDBSparseMatrix<Scalar>::neqs();
}

template<typename Scalar>
long CholmodSolver<Scalar>::size() const
{
	return GenDBSparseMatrix<Scalar>::size()
	       + (impl ? impl->memorySize() : 0);
}

template<typename Scalar>
void CholmodSolver<Scalar>::solve(const Scalar *rhs, Scalar *solution)
{
//	TODO
}

template<typename Scalar>
void CholmodSolver<Scalar>::factor()
{
	auto &cholmodImp = getImplementation();
	std::cout << "Setting the data" << std::endl;
	cholmodImp.setData(*this);
	std::cout << "Factorizing." << std::endl;
	cholmodImp.factorize();
	std::cout << "Done." << std::endl;
	throw "Blob";
}


template<typename Scalar>
void CholmodSolver<Scalar>::parallelFactor()
{
	factor(); // TODO make parallel.
}

template<typename Scalar>
CholmodImp &CholmodSolver<Scalar>::getImplementation()
{
	if ( !impl )
		impl = std::make_unique<CholmodImp>(*this, std::is_same<std::complex<double>, Scalar>::value);
	return *impl;
}

template class CholmodSolver<double>;
template class CholmodSolver<std::complex<double>>;

template <typename Scalar>
std::pair<GenSolver<Scalar> *, GenSparseMatrix<Scalar> *>
getCholmod(const Connectivity *cn, const EqNumberer *dsa)
{
	auto solver = new CholmodSolver<Scalar>(cn, dsa);
	return { solver, solver};
}

#endif // WITH_CHOLMOD

template std::pair<GenSolver<double> *, GenSparseMatrix<double> *>
getCholmod<double>(const Connectivity *cn, const EqNumberer *_dsa);
template std::pair<GenSolver<std::complex<double>> *, GenSparseMatrix<std::complex<double>> *>
getCholmod<std::complex<double>>(const Connectivity *cn, const EqNumberer *_dsa);