#include "BasisOrthoDriver.h"

#include "SvdOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "BasisFileStream.h"
#include "FileNameInfo.h"
#include "SimpleBuffer.h"

#include <Driver.d/Domain.h>
#include <Utils.d/dofset.h>
#include <Utils.d/DistHelper.h>

#include <utility>
#include <algorithm>

namespace Rom {

namespace { // anonymous

template <typename VecType>
class VectorTransform {
public:
  virtual void operator()(VecType &) const = 0;

  virtual ~VectorTransform() {}
};

template <typename VecType>
class NoOp : public VectorTransform<VecType> {
public:
  virtual void operator()(VecType &) const  { /* Nothing */ }
};

template <typename VecType>
class RefSubstraction : public VectorTransform<VecType> {
public:
  virtual void operator()(VecType &v) const;

  explicit RefSubstraction(const Domain *);
private:
  Vector ref_;
};

template <typename VecType>
void
RefSubstraction<VecType>::operator()(VecType &v) const {
  for (int i = 0; i < ref_.size(); ++i) {
    v[i] -= ref_[i];
  }
}

template <typename VecType>
RefSubstraction<VecType>::RefSubstraction(const Domain *domain) :
  ref_(const_cast<Domain *>(domain)->numUncon())
{
  Vector dummy(const_cast<Domain *>(domain)->numUncon());
  const_cast<Domain *>(domain)->initDispVeloc(ref_, dummy, dummy, dummy);
}

} // end anonymous namespace

BasisOrthoDriver::BasisOrthoDriver(Domain *domain) :
  SingleDomainDynamic(domain)
{}

void
BasisOrthoDriver::solve() {
  //preProcess() ;
  SingleDomainDynamic::preProcess();
  VecNodeDof6Conversion converter(*domain->getCDSA()) ;
  FileNameInfo fileInfo ;
  SvdOrthogonalization solver;

  std::vector<BasisId::Type> workload ;
       
  if(domain->solInfo().statevectPodRom) {
	workload.push_back(BasisId::STATE);
        fprintf(stderr," ... For State SVD, workload size = %zd ...\n", workload.size());}
  else if(domain->solInfo().residvectPodRom) {
	workload.push_back(BasisId::RESIDUAL);
	fprintf(stderr," ... For Residual SVD, workload size = %zd ...\n", workload.size());}
  else if(domain->solInfo().jacobvectPodRom) {
        workload.push_back(BasisId::JACOBIAN);
	fprintf(stderr," ... For Jacobian SVD, workload size = %zd ...\n", workload.size());}
  else if(domain->solInfo().forcevectPodRom) {
        workload.push_back(BasisId::FORCE);
	fprintf(stderr," ... For Force SVD, workload size = %zd ...\n", workload.size());}
  else if(domain->solInfo().accelvectPodRom) {
        workload.push_back(BasisId::ACCELERATION);
	fprintf(stderr," ... For Acceleration SVD, workload size = %zd ...\n", workload.size());}
  else if(domain->solInfo().velocvectPodRom) {
        workload.push_back(BasisId::VELOCITY);
        fprintf(stderr," ... For Velocity SVD, workload size = %zd ...\n", workload.size());}
  else { workload.push_back(BasisId::STATE);
	fprintf(stderr," ... For default SVD, workload size = %zd ...\n", workload.size());}

  typedef VectorTransform<double *> VecTrans;
  std::auto_ptr<VecTrans> transform(domain->solInfo().substractRefPodRom ?
                                    static_cast<VecTrans *>(new RefSubstraction<double *>(domain)) :
                                    static_cast<VecTrans *>(new NoOp<double *>));
  //Checking flags
  double beta = domain->solInfo().newmarkBeta;
  //Assembling mass matrix
  DynamMat * dummyDynOps = SingleDomainDynamic::buildOps(1.0,0.0,0.0) ;
  assert(dummyDynOps->M);
  GenSparseMatrix<double> *fullMass = dummyDynOps->M;

  for (std::vector<BasisId::Type>::const_iterator it = workload.begin(); it != workload.end(); ++it) {
    BasisId::Type type = *it;

    {
      BasisInputStream input(BasisFileId(fileInfo, type, BasisId::SNAPSHOTS), converter);
      filePrint(stderr, " ... Computation of a basis of size %d...\n", input.size());
      
      solver.matrixSizeIs(input.vectorSize(), input.size());

      int iCol = 0;
      for (int iCol = 0; iCol < solver.colCount(); ++iCol) {
        double *buffer = solver.matrixCol(iCol);
        input >> buffer;
        assert(input);
        if(domain->solInfo().normalize ==1) fullMass->squareRootMult(buffer); //executes for new method
        (*transform)(buffer);
      }
    }
    solver.solve();

    BasisOutputStream output(BasisFileId(fileInfo, type, BasisId::POD), converter, false); 

    const int orthoBasisDim = domain->solInfo().maxSizePodRom ?
                              std::min(domain->solInfo().maxSizePodRom, solver.singularValueCount()) :
                              solver.singularValueCount();

    //Output solution
    for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
      output << std::make_pair(solver.singularValue(iVec), solver.matrixCol(iVec));
    }

    //Check if explicit
    if(beta == 0 ) {
      //Read back in solution to renormalize basis
      VecBasis basis;
      //Read back in basis information
      BasisInputStream in(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD), converter);
      readVectors(in, basis);

      VecBasis normalizedBasis;
      //Check renormalize tag
      if(domain->solInfo().normalize == 0){
        renormalized_basis(*fullMass, basis, normalizedBasis) ;
      }
      if(domain->solInfo().normalize == 1){
        for(int col = 0 ; col < orthoBasisDim; col ++ ){
          fullMass->inverseSquareRootMult(basis[col].data());
        }
        normalizedBasis = basis;
      }
    
      std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
      fileName.append(".normalized");
      BasisOutputStream outputNormalized(fileName, converter, false); 
      for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
        outputNormalized << std::make_pair(solver.singularValue(iVec), normalizedBasis[iVec]);
        //outputNormalized << std::make_pair(normalizedBasis[iVec].data(), solver.singularValue(iVec));
      }
    
    }

  }
}

void
BasisOrthoDriver::preProcess() {
  domain->preProcessing();
 
  // Build the constrained DofSetArray incorporating the boundary conditions 
  const int numdof = domain->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain->make_bc(bc.array(), bcx.array());
  domain->make_constrainedDSA(bc.array());
}

} /* end namespace Rom */

Rom::DriverInterface *basisOrthoDriverNew(Domain *domain) {
  return new Rom::BasisOrthoDriver(domain);
}
