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
  domain_(domain)
{}

void
BasisOrthoDriver::solve() {
  preProcess();
 
  VecNodeDof6Conversion converter(*domain_->getCDSA());
  FileNameInfo fileInfo;
  SvdOrthogonalization solver;

  std::vector<BasisId::Type> workload;

       if(domain_->solInfo().statevectPodRom) {
	workload.push_back(BasisId::STATE);
        fprintf(stderr,"For State SVD, workload size = %d \n", workload.size());}
  else if(domain_->solInfo().residvectPodRom) {
	workload.push_back(BasisId::RESIDUAL);
	fprintf(stderr,"For Residual SVD, workload size = %d \n", workload.size());}
  else if(domain_->solInfo().jacobvectPodRom) {
        workload.push_back(BasisId::JACOBIAN);
	fprintf(stderr,"For Jacobian SVD, workload size = %d \n", workload.size());}
  else if(domain_->solInfo().forcevectPodRom) {
        workload.push_back(BasisId::FORCE);
	fprintf(stderr,"For Force SVD, workload size = %d \n", workload.size());}
  else if(domain_->solInfo().accelvectPodRom) {
        workload.push_back(BasisId::ACCELERATION);
	fprintf(stderr,"For Acceleration SVD, workload size = %d \n", workload.size());}
  else { workload.push_back(BasisId::STATE);
	fprintf(stderr,"For default SVD, workload size = %d \n", workload.size());}

  typedef VectorTransform<double *> VecTrans;
  std::auto_ptr<VecTrans> transform(domain_->solInfo().substractRefPodRom ?
                                    static_cast<VecTrans *>(new RefSubstraction<double *>(domain_)) :
                                    static_cast<VecTrans *>(new NoOp<double *>));

  for (std::vector<BasisId::Type>::const_iterator it = workload.begin(); it != workload.end(); ++it) {
    BasisId::Type type = *it;

    {
      BasisInputStream input(BasisFileId(fileInfo, type, BasisId::SNAPSHOTS), converter);
      filePrint(stderr, "Orthogonalization of a basis with %d vectors\n", input.size());
      
      solver.matrixSizeIs(input.vectorSize(), input.size());

      int iCol = 0;
      for (int iCol = 0; iCol < solver.colCount(); ++iCol) {
        double *buffer = solver.matrixCol(iCol);
        input >> buffer;
        assert(input);
        (*transform)(buffer);
      }
    }

    solver.solve();

    BasisOutputStream output(BasisFileId(fileInfo, type, BasisId::POD), converter);
    const int orthoBasisDim = domain->solInfo().maxSizePodRom ?
                              std::min(domain->solInfo().maxSizePodRom, solver.singularValueCount()) :
                              solver.singularValueCount();

    for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
      output << std::make_pair(solver.singularValue(iVec), solver.matrixCol(iVec));
    }
  }
}

void
BasisOrthoDriver::preProcess() {
  domain_->preProcessing();
 
  // Build the constrained DofSetArray incorporating the boundary conditions 
  const int numdof = domain_->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain_->make_bc(bc.array(), bcx.array());
  domain_->make_constrainedDSA(bc.array());
}

} /* end namespace Rom */

Rom::DriverInterface *basisOrthoDriverNew(Domain *domain) {
  return new Rom::BasisOrthoDriver(domain);
}
