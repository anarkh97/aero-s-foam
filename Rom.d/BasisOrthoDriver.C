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
  SingleDomainDynamic::preProcess();
  VecNodeDof6Conversion converter(*domain->getCDSA());
  FileNameInfo fileInfo;
  SvdOrthogonalization solver;

  std::vector<BasisId::Type> workload;
       
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

  double beta = domain->solInfo().newmarkBeta;
  // Assembling mass matrix
  DynamMat * dummyDynOps = SingleDomainDynamic::buildOps(1.0,0.0,0.0);
  assert(dummyDynOps->M);
  GenSparseMatrix<double> *fullMass = dummyDynOps->M;

  int vectorSize = 0; // size of vectors
  int sizeSnap = 0; // number of state snapshots
  int sizeROB = 0;
  int skipTime = domain->solInfo().skipPodRom;
  int skip;
  if(domain->solInfo().snapfiPodRom.empty() && domain->solInfo().robfi.empty()) {
    std::cerr << "*** Error: no files provided\n";
    exit(-1);
  }
 
  for (std::vector<BasisId::Type>::const_iterator it = workload.begin(); it != workload.end(); ++it) {
    BasisId::Type type = *it;
    // Loop over snapshots
    for(int i = 0; i < domain->solInfo().snapfiPodRom.size(); i++) {
      std::string fileName =  BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
      BasisInputStream input(fileName, converter);
      vectorSize = input.vectorSize();
      sizeSnap += input.size()/skipTime;
    }

    // Loop over rob files 
    for(int i = 0; i < domain->solInfo().robfi.size(); i++) {
      std::string fileName = BasisFileId(fileInfo,type,BasisId::ROB, i);
      BasisInputStream input(fileName, converter);
      vectorSize = input.vectorSize();
      sizeROB += input.size();
    }
  }
  solver.matrixSizeIs(vectorSize, sizeSnap+sizeROB);

  for (std::vector<BasisId::Type>::const_iterator it = workload.begin(); it != workload.end(); ++it) {
    BasisId::Type type = *it;
    filePrint(stderr, " ... Computation of a basis of size %d ...\n", sizeSnap+sizeROB);
    {
      int colCounter = 0 ; // Column counter for combined matrix
      for(int i = 0 ; i < domain->solInfo().snapfiPodRom.size(); i++) {
        std::string fileName =  BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
        BasisInputStream input(fileName, converter);
        filePrint(stderr, " ... Reading in snapshot file: %s ...\n", fileName.c_str());
        skip = 1;
        for (int iCol = 0; iCol < input.size(); ++iCol) {
          if(skip == skipTime) {
            double *buffer = solver.matrixCol(colCounter);
            input >> buffer;
            assert(input);
            colCounter++;
            // Multiply by weighting factor if given in input file
            if(!domain->solInfo().snapshotWeights.empty()) {
              for(int row = 0 ; row < vectorSize; row++) {
                buffer[row] *= domain->solInfo().snapshotWeights[i];
              }
            }
            if(beta == 0 && domain->solInfo().normalize == 1) fullMass->squareRootMult(buffer); // new method
            (*transform)(buffer);
            skip = 1;
          } else {
            SimpleBuffer<double> dummyVec;
            dummyVec.sizeIs(input.vectorSize());  
            double *dummyBuffer = dummyVec.array();
            input >> dummyBuffer;
            assert(input);
            ++skip;
          }
        }
      }
      
      for(int i = 0; i < domain->solInfo().robfi.size(); i++) {
        std::string fileName = BasisFileId(fileInfo, type, BasisId::ROB, i);
        BasisInputStream input(fileName,converter);
        filePrint(stderr, " ... Reading in ROB file: %s ...\n", fileName.c_str());
        for(int iCol = 0; iCol < input.size(); ++iCol) {
            double *buffer = solver.matrixCol(colCounter);
            std::pair<double, double *> data;
            assert(input);
            data.second = buffer;
            input >> data;
            colCounter++;
            for(int row = 0 ; row < vectorSize; row++) {
              data.second[row] *= data.first; //multiply in singluar values
              if(!domain->solInfo().snapshotWeights.empty()){
                data.second[row] *= domain->solInfo().snapshotWeights[domain->solInfo().snapfiPodRom.size()+i]; //multiply in weights if given
              }
            }
            if(beta == 0 && domain->solInfo().normalize == 1) fullMass->squareRootMult(buffer); // new method
            (*transform)(buffer);
        }
      }

    }
    solver.solve();

    BasisOutputStream output(BasisFileId(fileInfo, type, BasisId::POD), converter, false); 

    const int orthoBasisDim = domain->solInfo().maxSizePodRom ?
                              std::min(domain->solInfo().maxSizePodRom, solver.singularValueCount()) :
                              solver.singularValueCount();

    // Output solution
    if(beta != 0 || (beta == 0 && domain->solInfo().normalize == 0))
      filePrint(stderr, " ... Writing orthonormal basis to file %s ...\n", BasisFileId(fileInfo, type, BasisId::POD).name().c_str());
    for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
      output << std::make_pair(solver.singularValue(iVec), solver.matrixCol(iVec));
    }

    // Check if explicit
    if(beta == 0) {
      // Read back in output file to renormalize basis
      VecBasis basis;
      BasisInputStream in(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD), converter);
      readVectors(in, basis);

      VecBasis normalizedBasis;
      if(domain->solInfo().normalize == 0) {
        // Old method: renormalize the orthonormal basis
        renormalized_basis(*fullMass, basis, normalizedBasis);
      }
      else if(domain->solInfo().normalize == 1) {
        // New method: multiply by inverse square root of the mass matrix
        for(int col = 0; col < orthoBasisDim; col ++ ) {
          fullMass->inverseSquareRootMult(basis[col].data());
        }
        normalizedBasis = basis;
      }
    
      // Output the renormalized basis as separate file
      std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
      fileName.append(".normalized");
      BasisOutputStream outputNormalized(fileName, converter, false); 
      filePrint(stderr, " ... Writing mass-normalized basis to file %s ...\n", fileName.c_str());
      for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
        outputNormalized << std::make_pair(solver.singularValue(iVec), normalizedBasis[iVec]);
      }
    
      // Compute and output orthonormal basis if using new method
      if(domain->solInfo().normalize == 1) {
        std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
        MGSVectors(normalizedBasis.data(), normalizedBasis.numVec(), normalizedBasis.size());
        BasisOutputStream outputIdentityNormalized(fileName, converter, false); 
        filePrint(stderr, " ... Writing orthonormal basis to file %s ...\n", fileName.c_str());
        for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
          outputIdentityNormalized << std::make_pair(solver.singularValue(iVec), normalizedBasis[iVec]);
        }
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
