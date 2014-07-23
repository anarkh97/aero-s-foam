#include <Solvers.d/Spooles.C>

template<>
void
GenSpoolesSolver<double>::addImaginary(FullSquareMatrix &kel, int *dofs)
{
  fprintf(stderr, "GenSpoolesSolver<double> cannot addImaginary\n");
}

template<>
void
GenSpoolesSolver<complex<double> >
   ::addImaginary(FullSquareMatrix &kel, int *dofs)
{
  int i, j, m, mstart, mstop;
  int kndof = kel.dim();                       // Dimension of element stiff.

  for(i = 0; i < kndof; ++i) {                 // Loop over rows.
    if(unconstrNum[dofs[i]] == -1) continue;   // Skip constrained dofs
    for(j = 0; j < kndof; ++j) {               // Loop over columns.
      if(unconstrNum[dofs[j]] == -1) continue; // Skip constrained dofs
      if(unconstrNum[dofs[j]] < unconstrNum[dofs[i]]) continue;
      mstart = xunonz[unconstrNum[dofs[j]]];
      mstop  = xunonz[unconstrNum[dofs[j]]+1];
      for(m=mstart; m<mstop; ++m) {
        // if(rowu[m-1] > unconstrNum[dofs[j]]+1)
        //   fprintf(stderr, "Bigger: %d %d\n", rowu[m-1]-1, unconstrNum[dofs[j]]);
        if(rowu[m-1] == (unconstrNum[dofs[i]] + 1)) {
          unonz[m-1] += complex<double>(0.0, kel[i][j]);
          break;
        }
      }
    }
  }
}

template<>
void
GenSpoolesSolver<double>::add(FullSquareMatrixC &kel, int *dofs)
{
  fprintf(stderr,"GenSpoolesSolver<double>::add(FullSquareMatrixC &kel, int *dofs) is not implemented.\n");
}

template<>
void
GenSpoolesSolver<complex<double> >::add(FullSquareMatrixC &kel, int *dofs)
{
  int i, j, m, mstart, mstop;
  int kndof = kel.dim();                       // Dimension of element stiff.

  for(i = 0; i < kndof; ++i) {                 // Loop over rows.
    if(unconstrNum[dofs[i]] == -1) continue;   // Skip constrained dofs
    for(j = 0; j < kndof; ++j) {               // Loop over columns.
      if(unconstrNum[dofs[j]] == -1) continue; // Skip constrained dofs
      if(unconstrNum[dofs[j]] < unconstrNum[dofs[i]]) continue; 
      mstart = xunonz[unconstrNum[dofs[j]]];
      mstop  = xunonz[unconstrNum[dofs[j]]+1];
      for(m=mstart; m<mstop; ++m) {
        if(rowu[m-1] == (unconstrNum[dofs[i]] + 1)) {
          unonz[m-1] += kel[i][j];
          break;
        }
      }
    }
  }
}

#define SPOOLES_INSTANTIATION_HELPER(Scalar) \
template \
GenSpoolesSolver<Scalar>::GenSpoolesSolver(Connectivity*, EqNumberer*, int*); \
template \
GenSpoolesSolver<Scalar>::GenSpoolesSolver(Connectivity*, DofSetArray*, ConstrainedDSA*); \
template \
void \
GenSpoolesSolver<Scalar>::add(FullSquareMatrix&, int*); \
template \
void \
GenSpoolesSolver<Scalar>::add(int, int, Scalar); \
template \
void \
GenSpoolesSolver<Scalar>::add(GenFullM<Scalar>&, int*); \
template \
void \
GenSpoolesSolver<Scalar>::add(GenAssembledFullM<Scalar>&, int*); \
template \
void \
GenSpoolesSolver<Scalar>::add(GenFullM<Scalar>&, int, int); \
template \
void \
GenSpoolesSolver<Scalar>::addDiscreteMass(int, Scalar); \
template \
Scalar \
GenSpoolesSolver<Scalar>::diag(int) const; \
template \
Scalar & \
GenSpoolesSolver<Scalar>::diag(int); \
template \
void \
GenSpoolesSolver<Scalar>::symmetricScaling(); \
template \
void \
GenSpoolesSolver<Scalar>::applyScaling(Scalar*); \
template \
void \
GenSpoolesSolver<Scalar>::factor(); \
template \
void \
GenSpoolesSolver<Scalar>::parallelFactor(); \
template \
void \
GenSpoolesSolver<Scalar>::allFactor(bool); \
template \
void \
GenSpoolesSolver<Scalar>::solve(Scalar*); \
template \
void \
GenSpoolesSolver<Scalar>::solve(Scalar*, Scalar*); \
template \
double \
GenSpoolesSolver<Scalar>::getMemoryUsed(void); \
template \
void \
GenSpoolesSolver<Scalar>::print(); \
template \
long \
GenSpoolesSolver<Scalar>::size(); \
template \
void \
GenSpoolesSolver<Scalar>::unify(FSCommunicator*); \
template \
void \
GenSpoolesSolver<Scalar>::zeroAll();  \
template \
void \
GenSpoolesSolver<Scalar>::cleanUp(); \
template \
GenSpoolesSolver<Scalar>::~GenSpoolesSolver(); \
template \
void \
GenSpoolesSolver<Scalar>::init(); \

SPOOLES_INSTANTIATION_HELPER(double);
SPOOLES_INSTANTIATION_HELPER(complex<double>);

