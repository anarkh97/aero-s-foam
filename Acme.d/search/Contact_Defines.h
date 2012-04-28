// $Id$

#ifndef Contact_Defines_h_
#define Contact_Defines_h_

typedef double Real;
typedef int VariableHandle;
typedef int ScratchVarHandle;

#ifndef NULL
#define NULL 0
#endif

#ifndef FORTRAN

#if defined(LCFLINK)
# define FORTRAN(subname) subname
# define FORTRANNO_(subname) subname
#elif defined(LC_FLINK)
# define FORTRAN(subname) subname##_
# define FORTRANNO_(subname) subname##_
#elif defined(LC__FLINK)
# define FORTRAN(subname) subname##__
# define FORTRANNO_(subname) subname##_
#else
# define FORTRAN(subname) subname##_
# define FORTRANNO_(subname) subname##_
#endif

#endif

#ifndef CONTACT_NO_MPI
#define ZOLTAN_LID_SIZE 2
#define ZOLTAN_GID_SIZE 3
#else
typedef int MPI_Comm;
#endif
#if !defined(__PUMAGON__) && !defined(__sgi) && !defined(AIX) && !defined(OSF1)
#define MPI_INT_TYPE MPI_Fint
#define MPI_COMM_F2C(icomm) MPI_Comm_f2c(icomm)
#define MPI_COMM_C2F(comm) MPI_Comm_c2f(comm)
#else
#define MPI_INT_TYPE int
#define MPI_COMM_F2C(icomm) (MPI_Comm)icomm
#define MPI_COMM_C2F(comm) (MPI_INT_TYPE)comm
#endif

//---------------------------------------------------------------
//  
//                        N U M B E R S
//  
//---------------------------------------------------------------

#define BIGNUM 1.E308
#define TINYNUM 1.E-200
#define PI 3.14159265358979
#define MAX_DIMENSIONALITY 3
#define MAX_NODES_PER_FACE 9
#define MAX_NODE_ENTITY_INTERACTIONS_PER_NODE 3

//
//  Define the maximum depth of the search tree for various algorithms.  
//  The maximum allowable number of elements will be approximatly 
//  2**MAX_TREE_LEVELS
//
#define MAX_TREE_LEVELS 100

//
//  Facet array sizes for use in NODE-FACE searches
//
#define MAX_FACETS 32

//
//  Use either 1- or 2-step ghosting method
//
#define CONTACT_2_STEP_GHOSTING

//
// hash bin size paramater
//
#define BIN_FRACTION 0.5

//
// Select hash function
//
//#define JENKINS96_HASH_FUNC
//#define JENKINS32_HASH_FUNC
//#define WANG32_HASH_FUNC
//#define SEDGEWICK_HASH_FUNC
#define KNUTH_HASH_FUNC

//
// Select type of hash binning
//
#if defined (JENKINS96_HASH_FUNC) || defined (JENKINS32_HASH_FUNC) || defined (WANG32_HASH_FUNC)
#  define POWER2_HASH_BINS
#  define BIN_MINIMUM 128
#else
//#  define PRIME_HASH_BINS
#  define PRIME_LIKE_HASH_BINS
#  define BIN_MINIMUM 127
#endif

#define NO_CONTACT_ANALYZE_HASH

#define NO_CONTACT_ANALYZE_DATA_XFER

#define NO_CONTACT_HEARTBEAT

#define NO_CONTACT_OLD_XFER

#if defined (__osf__) && !defined(__PUMAGON__) && !defined(CONTACT_USE_BLOCKING_SEND)
#define CONTACT_USE_BLOCKING_SEND
#endif

//
//  Number of derivatives for Face-Face interactions
//  use 54 for dimensionality 3 and support for linear and quadratic faces
//  use 24 for dimensionality 3 and support for linear faces only
//  use 16 for dimensionality 2 and support for linear faces only
//
#define MAX_FFI_DERIVATIVES 24
//#define COMPUTE_FFI_SECOND_DERIVATIVES

#define FFI_AUTO_DIFF_SACADO

#if (MAX_FFI_DERIVATIVES > 0) && defined(USE_SACADO) && defined(FFI_AUTO_DIFF_SACADO)
#  include "Sacado.hpp"
   typedef Sacado::Fad::SFad<Real,MAX_FFI_DERIVATIVES> ActiveScalar;
   inline ActiveScalar InitActiveScalar(int nbDer, int derNumber, const Real& value) { return ActiveScalar(nbDer, derNumber, value); }
   inline Real GetActiveScalarValue(const ActiveScalar& s) { return s.val(); }
   inline Real GetActiveScalarDerivative(const ActiveScalar& s, int i) { return s.dx(i); }
   const int MAX_FFI_SECOND_DERIVATIVES = 0;
#elif (MAX_FFI_DERIVATIVES > 0) && defined(USE_EIGEN3) && !defined(COMPUTE_FFI_SECOND_DERIVATIVES)
#  include <Eigen/Core>
#  include <unsupported/Eigen/AutoDiff>
   typedef Eigen::Matrix<Real, MAX_FFI_DERIVATIVES, 1> DerivativeType;
   typedef Eigen::AutoDiffScalar<DerivativeType> ActiveScalar;
   inline ActiveScalar InitActiveScalar(int nbDer, int derNumber, const Real& value) { return ActiveScalar(value, nbDer, derNumber); }
   inline Real GetActiveScalarValue(const ActiveScalar& s) { return s.value(); }
   inline Real GetActiveScalarDerivative(const ActiveScalar& s, int i) { return s.derivatives()[i]; }
   const int MAX_FFI_SECOND_DERIVATIVES = 0;
#elif (MAX_FFI_DERIVATIVES > 0) && defined(USE_EIGEN3) && defined(COMPUTE_FFI_SECOND_DERIVATIVES)
#  include <Eigen/Core>
#  include <unsupported/Eigen/AutoDiff>
   typedef Eigen::Matrix<Real, MAX_FFI_DERIVATIVES, 1> DerivativeType1;
   typedef Eigen::Matrix<Eigen::AutoDiffScalar<DerivativeType1>, MAX_FFI_DERIVATIVES, 1> DerivativeType2;
   typedef Eigen::AutoDiffScalar<DerivativeType2> ActiveScalar;
   inline ActiveScalar InitActiveScalar(int nbDer, int derNumber, const Real& value) { 
     return ActiveScalar(Eigen::AutoDiffScalar<DerivativeType1>(value,nbDer,derNumber), nbDer, derNumber); 
   }
   inline Real GetActiveScalarValue(const ActiveScalar& s) { return s.value().value(); }
   inline Real GetActiveScalarDerivative(const ActiveScalar& s, int i) { return s.derivatives()[i].value(); }
   const int MAX_FFI_SECOND_DERIVATIVES = MAX_FFI_DERIVATIVES*(MAX_FFI_DERIVATIVES+1)/2;
   inline Real GetActiveScalarSecondDerivative(const ActiveScalar& s, int i, int j) { return s.derivatives()[i].derivatives()[j]; }
#else
#  error "Computation of Face-Face interaction derivatives requires Sacado or Eigen3"
#endif

#endif // #ifdef Contact_Defines_h_
