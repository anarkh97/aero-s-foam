#ifndef ROM_VECBASISUTILS_H
#define ROM_VECBASISUTILS_H

#include "VecBasis.h"

// Computes result^{T} := targetPod^{T} * originPod * originProjection^{T} 
template <typename Scalar>
const GenVecBasis<Scalar> &
combine_projections(const GenVecBasis<Scalar> &targetPod,
                    const GenVecBasis<Scalar> &originPod,
                    const GenVecBasis<Scalar> &originProjection,
                    GenVecBasis<Scalar> &result);

#endif /* ROM_VECBASISUTILS_H */
