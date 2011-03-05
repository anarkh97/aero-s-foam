#ifndef ROM_VECBASISFILE_H
#define ROM_VECBASISFILE_H

#include "VecBasis.h"

namespace Rom {

class BasisInputStream;
class BasisOutputStream;

class NodalRestrictionMapping;

// Output
// Output the full content of the basis
BasisOutputStream &
operator<<(BasisOutputStream &, const VecBasis &);

// Output a partial content of the basis
BasisOutputStream &
writeVectors(BasisOutputStream &, const VecBasis &, int countMax);

// Output the full extended content of the basis
BasisOutputStream &
writeExtendedVectors(BasisOutputStream &, const VecBasis &, const NodalRestrictionMapping &);

// Output the full extended content of the basis
BasisOutputStream &
writeRestrictedVectors(BasisOutputStream &, const VecBasis &, const NodalRestrictionMapping &);

// Input
// Reset basis with the full content of the stream (inverse of operator<<)
BasisInputStream &
operator>>(BasisInputStream &, VecBasis &);

// Reset basis with a partial content of the stream (inverse of writeVectors)
BasisInputStream &
readVectors(BasisInputStream &, VecBasis &, int countMax);

// Reset basis with the extended content of the stream (inverse of writeRestrictedVectors)
BasisInputStream &
readExtendedVectors(BasisInputStream &, VecBasis &, const NodalRestrictionMapping &);

// Reset basis with the restricted content of the stream (inverse of writeExtendedVectors)
BasisInputStream &
readRestrictedVectors(BasisInputStream &, VecBasis &, const NodalRestrictionMapping &);

} /* end namespace Rom */

#endif /* ROM_VECBASISFILE_H */
