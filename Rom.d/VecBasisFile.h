#ifndef ROM_VECBASISFILE_H
#define ROM_VECBASISFILE_H

#include "VecBasis.h"

class BasisInputStream;
class BasisOutputStream;

// Reset basis with the full content of the stream
BasisInputStream &
operator>>(BasisInputStream &, VecBasis &);

// Reset basis with a partial content of the stream
BasisInputStream &
readVectors(BasisInputStream &, VecBasis &, int countMax);

// Output the full content of the basis
BasisOutputStream &
operator<<(BasisOutputStream &, const VecBasis &);

// Output a partial content of the basis
BasisOutputStream &
writeVectors(BasisOutputStream &, const VecBasis &, int countMax);

#endif /* ROM_VECBASISFILE_H */
