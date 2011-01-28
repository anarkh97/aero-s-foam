#ifndef ROM_VECBASISFILE_H
#define ROM_VECBASISFILE_H

#include "VecBasis.h"

class BasisInputStream;
class BasisOutputStream;

BasisInputStream &
operator>>(BasisInputStream &, VecBasis &);

BasisOutputStream &
operator<<(BasisOutputStream &, const VecBasis &);

#endif /* ROM_VECBASISFILE_H */
