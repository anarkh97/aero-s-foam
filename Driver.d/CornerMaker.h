#ifndef _CORNER_MAKER_H_
#define _CORNER_MAKER_H_

#include <vector>
#include <Feti.d/CornerSelector.h>

class Elemset;

class SubCornerHandler : public FetiSubCornerHandler {
public:
	SubCornerHandler(int sub, int nn, CoordSet &n, Elemset &ele, Connectivity &nTn, DofSetArray &d,
		                 Connectivity &sh, int *nsb, ConstrainedDSA *c_dsa, FetiBaseSub *_subPre);
};

#endif
