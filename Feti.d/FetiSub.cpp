//
// Created by Michel Lesoinne on 12/6/17.
//
#include <Utils.d/dofset.h>
#include "FetiSub.h"

void
FetiBaseSub::markCornerDofs(int *glCornerDofs) const
{
	for(int i=0; i<numCRN; ++i)
		glCornerDofs[glCornerNodes[i]] |= cornerDofs[i].list();
}