// $Id$

#ifndef ContactEntity_H_
#define ContactEntity_H_

#include "Contact_Defines.h"
#include "contact_assert.h"
#include <cstring>

  enum ContactType { CT_UNKNOWN, CT_NODE, CT_SHELL_NODE,
		     CT_EDGE, CT_FACE, CT_ELEM, 
		     CT_ELEMENT, CT_ANALYTIC_SURFACE, 
                     CT_NFI, CT_NSI, CT_NNI,
                     CT_FFI, CT_FCI, CT_EEI,
		     CT_NUM_ENTITY_TYPES };
  
#endif  // ifdef ContactEntity_H_
