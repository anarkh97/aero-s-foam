// $Id$

#include "ContactInteractionEntity.h"
#include "ContactParOStream.h"
#include "contact_assert.h"
#include <cstddef>
#include <cstring>

ContactInteractionEntity::ContactInteractionEntity(Real *data_array_, ContactType base_type_)
  : data_array(data_array_), base_type(base_type_), index(-1), proc_index(-1)
{}

ContactInteractionEntity::~ContactInteractionEntity()
{}

void 
ContactInteractionEntity::Display(ContactParOStream& postream)
{
  postream<<"ContactEntity: \n";
  postream<<"               entity type:       ";
  switch (Base_Type()) {
  case CT_NNI:
    postream<<"CT_NNI\n";
    break;
  case CT_NFI:
    postream<<"CT_NFI\n";
    break;
  case CT_NSI:
    postream<<"CT_NSI\n";
    break;
  case CT_FFI:
    postream<<"CT_FFI\n";
    break;
  case CT_FCI:
    postream<<"CT_FCI\n";
    break;
  case CT_EEI:
    postream<<"CT_EEI\n";
    break;
  case CT_UNKNOWN:
    postream<<"CT_UNKNOWN\n";
    break;
  default:
    POSTCONDITION(0);
    break;
  }
  postream<<"               index:             "<<index<<"\n";
  postream<<"               proc_index:        "<<proc_index<<"\n";
}
