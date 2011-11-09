// $Id$

#include "ContactHostGlobalID.h"
#ifndef CONTACT_NO_MPI
#include "ContactZoltanID.h"
#endif
#include "Contact_Defines.h"
#include "ContactParOStream.h"
#include <iostream>

//ContactHostGlobalID::~ContactHostGlobalID(){}

std::ostream& operator<<( std::ostream& os, const ContactHostGlobalID& gid )
{
  os << "(" << gid.hi_int << "," << gid.lo_int << ")";
  return os;
}

ContactParOStream& operator<<( ContactParOStream& pos, 
			       const ContactHostGlobalID& gid )
{
  pos << "(" << gid.hi_int << "," << gid.lo_int << ")";
  return pos;
}
