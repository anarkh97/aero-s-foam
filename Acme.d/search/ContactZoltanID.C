// $Id$

#ifndef CONTACT_NO_MPI

#include "ContactZoltanID.h"
#include "ContactHostGlobalID.h"
#include "Contact_Defines.h"

// Classes to implement Zoltan local and global IDs.  
//
// An array of two unsigned ints is used to describe the local ID. 
// 
//    lid[0] = entity_type
//    lid[1] = index on processer
//
//An array of three unsigned ints is used to describe the global ID.
// 
//    gid[0] = entity_type
//    gid[1] = host_id[0]
//    gid[1] = host_id[1]


ContactZoltanLID::ContactZoltanLID() : 
  type(-1), index(-1)
{}

ContactZoltanLID::~ContactZoltanLID(){}

//=======================================================================

ContactZoltanGID::ContactZoltanGID() : 
  type(-1), hi(-1), lo(-1)
{}

ContactZoltanGID::~ContactZoltanGID(){}

void 
ContactZoltanGID::ZoltanGID(int id_type, ContactHostGlobalID* gid, 
                            LB_ID_PTR id)
{
  id[0] = id_type;
  id[1] = gid->HiInt();
  id[2] = gid->LoInt();
}

#endif
