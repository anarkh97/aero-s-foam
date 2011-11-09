// $Id: AcmeZoltanID.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef CONTACT_NO_MPI

#include "AcmeZoltanID.h"
#include "AcmeEntity.h"
#include "Contact_Defines.h"
#include <iostream.h>

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
//    gid[1] = host_id


AcmeZoltanLID::AcmeZoltanLID(){}

AcmeZoltanLID::~AcmeZoltanLID(){}

AcmeZoltanLID::AcmeZoltanLID(LB_ID_PTR id)
{
  type  = id[0];
  index = id[1];
}

void AcmeZoltanLID::ZoltanLID(int id_type, int indx, LB_ID_PTR id)
{
  id[0] = id_type;
  id[1] = indx;
}

bool AcmeZoltanLID::operator>( const AcmeZoltanLID& id ) const
{
  if( type > id.type )
    return true;
  else if( type < id.type )
    return false;
    
  if( index > id.index )
    return true;
  
  return false;
}


bool AcmeZoltanLID::operator<( const AcmeZoltanLID& id ) const
{
  if( type < id.type )
    return true;
  else if( type > id.type )
    return false;
    
  if( index < id.index )
    return true;
    
  return false;
}


bool AcmeZoltanLID::operator==( const AcmeZoltanLID& id ) const
{
  if( (id.type == type) && (id.index == index) ) return true;
  return false;
}

int AcmeZoltanLID::operator%( const int& divisor ) const
{
  return ( index % divisor );
}

//=======================================================================

AcmeZoltanGID::AcmeZoltanGID(){}

AcmeZoltanGID::~AcmeZoltanGID(){}

AcmeZoltanGID::AcmeZoltanGID(LB_ID_PTR id)
{
  type = id[0];
  gid  = id[1];
}

void AcmeZoltanGID::ZoltanGID(int id_type, int g_id,  LB_ID_PTR id)
{
  id[0] = id_type;
  id[1] = g_id;
}

bool AcmeZoltanGID::operator>( const AcmeZoltanGID& id ) const
{
  if( type > id.type )
    return true;
  else if( type < id.type )
    return false;
    
  if( gid > id.gid )
    return true;
    
  return false;
}


bool AcmeZoltanGID::operator<( const AcmeZoltanGID& id ) const
{
  if( type < id.type )
    return true;
  else if( type > id.type )
    return false;
    
  if( gid < id.gid )
    return true;
    
  return false;
}


bool AcmeZoltanGID::operator==( const AcmeZoltanGID& id ) const
{
  if( (id.type == type) && (id.gid == gid) ) return true;
  return false;
}

int AcmeZoltanGID::operator%( const int& divisor ) const
{
  return ( gid % divisor );
}

#endif
