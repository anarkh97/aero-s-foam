// $Id$

#ifndef ContactZoltanID_h_
#define ContactZoltanID_h_

#ifndef CONTACT_NO_MPI

#include "lbi_const.h"

class ContactHostGlobalID;

class ContactZoltanLID {

 public:
                    
  ContactZoltanLID();
  ContactZoltanLID(LB_ID_PTR);
  ~ContactZoltanLID();
  
  void ZoltanLID(int, int, LB_ID_PTR);
  bool operator> ( const ContactZoltanLID& ) const;
  bool operator< ( const ContactZoltanLID& ) const;
  bool operator==( const ContactZoltanLID& ) const;
  int  operator% ( const int& ) const;

  static int  Type(LB_ID_PTR id)  {return id[0];};
  static int  Index(LB_ID_PTR id) {return id[1];};

 private:
   int type;
   int index;
};

class ContactZoltanGID {

 public:
  enum Zoltan_ID_Type{ ZOLTAN_GID, ZOLTAN_LID };
                            
  ContactZoltanGID();
  ContactZoltanGID(LB_ID_PTR);
  ~ContactZoltanGID();
  
  void ZoltanGID(int, int, int, LB_ID_PTR);
  void ZoltanGID(int, ContactHostGlobalID*, LB_ID_PTR);
  bool operator> ( const ContactZoltanGID& ) const;
  bool operator< ( const ContactZoltanGID& ) const;
  bool operator==( const ContactZoltanGID& ) const;
  int  operator% ( const int& ) const;

  static int  Type(LB_ID_PTR id) {return id[0];};
  static int  Hi(LB_ID_PTR id)   {return id[1];};
  static int  Lo(LB_ID_PTR id)   {return id[2];};
  
  int Type() { return type; };
  int Hi()   { return hi; };
  int Lo()   { return lo; };

 private:
   int type;
   int hi;
   int lo;

};

inline
ContactZoltanLID::ContactZoltanLID(LB_ID_PTR id)
  : type(id[0]), index(id[1]) {}

inline void 
ContactZoltanLID::ZoltanLID(int id_type, int indx, LB_ID_PTR id)
{
  id[0] = id_type;
  id[1] = indx;
}

inline bool 
ContactZoltanLID::operator>( const ContactZoltanLID& id ) const
{
  return type != id.type ? type > id.type : index > id.index;
  //if( type > id.type )
  //  return true;
  //else if( type < id.type )
  //  return false;
  //if( index > id.index )
  //  return true;
  //return false;
}

inline bool 
ContactZoltanLID::operator<( const ContactZoltanLID& id ) const
{
  return type != id.type ? type < id.type : index < id.index;
  //if( type < id.type )
  //  return true;
  //else if( type > id.type )
  //  return false;
  //if( index < id.index )
  //  return true;
  //return false;
}

inline bool 
ContactZoltanLID::operator==( const ContactZoltanLID& id ) const
{
  if( (id.type == type) && (id.index == index) ) return true;
  return false;
}

inline int 
ContactZoltanLID::operator%( const int& divisor ) const
{
  return ( index % divisor );
}

inline
ContactZoltanGID::ContactZoltanGID(LB_ID_PTR id)
  : type(id[0]), hi(id[1]), lo(id[2]) {}

inline void 
ContactZoltanGID::ZoltanGID(int id_type, int hi_int, 
                            int lo_int, LB_ID_PTR id)
{
  id[0] = id_type;
  id[1] = hi_int;
  id[2] = lo_int;
}

//inline void 
//ContactZoltanGID::ZoltanGID(int id_type, ContactHostGlobalID* gid, 
//                            LB_ID_PTR id)
//{
//  id[0] = id_type;
//  id[1] = gid->HiInt();
//  id[2] = gid->LoInt();
//}

inline bool 
ContactZoltanGID::operator>( const ContactZoltanGID& id ) const
{
  if( type > id.type )
    return true;
  else if( type < id.type )
    return false;
  return hi != id.hi ? hi > id.hi : lo > id.lo;
  //if( hi > id.hi )
  //  return true;
  //else if( hi < id.hi )
  //  return false;
  //if( lo > id.lo )
  //  return true;
  //return false;
}


inline bool 
ContactZoltanGID::operator<( const ContactZoltanGID& id ) const
{
  if( type < id.type )
    return true;
  else if( type > id.type )
    return false;
  return hi != id.hi ? hi < id.hi : lo < id.lo;
  //if( hi < id.hi )
  //  return true;
  //else if( hi > id.hi )
  //  return false;
  //if( lo < id.lo )
  //  return true;
  //return false;
}


inline bool 
ContactZoltanGID::operator==( const ContactZoltanGID& id ) const
{
  if( (id.type == type) && (id.hi == hi) && (id.lo == lo) ) return true;
  return false;
}

inline int 
ContactZoltanGID::operator%( const int& divisor ) const
{
  return ( lo % divisor );
}

#endif

#endif
