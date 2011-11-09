// $Id: AcmeZoltanID.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeZoltanID_h_
#define _AcmeZoltanID_h_

#ifndef CONTACT_NO_MPI

#include "lbi_const.h"

class AcmeZoltanLID {

 public:
                    
  AcmeZoltanLID();
  AcmeZoltanLID(LB_ID_PTR);
  ~AcmeZoltanLID();
  
  void ZoltanLID(int, int, LB_ID_PTR);
  bool operator> ( const AcmeZoltanLID& ) const;
  bool operator< ( const AcmeZoltanLID& ) const;
  bool operator==( const AcmeZoltanLID& ) const;
  int operator% ( const int& ) const;

  static int  Type(LB_ID_PTR id) {return id[0];};
  static int  Index(LB_ID_PTR id) {return id[1];};

 private:
   int type;
   int index;

};

class AcmeZoltanGID {

 public:
  enum Zoltan_ID_Type{ ZOLTAN_GID, ZOLTAN_LID };
                            
  AcmeZoltanGID();
  AcmeZoltanGID(LB_ID_PTR);
  ~AcmeZoltanGID();
  
  void ZoltanGID(int, int, LB_ID_PTR);
  bool operator> ( const AcmeZoltanGID& ) const;
  bool operator< ( const AcmeZoltanGID& ) const;
  bool operator==( const AcmeZoltanGID& ) const;
  int operator% ( const int& ) const;

  static int  Type(LB_ID_PTR id) {return id[0];};
  static int  ID(LB_ID_PTR id) {return id[1];};

 private:
   int type;
   int gid;

};

#endif

#endif
