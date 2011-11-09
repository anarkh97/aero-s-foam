// $Id: AcmeNode.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeNode_h_
#define _AcmeNode_h_

#include "AcmeEntity.h"

class AcmeFace;
class AcmeElem;

class AcmeNode : public AcmeEntity {

 public:

  AcmeNode( AcmeEntity::AcmeNodeType,
            int Block_Index=-1,
	    int Index_in_Block=-1,
	    int Exo_Id=-1 );
                 
  ~AcmeNode();
  
  inline AcmeEntity::AcmeNodeType NodeType() { return type; };
  
  void Coordinates(Real x, Real y, Real z) 
    { coordinates[0]=x; coordinates[1]=y; coordinates[2]=z; };
  Real* Coordinates() { return coordinates; };

  inline AcmeFace*  Face( int i ) {return( faces[i] );};
  inline AcmeFace** Faces() { return faces; };
  inline AcmeElem*  Elem( int i ) {return( elems[i] );};
  inline AcmeElem** Elems() { return elems; };

  // Functions to connect up the topology
  void Delete_Face_Connections( );
  void Allocate_Face_Connections( int number );
  void Connect_Face( int, AcmeFace* );
  void Delete_Elem_Connections( );
  void Allocate_Elem_Connections( int number );
  void Connect_Elem( int, AcmeElem* );

  inline void Number_Face_Connections(int n) {number_face_connections=n;};
  inline int  Number_Face_Connections() {return number_face_connections;};
  inline void Number_Elem_Connections(int n) {number_elem_connections=n;};
  inline int  Number_Elem_Connections() {return number_elem_connections;};
         
  // Packing/Unpacking Functions
  inline virtual int  Size( );
  inline virtual void Pack( char* );
  inline virtual void Unpack( char* );

 private:
  AcmeEntity::AcmeNodeType type;
  int exodus_id;
  int number_face_connections;
  int number_elem_connections;
  AcmeFace** faces;
  AcmeElem** elems;
  
  Real  coordinates[3];
  
};

inline int AcmeNode::Size()
{
  int size = AcmeEntity::Size() + 2*sizeof(int) + 3*sizeof(Real);
  return size;
}

inline void AcmeNode::Pack(char* buffer)
{
  char* buff = buffer;
  AcmeEntity::Pack( buff );
  buff += AcmeEntity::Size();
  
  int cnt        = 0;
  int* i_buf     = REINTERPRET_CAST(int*) (buff);
  i_buf[cnt++]   = exodus_id;
  i_buf[cnt++]   = type;
  Real* data_loc = REINTERPRET_CAST(Real*)(buff+2*sizeof(int));
  memcpy(data_loc,coordinates,3*sizeof(Real));
  
}

inline void AcmeNode::Unpack( char* buffer )
{
  char* buff = buffer;
  AcmeEntity::Unpack( buff );
  buff += AcmeEntity::Size();
  
  int cnt        = 0;
  int* i_buf     = REINTERPRET_CAST(int*) ( buff );
  exodus_id      = i_buf[cnt++];
  type           = (AcmeEntity::AcmeNodeType) i_buf[cnt++];
  Real* data_loc = REINTERPRET_CAST(Real*)(buff+2*sizeof(int));
  memcpy(coordinates,data_loc,3*sizeof(Real));
}

#endif  // #ifdef _AcmeNode_h_
