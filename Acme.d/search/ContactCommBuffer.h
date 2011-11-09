// $Id$

#ifndef ContactCommBuffer_h_
#define ContactCommBuffer_h_

#include "Contact_Defines.h"
#include "Contact_Communication.h"

#ifndef CONTACT_NO_MPI
class ContactSymComm;
class ContactAsymComm;
#endif

class ContactCommBuffer {
  
 public:

  ContactCommBuffer(); 
  ~ContactCommBuffer();

#ifndef CONTACT_NO_MPI
  void Buffers( ContactSymComm&,    int, char**, char**, 
                RequestHandle**, RequestHandle**);
  void Buffers( ContactAsymComm&,   int, char**, char**, RequestHandle** );
  void Buffers_Export( ContactAsymComm&,   int, char**, char**, 
                       RequestHandle**, RequestHandle** );
  void Buffers_Import( ContactAsymComm&,   int, char**, char**, 
                       RequestHandle**, RequestHandle** );
  void Buffers( int size_send_buf, int size_recv_buff, int size_request_handle,
		char** send_buf, char** recv_buf, 
		RequestHandle** request_handles );
#endif

 private:
  
  ContactCommBuffer(ContactCommBuffer&);
  ContactCommBuffer& operator=(ContactCommBuffer&);

  int buffer_size;
  char* buffer;
  int handles_size;
  RequestHandle* handles;
};

#endif
