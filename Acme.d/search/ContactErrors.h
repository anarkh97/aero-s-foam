// $Id$

#ifndef ContactErrors_h_
#define ContactErrors_h_

#include "Contact_Defines.h"
#include "CString.h"

class ContactErrors{

 public:
  ContactErrors();
  ~ContactErrors();

  void Add_Error_Message(const char* message);

  inline int Number_of_Errors() {return number_of_errors;};
  const char* Error_Message( int i );

 private:
  
  ContactErrors(ContactErrors&);
  ContactErrors& operator=(ContactErrors&);
  
  int errors_array_size;
  int number_of_errors;
  CString** errors;

  void Reallocate();

};

#endif
