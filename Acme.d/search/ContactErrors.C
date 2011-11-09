// $Id$

#include "ContactErrors.h"
#include <cstring>

ContactErrors::ContactErrors()
{
  number_of_errors = 0;
  errors_array_size = 10;
  errors = new CString*[errors_array_size];
  std::memset( errors, 0, errors_array_size*sizeof(CString*) );
}

ContactErrors::~ContactErrors()
{
  for( int i=0 ; i<number_of_errors ; ++i )
    delete errors[i];
  delete [] errors;
}


const char* ContactErrors::Error_Message( int i )
{
  return errors[i]->data();
}

void ContactErrors::Add_Error_Message( const char* Message )
{
  if( number_of_errors == errors_array_size ) Reallocate();
  // The last error message wasn't deleted in Error_Message so delete it
  // here to avoid a memory lead (if needed)
  errors[number_of_errors++] = new CString(Message);
}

void ContactErrors::Reallocate()
{
  CString** temp = errors;
  errors_array_size += 10;
  errors = new CString*[errors_array_size];
  std::memset( errors, 0, errors_array_size*sizeof(CString*) );
  for( int i=0 ; i<number_of_errors ; ++i )
    errors[i] = temp[i];
  delete [] temp;
}
