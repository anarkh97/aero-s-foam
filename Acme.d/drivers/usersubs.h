// $Id: usersubs.h,v 2002.2 2004/03/30 18:08:45 rjones Exp $

#ifndef _usersub_h_
#define _usersub_h_

typedef struct {
  void *fn;
  char *label;
} UserSubsList;

#include "Contact_Defines.h"

extern "C" {

void FORTRAN(reg_usersub)( void* user_sub,
                           const char* global_label,
                           int&        global_label_length );

void FORTRAN(get_usersub)( void** user_sub,
                           const char* global_label,
                           int&        global_label_length );

void FORTRAN(str_len) (const char* string,
                               int& length);

}
#endif
