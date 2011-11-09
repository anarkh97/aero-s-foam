// $Id: usersubs.C,v 2002.3 2004/05/10 16:17:04 mwglass Exp $

#include "usersubs.h"
#include <string.h>

extern "C" {

int          nsubs = 0;
UserSubsList user_subs[100];

void FORTRAN(reg_usersub)( void* user_sub ,
                           const char* label ,
                           int&        label_length )
{
  int i=0;
  for (i=0; i<nsubs; i++) {
    if (strcmp(label,user_subs[i].label)==0) {
    }
  }
  user_subs[i].fn    = user_sub;
  user_subs[i].label = new char[label_length+1];
  strncpy(user_subs[i].label,label,label_length);
  user_subs[i].label[label_length] = (char)NULL;
  nsubs++;
}

void FORTRAN(get_usersub)( void** user_sub ,
                           const char* label ,
                           int&        label_length )
{
  int i;
  for (i=0; i<nsubs; i++) {
    if (strcmp(label,user_subs[i].label)==0) break;
  }
  *user_sub = user_subs[i].fn;
}

void FORTRAN(str_len) (const char* string,
                             int& length)
{ length = strlen(string); 
/*
  int i =0;
  while (string[i] != '\0' && i < 256) i++;
  length = i;
*/
}

}
