// $Id$

#ifndef ContactShellTriFaceL3_h_
#define ContactShellTriFaceL3_h_

#include "ContactTriFaceL3.h"
#include <cstring>

/* This class represents the linear three node triangle shell face with the
   following fortran numbering convention.

                              2
                              /\
                             /  \
                            /    \
                        E2 /      \ E1
                          /        \
                         /          \
                        0------------1
                              E0
*/

class ContactShellTriFaceL3 : public ContactTriFaceL3<Real> {

  public :
    ContactShellTriFaceL3( ContactFixedSizeAllocator*, 
                           int blk_indx=-1, int indx_in_block=-1, int key=-1);
  static 
    ContactShellTriFaceL3* new_ContactShellTriFaceL3(ContactFixedSizeAllocator*,
                                                     int blk_indx=-1, 
                                                     int indx_in_block=-1, 
                                                     int key=-1);
  ~ContactShellTriFaceL3();

 private:

#include "shell_functions.h"

};

#endif // ContactShellTriFaceL3_h_
