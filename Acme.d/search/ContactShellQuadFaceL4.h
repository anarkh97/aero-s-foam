// $Id$

#ifndef ContactShellQuadFaceL4_h_
#define ContactShellQuadFaceL4_h_

#include "ContactQuadFaceL4.h"
#include <cstring>

class ContactFixedSizeAllocator;

/* This class represents the bilinear four node quadrilateral shell face with 
   the following fortran numbering convention.

                             E2
                       3-------------2
                       |             |
                       |             |
                    E3 |             | E1
                       |             |
                       0-------------1
                             E0
*/

class ContactShellQuadFaceL4 : public ContactQuadFaceL4<Real> {

  public :
    ContactShellQuadFaceL4( ContactFixedSizeAllocator*, 
                            int blk_indx=-1, int indx_in_block=-1, int key=-1);
  static ContactShellQuadFaceL4* new_ContactShellQuadFaceL4(
		    ContactFixedSizeAllocator*,
                    int blk_indx=-1, int indx_in_block=-1, int key=-1);

  ~ContactShellQuadFaceL4();

  private :
    
#include "shell_functions.h"

};

#endif // ContactShellQuadFaceL4_h_
