// $Id$
  
#include <ContactBoundingBox.h>
#include <Contact_Communication.h>
#include "ContactTimer.h"
#include "CString.h"
#include "contact_sorting.h"
#include <iostream>
#include "cmath"
#include "ContactParOStream.h"

using namespace std;
//
//  Combine a set of boxes that exist on several processors into a single
//  box that encompasses all of the object boxes
//
void ObjectBoundingBox::global_box_combine(
                           ObjectBoundingBox *box_array, 
                           const int &num_boxes, 
                           MPI_Comm &communicator) {
  //
  //  Allocate a common set of arrays to perform the reductions on
  //
  int array_length = num_boxes * 3;
  Real *all_box_min_local = new Real[array_length];
  Real *all_box_max_local = new Real[array_length];
  Real *all_box_min_global = new Real[array_length];
  Real *all_box_max_global = new Real[array_length];
  //
  //  Fill the local arrays
  //
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    all_box_min_local[ibox * 3 + 0] = box_array[ibox].get_x_min();
    all_box_min_local[ibox * 3 + 1] = box_array[ibox].get_y_min();
    all_box_min_local[ibox * 3 + 2] = box_array[ibox].get_z_min();
    all_box_max_local[ibox * 3 + 0] = box_array[ibox].get_x_max();
    all_box_max_local[ibox * 3 + 1] = box_array[ibox].get_y_max();
    all_box_max_local[ibox * 3 + 2] = box_array[ibox].get_z_max();
  }
  //
  //  Perform the global MPI reductions
  //
  contact_global_minimum(all_box_min_local, 
                         all_box_min_global, 
                         array_length, 
                         communicator);
  contact_global_maximum(all_box_max_local, 
                         all_box_max_global, 
                         array_length, 
                         communicator);
  //
  //  Scatter the local arrays back to the boxes
  //
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    box_array[ibox].set_x_min(all_box_min_global[ibox * 3 + 0]);
    box_array[ibox].set_y_min(all_box_min_global[ibox * 3 + 1]);
    box_array[ibox].set_z_min(all_box_min_global[ibox * 3 + 2]);
    box_array[ibox].set_x_max(all_box_max_global[ibox * 3 + 0]);
    box_array[ibox].set_y_max(all_box_max_global[ibox * 3 + 1]);
    box_array[ibox].set_z_max(all_box_max_global[ibox * 3 + 2]);
  }
  delete [] all_box_min_local;
  delete [] all_box_max_local;
  delete [] all_box_min_global;
  delete [] all_box_max_global;
}

int ContactBoundingBox::longest_dimension() {
  Real length_0 = x_max-x_min;
  Real length_1 = y_max-y_min;
  Real length_2 = z_max-z_min;
  if(length_0 > length_1 && length_0 > length_2) {
    return 0;
  } else if(length_1 > length_2) {
    return 1;
  } else {
    return 2;
  }
}

ContactParOStream& operator<<( ContactParOStream& pos, 
                               const ContactBoundingBox &box) {
  pos <<"Min corner "<<box.get_x_min()<<" "<<box.get_y_min()<<" "<<box.get_z_min();
  pos <<" Max corner "<<box.get_x_max()<<" "<<box.get_y_max()<<" "<<box.get_z_max();
  return pos;
}

