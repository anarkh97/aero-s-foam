// $Id: output_interactions.h,v 2002.4 2003/12/24 22:09:58 mwglass Exp $

#include "Contact_Defines.h"
#include <iostream.h>

#ifndef _output_interactions_h_
#define _output_interactions_h_

extern "C" {

  void FORTRAN(output_nodenode_interactions)( int& num_interactions,
					      int* slave_node_block_id,
					      int* slave_node_indexes_in_block,
					      int* master_node_block_id,
					      int* master_node_indexes_in_block,
					      int* master_node_proc,
					      Real* interaction_Data,
                                              int& data_size);

  void FORTRAN(output_nodeface_interactions)( int& num_interactions,
					      int* node_block_id,
					      int* node_indexes_in_block,
					      int* face_block_id,
					      int* face_indexes_in_block,
					      int* face_proc,
					      Real* interaction_Data,
                                              int& data_size);

  void FORTRAN(output_nodesurface_interactions)( int& num_interactions,
						 int* n_blk_id,
						 int* n_ind,
						 int* surf_id,
						 Real* data, 
                                                 int* data_size );
  
  void FORTRAN(print_message)( char* message );
  
}						

       
#endif
