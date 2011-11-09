// $Id: output_interactions.C,v 2002.4 2003/12/24 22:09:58 mwglass Exp $

#include "Contact_Defines.h"
#include "output_interactions.h"
#include "ContactParOStream.h"
#include <iostream.h>
#ifndef CONTACT_NO_MPI
#include <mpi.h>
#endif

void FORTRAN(output_nodenode_interactions)( int& num_interactions,
					    int* slave_node_block_id,
					    int* slave_node_indexes_in_block,
					    int* master_node_block_id,
					    int* master_node_indexes_in_block,
					    int* master_node_proc,
					    Real* interaction_data,
                                            int& data_size)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=3
#ifndef CONTACT_NO_MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  ContactParOStream postream( comm );
#else
  int dummy = 0;
  ContactParOStream postream( dummy );
#endif
  postream<<"***** Host Code Node-Node Interactions *****\n";
  if (num_interactions==0) {
    postream<<"** NONE **\n";
  }
  for( int i=0 ; i<num_interactions ; i++ ){
    postream << "Node-Node Interaction " << i+1 << "\n";
    postream << "   Slave Node:   Block = " << slave_node_block_id[i] 
             << "   Index = " << slave_node_indexes_in_block[i] << "\n";
    postream << "   Master Node:  Block = " << master_node_block_id[i] 
             << "   Index = " << master_node_indexes_in_block[i] 
	     << "    Proc = " << master_node_proc[i] << "\n";
    postream << "   Distance: " << interaction_data[data_size*i] << "\n";
  }
  postream.flush();
#endif
}				    

void FORTRAN(output_nodeface_interactions)( int& num_interactions,
					    int* node_block_id,
					    int* node_indexes_in_block,
					    int* face_block_id,
					    int* face_indexes_in_block,
					    int* face_proc,
					    Real* interaction_data,
                                            int& data_size)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=3
#ifndef CONTACT_NO_MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  ContactParOStream postream( comm );
#else
  int dummy = 0;
  ContactParOStream postream( dummy );
#endif
  postream<<"***** Host Code Node-Face Interactions *****\n";
  if (num_interactions==0) {
    postream<<"** NONE **\n";
  }
  for( int i=0 ; i<num_interactions ; i++ ){
    postream << "Node-Face Interaction " << i+1 << "\n";
    postream << "   Node:  Block = " << node_block_id[i] << "   Index = " 
	     << node_indexes_in_block[i] << "\n";
    postream << "   Face:  Block = " << face_block_id[i] << "   Index = "
	     << face_indexes_in_block[i] << "    Proc = " 
	     << face_proc[i] << "\n";
    postream << "   Contact Point: " << interaction_data[data_size*i+0] << " "
	     << interaction_data[data_size*i+1] << "\n";
    postream << "   Gap: " << interaction_data[data_size*i+2] << "\n";
    postream << "   Interaction Vector: " << interaction_data[data_size*i+3] << " "
	     << interaction_data[data_size*i+4] << " "
	     << interaction_data[data_size*i+5] << "\n";
    postream << "   Surface Normal: " << interaction_data[data_size*i+6] << " "
	     << interaction_data[data_size*i+7] << " " 
	     << interaction_data[data_size*i+8] << "\n";
    postream << "   Algorithm: " << interaction_data[data_size*i+9] << "\n";
  }
  postream.flush();
#endif
}				    

void FORTRAN(output_nodesurface_interactions)( int& num_interactions,
						int* n_blk_id,
						int* n_ind,
						int* surf_id,
						Real* data, int *data_size )
{
#if CONTACT_DEBUG_PRINT_LEVEL>=3
#ifndef CONTACT_NO_MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  ContactParOStream postream( comm );
#else
  int dummy = 0;
  ContactParOStream postream( dummy );
#endif
  postream<<"***** Host Code Node-Surface Interactions *****\n";
  if (num_interactions==0) {
    postream<<"** NONE **\n";
  }
  int n = *data_size;
  for( int i=0 ; i<num_interactions ; i++ ){
    int nn = n*i;
    postream << "Node-Surface_Interaction " << i << "\n";
    postream << "   Node:  Block = " << n_blk_id[i] << "   Index = " 
	     << n_ind[i] << "\n";
    postream << "   Surface ID = " << surf_id[i] << "\n";
    postream << "   Gap = " << data[nn+3] << "\n";
    postream << "   Contact Point: " 
             << data[nn+0] << " " 
	     << data[nn+1] << " "
	     << data[nn+2] << "\n";
    postream << "   Surface Normal: " 
             << data[nn+4] << " "
	     << data[nn+5] << " " 
	     << data[nn+6] << "\n";
  }
  postream.flush();
#endif
}

void FORTRAN(print_message)(char* message)
{
  cout << message << endl;
}
