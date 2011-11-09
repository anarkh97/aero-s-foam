// $Id: make_unit_test.C,v 2002.1 2002/01/04 16:24:31 khbrown Exp $

// This program reads a plot file output from the contact search and
// creates a unit test from it that can then be run from the drivers.

#include "ContactSearch.h"
#include "Contact_Defines.h"

#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <exodusII.h>
#include <string.h>

#ifdef CONTACT_MPI
#include "ne_nemesisI.h"
#endif

int process_file(char*, char*, char*, int, int);

int main(int argc, char*argv[] )
{
  if( argc != 4 && argc !=5){
    cout << "\n\nUsage: make_unit_test contact_mesh.exo mesh_inp.exo "
	 << "mesh.data [num_procs] \n";
    exit(1);
  }
  
  int num_procs;
  char* Exodus_in_File = argv[1];
  char* Exodus_out_File = argv[2];
  char* UnitTest_File = argv[3];
  if (5 == argc)
    num_procs = atoi(argv[4]);
  else
    num_procs = 0;

  // if there is more than one exodus file, then loop over each file

  if ( num_procs > 0 ) {
    for (int i = 0; i < num_procs; i++) {
      char par_in_file_name[1024];
      char par_out_file_name[1024];
      char trans[81];
      char trans_tmp[81];
      strcpy(par_in_file_name,Exodus_in_File);
      strcpy(par_out_file_name,Exodus_out_File);
      sprintf(trans_tmp,".%d.%d",num_procs,i);
      sprintf(trans,".%d.",num_procs);
      strcat(par_out_file_name,trans);
      if (num_procs >= 1000) {
 	sprintf(trans,"%04d",i);
      }
      else if (num_procs >= 100) {
 	sprintf(trans,"%03d",i);
      }
      else if (num_procs >= 10) {
 	sprintf(trans,"%02d",i);
      }
      else {
	sprintf(trans,"%d",i);
      }
      strcat(par_in_file_name,trans_tmp);
      strcat(par_out_file_name,trans);
      
      process_file(par_in_file_name, par_out_file_name,
		   UnitTest_File, num_procs, i);
    }
  }
  else {
    process_file(Exodus_in_File, Exodus_out_File, UnitTest_File, 0, 0);
  }
  
}

int process_file(char* Exodus_in_File, char* Exodus_out_File,
		 char* UnitTest_File, int num_procs, int my_proc_id) {

  int i;

  // Open the exodus file
  int cpu_ws = 8;
  int io_ws = 0;
  float ver;
  char Title[81];
  int Exodus_in_ID = ex_open( Exodus_in_File, EX_READ, &cpu_ws, &io_ws, &ver );
  if( Exodus_in_ID < 0 ){
    cout << "\n\n  ERROR: Opening Exodus Input File: " << Exodus_in_File << "\n\n";
    exit(1);
  }

  // Get the size of the problem
  int Number_of_Dimensions;
  int Number_of_Nodes;
  int Number_of_Elements;
  int Number_of_Element_Blocks;
  int Number_of_Node_Sets;
  int Number_of_Side_Sets;
  if( ex_get_init( Exodus_in_ID, Title, &Number_of_Dimensions , 
		   &Number_of_Nodes ,
                   &Number_of_Elements , &Number_of_Element_Blocks ,
                   &Number_of_Node_Sets , &Number_of_Side_Sets ) < 0 ) {
    cerr << "\n\n**ERROR** EX_GET_INIT \n";
    exit(1);
  }

  // Read the coordinates

  Real* x_coord = new Real[Number_of_Nodes];
  Real* y_coord = new Real[Number_of_Nodes];
  Real* z_coord = new Real[Number_of_Nodes];

  if( ex_get_coord( Exodus_in_ID, x_coord, y_coord, z_coord ) ){
    cout << "Error Reading coordinates" << endl;
    exit(1);
  }

  // Read the displacements and compute the predicted configuration

  Real* x_disp = new Real[Number_of_Nodes];
  Real* y_disp = new Real[Number_of_Nodes];
  Real* z_disp = new Real[Number_of_Nodes];

  if( ex_get_nodal_var( Exodus_in_ID, 1, 1, Number_of_Nodes, x_disp ) ){
    cout << "Error Reading displacements" << endl;
    exit(1);
  }
  if( ex_get_nodal_var( Exodus_in_ID, 1, 2, Number_of_Nodes, y_disp ) ){
    cout << "Error Reading displacements" << endl;
    exit(1);
  }
  if( Number_of_Dimensions == 3){
    if( ex_get_nodal_var( Exodus_in_ID, 1, 3, Number_of_Nodes, z_disp ) ){
      cout << "Error Reading displacements" << endl;
      exit(1);
    }
  }

  // read in the info data -- number of face blocks and search data
  // NOTE: here we assume that the info records first hold the number of
  //   face blocks, then the search data.  If this changes, then 
  //   this line will not work.
  char* cdum = NULL;
  Real fdum;
  int num_info_records;
  ex_inquire( Exodus_in_ID, EX_INQ_INFO, &num_info_records, &fdum, cdum );
  int num_entity_keys = (int) sqrt(double((num_info_records-1)/ContactSearch::NSIZSD));
  if( ! num_info_records ){
    cout << "Error reading info records" << endl;
    exit(1);
  }
  char** info_strings = new char*[num_info_records];
  int info_string_len = 81;
  for( i=0 ; i<num_info_records ; i++ )
    info_strings[i] = new char[info_string_len];
  ex_get_info( Exodus_in_ID, info_strings );
  int num_face_blocks = atoi(info_strings[0]);
  
  // read in the faces. Get rid of bar elements for now.
  //  NOTE: THIS WILL NEED TO BE CHANGED WHEN WE GO TO FACE/EDGE CONTACT
  int number_of_elems_per_block;
  int nodes_per_element,attributes_per_element;
  char element_name[81];
  
  const int char_size = 81;
  int Total_Number_of_Faces = 0;
  int * num_faces_per_block = new int[num_face_blocks];
  char ** face_names = new char*[num_face_blocks];
  for (i = 0; i < num_face_blocks; i++) 
    face_names[i] = new char[char_size];
  int * nodes_per_face = new int[num_face_blocks];
  int * attributes_per_face = new int[num_face_blocks];
  int * face_block_index = new int[num_face_blocks];
  
  for( i=0 ; i<num_face_blocks ; i++ ){
    if( ex_get_elem_block( Exodus_in_ID, i+1, element_name,
			   &number_of_elems_per_block, &nodes_per_element,
			   &attributes_per_element )){
      cout << "Error Reading Element Block Parameters.\n" << endl;
      exit(1);
    }
    
    // copy element data into face structures we will output
    for ( int j = 0; j < 81; j++)
      face_names[i][j] = element_name[j];
    num_faces_per_block[i] = number_of_elems_per_block;
    nodes_per_face[i] = nodes_per_element;
    attributes_per_face[i] = attributes_per_element;
    Total_Number_of_Faces += num_faces_per_block[i];
    face_block_index[i] = i+1;
  }
  
  //
  //  Now create and write the new exodus file
  //

  // create new file
  int Exodus_out_ID = ex_create( Exodus_out_File, EX_CLOBBER, &cpu_ws, &io_ws);
  if( Exodus_out_ID < 0 ){
    cout << "\n\n  ERROR: Creating Exodus Output File: " << Exodus_out_File 
	 << "\n\n";
    exit(1);
  }

  // write header info
  int num_node_sets = 0;
  int num_side_sets = 0;
  if (ex_put_init( Exodus_out_ID, Title, Number_of_Dimensions, 
		   Number_of_Nodes, Total_Number_of_Faces, 
		   num_face_blocks, num_node_sets, num_side_sets ) < 0 ) {
    cout << "\n\n  ERROR: Writing Input File Header " << "\n\n";
    exit(1);
  }
  
  // Output the coordinate names
  char* coord_names[3];
  coord_names[0] = (char*) "X";
  coord_names[1] = (char*) "Y";
  coord_names[2] = (char*) "Z";
  if ( ex_put_coord_names( Exodus_out_ID, coord_names ) <0 ) {
    cout << "ERROR in writing coordinate names" << endl;
    exit(1);
  }

  // write coords
  if( ex_put_coord( Exodus_out_ID, x_coord, y_coord, z_coord ) < 0 ){
    cout << "Error Writing coordinates" << endl;
    exit(1);
  }

  // write displacements
  if( ex_put_var_param( Exodus_out_ID, "n", 3 ) < 0 ) {
    cout << "Error setting up variables" << endl;
    exit(1);
  }

  char* var_names[3];
  var_names[0] = (char *) "dispx";
  var_names[1] = (char *) "dispy";
  var_names[2] = (char *) "dispz";

  if( ex_put_var_names( Exodus_out_ID, (char *) "n", 3, var_names) < 0 ){
    cout << "Error registering parameters" << endl;
    exit(1);
  }

  Real time = 1.0;
  if( ex_put_time( Exodus_out_ID, 1, &time) < 0 ){
    cout << "Error setting up time" << endl;
    exit(1);
  }

  if( ex_put_nodal_var( Exodus_out_ID, 1, 1, Number_of_Nodes, x_disp ) < 0 ){
    cout << "Error Writing displacements" << endl;
    exit(1);
  }
  if( ex_put_nodal_var( Exodus_out_ID, 1, 2, Number_of_Nodes, y_disp ) < 0 ){
    cout << "Error Writing displacements" << endl;
    exit(1);
  }
  if( Number_of_Dimensions == 3){
    if( ex_put_nodal_var( Exodus_out_ID, 1, 3, Number_of_Nodes, z_disp ) < 0 ){
      cout << "Error Writing displacements" << endl;
      exit(1);
    }
  }

  // write face blocks
  for( i=0 ; i<num_face_blocks ; i++ ){
    if( ex_put_elem_block( Exodus_out_ID, i+1, face_names[i],
			   num_faces_per_block[i], nodes_per_face[i],
			   attributes_per_face[i] )){
      cout << "Error Writing Element Block Parameters.\n" << endl;
      exit(1);
    }
  }
  
  // read and write connectivity
  int* connectivity = new int[8*Number_of_Elements];
  for( i=0 ; i<num_face_blocks ; i++ ){
    if (num_faces_per_block[i]>0){
      if( ex_get_elem_conn( Exodus_in_ID, face_block_index[i], connectivity) ){
	cout << "Error Reading Element Block Connectivity" << endl;
	exit(1);
      }
      if( ex_put_elem_conn( Exodus_out_ID, i+1, connectivity ) ){
	cout << "Error writing Element Block Connectivity" << endl;
	exit(1);
      }
    }
  }

  // write qa records
  char *qa_record[1][4];
  qa_record[0][0]=(char *) "Contact Test";
  qa_record[0][1]=(char *) "contact test";
  qa_record[0][2]=(char *) "no date recorded";
  qa_record[0][3]=(char *) "no time recorded";
  if( ex_put_qa( Exodus_out_ID, 1, qa_record) < 0 ){
    cout << "Error setting up qa records" << endl;
    exit(1);
  }

  // write the output file if this is the first processor (or only processor)
  if ( 0 == my_proc_id ) {

    // Open the output file
    ofstream UnitTest( UnitTest_File, ios::out );

    // write header info
    UnitTest << "Number of Load Steps" 
	     << endl;
    UnitTest << "1" << endl;
    
    UnitTest << "Search Type (0=static 1 state, 1= static 2 state, 2=dynamic)" 
	     << endl;
    UnitTest << "1" << endl;
    
    UnitTest << "Restart Flag (0=off,1=read only,2=write only,3=read/write)" 
	     << endl;
    UnitTest << "0" << endl;
    
    // write analytic surfaces and search data to output file
    UnitTest << "Number of Analytic Surfaces" << endl;
    UnitTest << "0" << endl;
    
    // write number of entity keys to output file
    UnitTest << "Number of Entity Keys" << endl;
    UnitTest << num_entity_keys << endl;
    
    // Output the search Data from the info records read earlier -- skip 
    // first entry, because it is just the number of face blocks 
    UnitTest << "Search Data" << endl;
    for( i=1 ; i<num_info_records ; i++ )
      UnitTest << info_strings[i] << endl;
    
    // Get global variables for multiple interaction and normal
    // smoothing status and data
    int num_glob_vars;
    if ( ex_get_var_param (Exodus_in_ID, (char *) "g", &num_glob_vars)) {
      cout << "ERROR in reading global variables" << endl;
      exit(1);
    }
    if ( num_glob_vars != 9 ) {
      cout << "ERROR: make_unit_test knows about 9 global variables, " << endl
	   << "   but the given exodusII file has " << num_glob_vars << endl;
      exit(1);
    }
    char** glob_var_names = new char*[num_glob_vars];
    for (i = 0; i< num_glob_vars; i++) {
      char * temp = new char[MAX_STR_LENGTH];
      glob_var_names[i] = temp;
    }
    Real * glob_var_vals = new Real[num_glob_vars];
    if ( ex_get_var_names (Exodus_in_ID, (char *) "g", num_glob_vars, 
			   glob_var_names)) {
      cout << "ERROR in reading global variable names" << endl;
      exit(1);
    }
    if ( ex_get_glob_vars (Exodus_in_ID, 1, num_glob_vars, glob_var_vals)) {
      cout << "ERROR in reading global variable names" << endl;
      exit(1);
    }
    if ( strcmp(glob_var_names[4],"mult_interaction_status") ) {
      cout << "ERROR: cannot find global variable mult_interaction_status"
	   << endl;
      exit(1);
    }
    if ( strcmp(glob_var_names[5], "norm_smoothing_status" )) {
      cout << "ERROR: cannot find global variable norm_smoothing_status"
	   << endl;
      exit(1);
    }
    if ( strcmp(glob_var_names[6], "smoothing_angle" )) {
      cout << "ERROR: cannot find global variable smoothing_angle"
	   << endl;
      exit(1);
    }
    if ( strcmp(glob_var_names[7], "smoothing_length" )) {
      cout << "ERROR: cannot find global variable smoothing_length"
	   << endl;
      exit(1);
    }
    if ( strcmp(glob_var_names[8], "smoothing_resolution" )) {
      cout << "ERROR: cannot find global variable smoothing_resolution"
	   << endl;
      exit(1);
    }
    
    UnitTest << "Multiple Interaction Data (1=on)" << endl;
    UnitTest << (int) glob_var_vals[4] << endl;
    if ( glob_var_vals[4] != 0.0 ) {
      UnitTest << "  Sharp-Smooth Angle" << endl;
      UnitTest << "  " << glob_var_vals[3]<< endl;
    }
    UnitTest << "Normal Smoothing Option (0=off)" << endl;
    UnitTest << (int) glob_var_vals[5] << endl;
    if ( glob_var_vals[5] != 0.0 ) {
      UnitTest << "  Sharp-Smooth Angle" << endl;
      UnitTest << "  " << glob_var_vals[6]<< endl;
      UnitTest << "  Smoothing Length" << endl;
      UnitTest << "  " << glob_var_vals[7]<< endl;
      UnitTest << "  Smoothing Resolution" << endl;
      UnitTest << "  " << (int) glob_var_vals[8]<< endl;
    }
    
    for (i = 0; i < num_glob_vars; i++) 
      delete [] glob_var_names[i];
    delete [] glob_var_names;
    delete [] glob_var_vals;

    UnitTest.close();
  }    

#ifdef CONTACT_MPI
  // if we are in parallel, then copy over the communication lists
  if (num_procs > 1 ) {

    // initialize the nemesis data
    char* file_type = (char*)"p";
    ne_put_init_info( Exodus_out_ID, num_procs, 1, file_type);
    int num_nodes_global, num_elems_global, num_elem_blk_global, 
      num_node_sets_global, num_side_sets_global;
    ne_get_init_global( Exodus_in_ID, &num_nodes_global, 
			&num_elems_global, &num_elem_blk_global, 
			&num_node_sets_global, &num_side_sets_global);

    // global element block information. exclude edge elements
    int new_num_elems_global = 0;
    int new_num_elem_blk_global = 0;
    int* block_ids = new int[num_elem_blk_global];
    int* global_counts = new int[num_elem_blk_global];
    int* new_block_ids = new int[num_elem_blk_global];
    int* new_global_counts = new int[num_elem_blk_global];
    ne_get_eb_info_global( Exodus_in_ID, block_ids, global_counts );
    for ( i = 0; i < num_elem_blk_global; i++) {
      int block_id = block_ids[i];
      int found = 0;
      for ( j = 0; j < num_face_blocks; j++) {
	if (face_block_index[j] == block_id) {
	  found = 1;
	  break;
	} 
      }
      if (found) {
	new_block_ids[new_num_elem_blk_global] = block_ids[i];
	new_global_counts[new_num_elem_blk_global] = global_counts[i];
	new_num_elems_global += global_counts[i];
	new_num_elem_blk_global++;
      }
    }
    ne_put_init_global( Exodus_out_ID, num_nodes_global, 
			new_num_elems_global, new_num_elem_blk_global, 
			num_node_sets_global, num_side_sets_global);

    ne_put_eb_info_global( Exodus_out_ID, new_block_ids, new_global_counts );
    delete [] block_ids;
    delete [] global_counts;

    // decomposition parameters -- subtract off the edge elements out
    // of the number of internal elements 
    int num_internal_nodes, num_border_nodes, num_external_nodes;
    int num_internal_elems, num_border_elems, num_elem_comm_procs;
    int num_comm_procs;
    ne_get_loadbal_param( Exodus_in_ID, &num_internal_nodes, &num_border_nodes,
			  &num_external_nodes, &num_internal_elems,
			  &num_border_elems, &num_comm_procs,
			  &num_elem_comm_procs, my_proc_id );
    int new_num_internal_elems = 
      num_internal_elems - (Number_of_Elements - Total_Number_of_Faces);
    ne_put_loadbal_param( Exodus_out_ID, num_internal_nodes, num_border_nodes,
			  num_external_nodes, new_num_internal_elems,
			  num_border_elems, num_comm_procs,
			  num_elem_comm_procs, my_proc_id );

    // node maps
    int* node_mapi = new int[num_internal_nodes];
    int* node_mapb = new int[num_border_nodes];
    int* node_mape = new int[num_external_nodes];
    ne_get_node_map( Exodus_in_ID, node_mapi, node_mapb, node_mape, 
		     my_proc_id );
    ne_put_node_map( Exodus_out_ID, node_mapi, node_mapb, node_mape, 
		     my_proc_id );
    delete [] node_mapi;
    delete [] node_mapb;
    delete [] node_mape;

    // element maps
    //    NOTE: here we assume that the last elements in the element map
    //     correspond to the edge elements we are trying to strip out
    int* elem_mapi = new int[num_internal_elems];
    int* new_elem_mapi = new int[new_num_internal_elems];
    int* elem_mapb = new int[num_border_elems];
    ne_get_elem_map( Exodus_in_ID, elem_mapi, elem_mapb, my_proc_id );
    for (i =0; i < new_num_internal_elems; i++) 
      new_elem_mapi[i] = elem_mapi[i];
    ne_put_elem_map( Exodus_out_ID, new_elem_mapi, elem_mapb, my_proc_id );
    delete [] elem_mapi;
    delete [] new_elem_mapi;
    delete [] elem_mapb;

    // basic communication parameters
    int * comm_node_proc_ids = new int[num_comm_procs];
    int * comm_elem_proc_ids = new int[num_elem_comm_procs];
    int * num_nodes_to_proc = new int[num_comm_procs];
    int * num_elems_to_proc = new int[num_elem_comm_procs];
    ne_get_cmap_params( Exodus_in_ID, comm_node_proc_ids, num_nodes_to_proc,
			comm_elem_proc_ids, num_elems_to_proc, my_proc_id );
    ne_put_cmap_params( Exodus_out_ID, comm_node_proc_ids, num_nodes_to_proc,
			comm_elem_proc_ids, num_elems_to_proc, my_proc_id );

    // node and element communication maps
    for( i=0 ; i<num_comm_procs ; i++ ){
      int * node_comm_ids = new int[num_nodes_to_proc[i]];
      int * comm_nodes = new int[num_nodes_to_proc[i]];
      int * elem_ids = new int[num_elems_to_proc[i]];
      int * side_ids = new int[num_elems_to_proc[i]];
      int * elem_proc_ids = new int[num_elems_to_proc[i]];

      ne_get_node_cmap( Exodus_in_ID, comm_node_proc_ids[i], comm_nodes,
			node_comm_ids, my_proc_id );
      ne_put_node_cmap( Exodus_out_ID, comm_node_proc_ids[i], comm_nodes,
			node_comm_ids, my_proc_id );

      ne_get_elem_cmap( Exodus_in_ID, comm_elem_proc_ids[i], elem_ids, 
			side_ids, elem_proc_ids, my_proc_id );
      ne_put_elem_cmap( Exodus_out_ID, comm_elem_proc_ids[i], elem_ids, 
			side_ids, elem_proc_ids, my_proc_id );

      delete [] node_comm_ids;
      delete [] comm_nodes;
      delete [] elem_ids;
      delete [] side_ids;
      delete [] elem_proc_ids;
    }
    delete [] comm_node_proc_ids;
    delete [] comm_elem_proc_ids;
    delete [] num_nodes_to_proc;
    delete [] num_elems_to_proc;
  }

  // transfer node map
  int * node_ids = new int[Number_of_Nodes];
  int * elem_ids = new int[Number_of_Elements];
  int * new_elem_ids = new int[Total_Number_of_Faces];
  ex_get_node_num_map(Exodus_in_ID, node_ids);
  ex_put_node_num_map(Exodus_out_ID, node_ids);
  ex_get_elem_num_map( Exodus_in_ID, elem_ids );
  for ( i = 0; i < Total_Number_of_Faces; i++) {
    new_elem_ids[i] = elem_ids[i];
  }
  ex_put_elem_num_map( Exodus_out_ID, new_elem_ids );
  delete [] node_ids;
  delete [] elem_ids;
  delete [] new_elem_ids;
#endif

  // cleanup

  ex_close (Exodus_in_ID);
  ex_close (Exodus_out_ID);

  delete [] x_coord;
  delete [] y_coord;
  delete [] z_coord;
  delete [] x_disp;
  delete [] y_disp;
  delete [] z_disp;
  delete [] connectivity;
  delete [] num_faces_per_block;
  delete [] nodes_per_face;
  delete [] attributes_per_face;
  delete [] face_block_index;
  for (i = 0 ; i < num_face_blocks; i++) 
    delete [] face_names[i];
  delete [] face_names;
  for (i = 0; i < num_info_records; i++) 
    delete [] info_strings[i];
  delete [] info_strings;

  return(0);
}
    
