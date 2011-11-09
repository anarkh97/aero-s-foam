// $Id: test_init.C,v 2002.43 2004/06/18 17:47:04 mwglass Exp $

#include "params.h"

#include "ContactSearch.h"
#include "Contact_Defines.h"
#include "output_interactions.h"
#include <stdlib.h>
#include <unistd.h>
#include <iostream.h>
#include <fstream.h>
#include "contact_assert.h"
#ifndef CONTACT_NO_EXODUS_INPUT
#include "exodusII.h"
#endif
#include <string.h>
#include "read_data.h"
#include "test_init.h"

#include <stdio.h>
#include "Contact_Communication.h"

#include "ContactTDEnforcement.h"
#include "ContactGapRemoval.h"
#include "ContactVolumeTransfer.h"
#include "DriverTopology.h"

DriverTopology* topology;

void FORTRAN(test_init) (int&    idim,
			 int&    num_load_steps,
			 int&    use_config,
			 int&    search_type, 
			 int&    restart_flag,
			 int&    start_step, 
			 int&    multiple_interaction_status,
			 double* interaction_data, 
			 int&    normal_smoothing_status,
			 double* smoothing_data, 
			 int&    compute_node_areas,
			 int&    num_as, 
			 int*    as_type,
			 double* as_data, 
			 double* search_data, 
			 int&    number_of_tables,
			 int*    table_ids,
			 int*    table_num_points,
			 double* table_abscissas,
			 double* table_ordinates,
			 int&    enforcement_type, 
			 double* enforcement_data,
			 double* enf_data_vars,
			 int&    num_enf_models, 
			 int*    enf_model_types,
			 int*    enf_model_ids, 
  	                 int*    enf_model_sub_num,
  	                 int*    enf_model_sub_ids,
 	                 char*   enf_model_sub_names,  
			 double* enf_model_data,
			 int*    num_node_ignore, 
			 int*    node_ignore_list, 
			 int*    num_face_ignore, 
			 int*    face_ignore_list,
			 int*    num_elem_ignore, 
			 int*    elem_ignore_list,
			 int&    update_topology,                       
			 int&    my_proc_id,
			 int&    num_procs,
                         int&    icomm) 
{
  int i;
#ifndef CONTACT_NO_MPI
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#else
  int comm = icomm;
#endif
#ifdef CONTACT_DEBUG
  //
  // enable ddd debugging of MPI processes
  //
  char* break_var = getenv("CONTACT_BREAK");
  int break_int=0;
  if( my_proc_id == 0 && break_var ) break_int = atoi(break_var);
  if (contact_global_sum(break_int,comm)) {
#ifdef __sun
    char buf[80];
    char hostname[80]; gethostname(hostname, sizeof(hostname));
    sprintf(buf, "Host: %s   PID: %d", hostname, getpid());
    for( int i=0 ; i<num_procs ; i++ ){
      if( i==my_proc_id )
	cout << "Processor " << my_proc_id << ":" << buf << endl;
      contact_global_sync(comm);
    }
#endif
    if (my_proc_id == 0) {
      cout << "You have requested that the code to halt with the contact"
	   << endl
	   << "environment variable CONTACT_BREAK" << endl;
      cout << "Enter any character to continue." << endl;
      cin.get();
    }
    contact_global_sync(comm);
  }
#endif

  char current_dir[1024], base_name[1024], inp_file[1024];
  // variables for configuration file data
  char cfg_root_dir[1024], cfg_sub_dir[1024];
  int cfg_num_procs(-1), cfg_num_raids(-1), cfg_raid_offset(-1),
      cfg_zero(-1);
    
  get_config( icomm, my_proc_id, current_dir, base_name, cfg_root_dir, 
	      cfg_sub_dir, cfg_num_procs, cfg_num_raids, cfg_raid_offset,
	      cfg_zero );

  if( my_proc_id == 0 ){
    // put together the full path for the data file
    strcpy(inp_file,current_dir);
    strcat(inp_file,"/");
    strcat(inp_file,base_name);
    strcat(inp_file,".inp");
  }

  //
  //  READ ASCII DATA FILE
  //
  
  // read data file on processor 0
  int    number_face_blocks;
  int    face_block_types[MAX_BLOCKS];
  double shell_block_offsets[MAX_BLOCKS];
  if ( my_proc_id == 0 ) {
    FORTRAN(read_data)( num_load_steps, use_config,
                        search_type, restart_flag,
                        start_step, multiple_interaction_status, 
                        interaction_data, normal_smoothing_status, 
                        smoothing_data, compute_node_areas, num_as, as_type, 
                        as_data, number_face_blocks, face_block_types, 
                        shell_block_offsets, search_data, number_of_tables, 
                        table_ids, table_num_points, table_abscissas,
                        table_ordinates, inp_file,
                        enforcement_type, enforcement_data,
                        enf_data_vars,
                        num_enf_models, 
                        enf_model_types,
                        enf_model_ids, 
  	                enf_model_sub_num,
  	                enf_model_sub_ids,
 	                enf_model_sub_names, 
                        enf_model_data,
                        update_topology, face_ignore_list, num_face_ignore);
  }
  char gen_file[1024];
  if ( num_procs == 1 ) {
    // create name for file to open
    strcpy(gen_file,current_dir);
    strcat(gen_file,"/");
    strcat(gen_file,base_name);
    strcat(gen_file,".gen");
#ifndef CONTACT_NO_MPI
  } else {
    contact_global_sync(comm);
    // broadcast data to all procs
    contact_broadcast(&num_load_steps,1,0,comm);
    if( num_load_steps > 1 )
      contact_broadcast(&use_config,1,0,comm);
    contact_broadcast(&search_type,1,0,comm);
    contact_broadcast(&restart_flag,1,0,comm);
    contact_broadcast(&start_step,1,0,comm);
    contact_broadcast(&update_topology,1,0,comm);
    if( restart_flag != 1 && restart_flag != 3 ) {
      // broadcast data to all procs
      contact_broadcast(&multiple_interaction_status,1,0,comm);
      if( multiple_interaction_status == 1 )
	contact_broadcast(interaction_data,1,0,comm);
      contact_broadcast(&normal_smoothing_status,1,0,comm);
      if( normal_smoothing_status == 1 )
	contact_broadcast(smoothing_data,3,0,comm);
      contact_broadcast(&compute_node_areas,1,0,comm);
      contact_broadcast(&num_as,1,0,comm);
      contact_broadcast(as_type,num_as,0,comm);
      contact_broadcast(as_data,num_as*8,0,comm);
      contact_broadcast(&number_face_blocks,1,0,comm);
      contact_broadcast(face_block_types,number_face_blocks,0,comm);
      contact_broadcast(shell_block_offsets,number_face_blocks,0,comm);
      contact_broadcast(&number_of_tables,1,0,comm);
      if( number_of_tables ){
	contact_broadcast(table_ids,number_of_tables,0,comm);
	contact_broadcast(table_num_points,number_of_tables,0,comm);
	int bsize = 0;
	for( i=0 ; i<number_of_tables ; i++ ) 
	  bsize += table_num_points[i];
	contact_broadcast(table_abscissas,bsize,0,comm);
	contact_broadcast(table_ordinates,bsize,0,comm);
      }
    } else {
      contact_broadcast(&number_face_blocks,1,0,comm);
      contact_broadcast(face_block_types,number_face_blocks,0,comm);
      contact_broadcast(shell_block_offsets,number_face_blocks,0,comm);
    }
    contact_broadcast(&enforcement_type,1,0,comm);
    if( enforcement_type==1 )
      contact_broadcast(enf_data_vars,6,0,comm);
    else if( enforcement_type == 2 )
      contact_broadcast(enf_data_vars,2,0,comm);

    contact_broadcast(&num_enf_models,1,0,comm);
    if( restart_flag != 1 && restart_flag != 3 ){
      if( num_enf_models ){
	contact_broadcast( enf_model_ids, num_enf_models, 0, comm );
	contact_broadcast( enf_model_types, num_enf_models, 0, comm );
	contact_broadcast( enf_model_sub_num, num_enf_models, 0, comm );
	contact_broadcast( enf_model_sub_ids, 
                           num_enf_models*MAX_ENF_MODEL_SUBS, 0, comm );
	contact_broadcast( enf_model_sub_names, 
                           num_enf_models*MAX_ENF_MODEL_SUBS*64, 0, comm );
	contact_broadcast( enf_model_data, 
			   num_enf_models*MAX_ENF_MODEL_DATA_VARS, 0, comm );
      }
    }
    if ( update_topology == 1 || update_topology == 3) {
      contact_broadcast(num_face_ignore,MAX_STEPS,0,comm);
      contact_broadcast(face_ignore_list,MAX_ELEMENTS*(MAX_STEPS+1),0,comm);
    }
    //
    //  READ EXODUS INPUT FILE
    //
    // from the config file, determine the directory that we will get
    // the spread mesh files from
    char exodus_dir[1024];
    char trans[81];
    strcpy(exodus_dir,cfg_root_dir);
    int on_raid = my_proc_id%cfg_num_raids+cfg_raid_offset;
    if (on_raid < 10 && cfg_zero == 1) strcat(exodus_dir,"0");
    sprintf(trans,"%d",on_raid);
    strcat(exodus_dir,trans);
    strcat(exodus_dir,"/");
    strcat(exodus_dir,cfg_sub_dir);
    strcat(exodus_dir,"/");
    
    // get name for parallel mesh file (absolute path)
    strcpy(gen_file,exodus_dir);
    strcat(gen_file,base_name);
    strcat(gen_file,".par.");
    sprintf(trans,"%d",num_procs);
    strcat(gen_file,trans);
    strcat(gen_file,".");
    PRECONDITION(num_procs < 10000);
    if (num_procs >= 1000) {
      sprintf(trans,"%04d",my_proc_id);
    }
    else if (num_procs >= 100) {
      sprintf(trans,"%03d",my_proc_id);
    }
    else if (num_procs >= 10) {
      sprintf(trans,"%02d",my_proc_id);
    }
    else {
      sprintf(trans,"%d",my_proc_id);
    }
    strcat(gen_file,trans);
#endif  
  }
  
  topology = new DriverTopology(gen_file, idim, 
                              face_block_types,
                              shell_block_offsets, 
                              my_proc_id, num_procs, comm);
#ifndef CONTACT_NO_MPI
  if ( num_procs > 1 ) {
    if( restart_flag != 1  && restart_flag != 3 ){
      // broadcast search data to all procs
      int ikeys = num_as-1+
                  topology->NumNodeBlocks()+
                  topology->NumFaceBlocks()+
                  topology->NumElemBlocks();  
      contact_broadcast(search_data,ContactSearch::NSIZSD*ikeys*ikeys,0,comm);
      switch( enforcement_type ){
      case(1):{
	contact_broadcast(enforcement_data,
			  ContactTDEnforcement::NSIZED*ikeys*ikeys,0,comm);
	break;
      }
      case(2):{
	contact_broadcast(enforcement_data,
			  ContactGapRemoval::NSIZED*ikeys*ikeys,0,comm);
	break;
      }
      case(4):{
	contact_broadcast(enforcement_data,
			  ContactVolumeTransfer::NSIZED*ikeys*ikeys,0,comm);
	break;
      }
      }
    }
  }
#endif 
}

void FORTRAN(start_topology)(int & num_node_ignore,
                             int * node_ignore_list,
                             int & num_face_ignore,
                             int * face_ignore_list,
                             int & num_elem_ignore,
                             int * elem_ignore_list,
                             
                             int & number_node_blocks,
                             int * node_block_types,
                             int * nodes_per_block,
                             int * node_eids,
                             int * node_ids,
                             
                             int & number_face_blocks,
                             int * face_block_types,
                             int * faces_per_block,
                             int * face_ids,
                             int * face_conn,
                             
                             int & number_element_blocks,
                             int * element_block_types,
                             int * elements_per_block,
                             int * element_ids,
                             int * element_conn,
                                                         
                             double* position,
                             double* nodal_vars,
                             int&    num_nodal_vars,
                             double* element_vars,
                             int&    num_element_vars,
                             double* attributes,
                             int*    num_attributes,
                             double* node_radius,
                             double* shell_thickness,
                             double* shell_offset,
                             
                             int & num_comm_partners,
                             int * comm_proc_ids,
                             int * number_nodes_to_partner,
                             int * comm_node,
                             
                             int & my_proc_id,
                             int & num_procs)
{
  topology->GetInitialTopology(num_node_ignore,
                               node_ignore_list,
                               num_face_ignore,
                               face_ignore_list,
                               num_elem_ignore,
                               elem_ignore_list,
                               
                               number_node_blocks,
                               node_block_types,
                               nodes_per_block,
                               node_eids,
                               node_ids,
                               
                               number_face_blocks,
                               face_block_types,
                               faces_per_block,
                               face_ids,
                               face_conn,
                               
                               number_element_blocks,
                               element_block_types,
                               elements_per_block,
                               element_ids,
                               element_conn,
                                                     
                               position,
                               nodal_vars,
                               num_nodal_vars,
                               element_vars,
                               num_element_vars,
                               attributes,
                               num_attributes,
                               node_radius,
                               shell_thickness,
                               shell_offset,
                               
                               num_comm_partners,
                               comm_proc_ids,
                               number_nodes_to_partner,
                               comm_node,
                               
                               my_proc_id,
                               num_procs);
}

void FORTRAN(get_topology_birthdeath)(
                             int& num_node_ignore,
                             int* node_ignore_list,
                             int& num_face_ignore,
                             int* face_ignore_list,
                             int& num_elem_ignore,
                             int* elem_ignore_list,
                         
                             int* num_node_deaths_per_block, 
	                     int* node_deaths_global_ids,
                             int* num_face_deaths_per_block, 
	                     int* face_deaths_global_ids,
                             int* num_element_deaths_per_block, 
	                     int* element_deaths_global_ids,
                             
                             int* num_node_births_per_block, 
	                     int* node_births_exodus_ids,
	                     int* node_births_global_ids,
	                     int* number_face_births_per_block, 
	                     int* face_births_global_ids,
	                     int* face_births_connectivity,
	                     int* number_element_births_per_block, 
	                     int* element_births_global_ids,
	                     int* element_births_connectivity,
                             
                             int & my_proc_id,
                             int & num_procs)
{
  topology->GetModsForBirthDeath(num_node_ignore,
                                 node_ignore_list,
                                 num_face_ignore,
                                 face_ignore_list,
                                 num_elem_ignore,
                                 elem_ignore_list,
                               
                                 num_node_deaths_per_block, 
	                         node_deaths_global_ids,
                                 num_face_deaths_per_block, 
	                         face_deaths_global_ids,
                                 num_element_deaths_per_block, 
	                         element_deaths_global_ids,
                                
                                 num_node_births_per_block, 
	                         node_births_exodus_ids,
	                         node_births_global_ids,
	                         number_face_births_per_block, 
	                         face_births_global_ids,
	                         face_births_connectivity,
	                         number_element_births_per_block, 
	                         element_births_global_ids,
	                         element_births_connectivity,
                                  
                                 my_proc_id,
                                 num_procs);
}

void FORTRAN(get_topology_comm_plan)(int& num_comm_partners,
                                     int* comm_proc_ids,
                                     int* number_nodes_to_partner,
                                     int* comm_node,
                                     int& my_proc_id,
                                     int& num_procs)
{
  topology->GetCommPlan(num_comm_partners,
                        comm_proc_ids,
                        number_nodes_to_partner,
                        comm_node,
                        my_proc_id,
                        num_procs);
}

void FORTRAN(get_topology_host_ids)(int& number_of_nodes,
                                    int* node_eids,
                                    int* node_ids,
                                    int& number_of_faces,
                                    int* face_ids,
                                    int& number_of_elems,
                                    int* elem_ids)
{
  topology->GetHostIDs(number_of_nodes,
                       node_eids,
                       node_ids,
                       number_of_faces,
                       face_ids,
                       number_of_elems,
                       elem_ids);
}

void FORTRAN(get_topology_variables)(double* position,
                                     double* nodal_vars,
                                     int&    num_nodal_vars,
                                     double* element_vars,
                                     int&    num_element_vars,
                                     double* attributes,
                                     int*    num_attributes,
                                     double* node_radius,
                                     double* shell_thickness,
                                     double* shell_offset)
{
  topology->GetVariables(position,
                         nodal_vars,
                         num_nodal_vars,
                         element_vars,
                         num_element_vars,
                         attributes,
                         num_attributes,
                         node_radius,
                         shell_thickness,
                         shell_offset);
}

void FORTRAN(get_topology_dlb)(int& num_node_exports,
                               int* node_export_id_list,
                               int* node_export_pids,
		               int& num_face_exports,
                               int* face_export_id_list,
                               int* face_export_pids,
		               int& num_elem_exports,
                               int* elem_export_id_list,
                               int* elem_export_pids)
{
  topology->GetModsForDLB(num_node_exports,
                          node_export_id_list,
                          node_export_pids,
		          num_face_exports,
                          face_export_id_list,
                          face_export_pids,
		          num_elem_exports,
                          elem_export_id_list,
                          elem_export_pids);
}
                                     
extern "C" void FORTRAN(test_close)(int& exodus_id)
{
#ifndef CONTACT_NO_EXODUS_INPUT
  ex_close( exodus_id );
#endif
}

extern "C" void FORTRAN(create_plot_file)( int& icomm, int& my_proc_id,
					   int& plot_step, int& exodus_out )
{
  char current_dir[1024], base_name[1024], cfg_root_dir[1024];
  char cfg_sub_dir[1024];
  int cfg_num_procs, cfg_num_raids, cfg_raid_offset, cfg_zero;

#ifndef CONTACT_NO_MPI
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#else
  int comm = icomm;
#endif
  int num_procs = contact_number_of_processors( comm );
  
  get_config( icomm, my_proc_id, current_dir, base_name, cfg_root_dir,
	      cfg_sub_dir, cfg_num_procs, cfg_num_raids,
	      cfg_raid_offset, cfg_zero );

  int iows = 8;
  int compws = 8;
  char out_file[1024],exodus_dir[1024];
  char trans[81];
  if( num_procs == 1 ){       
    strcpy(out_file,current_dir);
    strcat(out_file,"/");
    strcat(out_file,base_name);
    strcat(out_file,".exo");
    if( plot_step > 1 ){
      sprintf(trans,"_%d",plot_step);
      strcat(out_file,trans);
    }
  } else {
    strcpy(exodus_dir,cfg_root_dir);
    int on_raid = my_proc_id%cfg_num_raids+cfg_raid_offset;
    if (on_raid < 10 && cfg_zero == 1) strcat(exodus_dir,"0");
    sprintf(trans,"%d",on_raid);
    strcat(exodus_dir,trans);
    strcat(exodus_dir,"/");
    strcat(exodus_dir,cfg_sub_dir);
    strcat(exodus_dir,"/");
    strcpy(out_file,exodus_dir);
    strcat(out_file,base_name);
    strcat(out_file,".exo");
    if( plot_step > 1 ){
      sprintf(trans,"_%d",plot_step);
      strcat(out_file,trans);
    }
    strcat(out_file,".");
    sprintf(trans,"%d",num_procs);
    strcat(out_file,trans);
    strcat(out_file,".");
    PRECONDITION(num_procs < 10000);
    if (num_procs >= 1000) {
      sprintf(trans,"%04d",my_proc_id);
    }
    else if (num_procs >= 100) {
      sprintf(trans,"%03d",my_proc_id);
    }
    else if (num_procs >= 10) {
      sprintf(trans,"%02d",my_proc_id);
    }
    else {
      sprintf(trans,"%d",my_proc_id);
    }
    strcat(out_file,trans);
  }
#ifndef CONTACT_NO_EXODUS_OUTPUT
  exodus_out = ex_create( out_file, EX_CLOBBER, &compws, &iows );
#else
  exodus_out = 0;
#endif
}

extern "C" void get_config(int icomm, int my_proc_id, char current_dir[1024],
			   char base_name[1024], char cfg_root_dir[1024],
			   char cfg_sub_dir[1024], int& cfg_num_procs,
			   int& cfg_num_raids,
			   int& cfg_raid_offset, int& cfg_zero )
{
    // get current working directory and base name for this test
  
#ifndef CONTACT_NO_MPI
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#else
  int comm = icomm;
#endif

  int len_current_dir, len_root_dir, len_sub_dir, len_base_name;
  int num_procs = contact_number_of_processors( comm );
  
  if ( my_proc_id == 0 ) {
    
    // DEBUG -- print out current working directory
    // get current working directory
    char * dir_pointer = getenv("CONTACT_CWD");
    if ( NULL == dir_pointer ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      cerr << "FATAL ERROR: the environment variable CONTACT_CWD" << endl
	   << "   must be set to the test directory." << endl;
#endif
      exit(1);
    }
    strcpy(current_dir,dir_pointer);
    len_current_dir = strlen(current_dir)+1;
    
    // strip out the base name from the directory -- find last "/" and
    // then copy remainer of string
    int char_index = len_current_dir; 
    do { 
      char_index--;
    } while (char_index >= 0 && current_dir[char_index] != '/');
    if (char_index < 0) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      cerr << "ERROR: problem with parsing out base name" << endl;
#endif
      exit(1);
    }
    else {
      for ( int i = 0; i<=len_current_dir-char_index-1; i++)
	base_name[i] = current_dir[i+char_index+1];
    }
  }

  // have processor zero get all the data needed to send to the other 
  // processors
  if ( my_proc_id == 0 ) {
    
    // get current working directory
    char * dir_pointer = getenv("CONTACT_CWD");
    strcpy(current_dir,dir_pointer);
    len_current_dir = strlen(current_dir)+1;
    
    // strip out the base name from the directory -- find last "/" and
    // then copy remainer of string
    int char_index = len_current_dir; 
    do { 
      char_index--;
    } while (char_index >= 0 && current_dir[char_index] != '/');
    if (char_index < 0) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      cerr << "ERROR: problem with parsing out base name" << endl;
#endif
      exit(1);
    }
    else {
      for ( int i = 0; i<=len_current_dir-char_index-1; i++)
	base_name[i] = current_dir[i+char_index+1];
    }
    len_base_name = strlen(base_name)+1;
    
    if( num_procs > 1 ){
      // find and open configuration file
      ifstream config_file;
      char cfg_file[1024];
      strcpy(cfg_file,current_dir);
      strcat(cfg_file,"/");
      strcat(cfg_file,base_name);
      strcat(cfg_file,".cfg");
      
      // process configuration file 
      config_file.open(cfg_file);
      if ( ! config_file ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	cerr << "ERROR: Cannot open " << base_name << ".cfg." << endl;
#endif
	exit(1);
      }
      
      char word[1024];
      for (;;) {
	if (config_file.fail()) break;
	
	// get first word -- skip if not "set"
	config_file >> word;
	// get thing to set
	config_file >> word;
	if (strcmp(word,"num_procs") == 0) {
	  config_file >> word;
	  config_file >> cfg_num_procs;
	}
	else if (strcmp(word,"root_dir") == 0) {
	  config_file >> word;
	  config_file >> cfg_root_dir;
	}
	else if (strcmp(word,"sub_dir") == 0) {
	  config_file >> word;
	  config_file >> cfg_sub_dir;
	}
	else if (strcmp(word,"num_raids") == 0) {
	  config_file >> word;
	  config_file >> cfg_num_raids;
	}
	else if (strcmp(word,"raid_offset") == 0) {
	  config_file >> word;
	  config_file >> cfg_raid_offset;
	}
	else if (strcmp(word,"zero") == 0) {
	  config_file >> word;
	  int test;
	  if ((config_file >> test).fail()) {
	    cfg_zero = 0;
	    config_file.clear();
	  }
	  else {
	    cfg_zero = 1;
	  }
	}
	else {
	  // We don't want this line (e.g., MACHINE)
	  while ( ( config_file.get())!='\n' && !config_file.fail() );
	}	
      }
      
      // now we have read config file, check some values    
      if (cfg_num_procs == -1 || cfg_num_raids == -1 || 
	  cfg_raid_offset == -1 || cfg_zero == -1) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	cerr << "ERROR in configuration file." << endl;
#endif
#ifndef CONTACT_NO_MPI
	MPI_Abort(comm, 1);
#endif
      }
      
      if (cfg_num_procs != num_procs ) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	cerr << "ERROR: decomposition is for " << cfg_num_procs 
	     << " processors, but " << num_procs 
	     << " processors are executing." << endl;
#endif
#ifndef CONTACT_NO_MPI
	MPI_Abort(comm, 1);
#endif
      }
      
      len_sub_dir = strlen(cfg_sub_dir)+1;
      len_root_dir = strlen(cfg_root_dir)+1;
      
    }
  }

  // now broadcast the info read in the configuration file to all processors
  contact_broadcast(&len_current_dir,1,0,comm);
  contact_broadcast(current_dir,len_current_dir,0,comm);
  contact_broadcast(&len_base_name,1,0,comm);
  contact_broadcast(base_name,len_base_name,0,comm);
  if( num_procs > 1 ){
    contact_broadcast(&len_sub_dir,1,0,comm);
    contact_broadcast(cfg_sub_dir,len_sub_dir,0,comm);
    contact_broadcast(&len_root_dir,1,0,comm);
    contact_broadcast(cfg_root_dir,len_root_dir,0,comm);
    contact_broadcast(&cfg_num_procs,1,0,comm);
    contact_broadcast(&cfg_num_raids,1,0,comm);
    contact_broadcast(&cfg_raid_offset,1,0,comm);
    contact_broadcast(&cfg_zero,1,0,comm);
  }
}
