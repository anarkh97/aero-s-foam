// $Id: test_init.h,v 2002.20 2004/06/18 17:47:04 mwglass Exp $

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#endif

#include "Contact_Defines.h"

#ifndef _test_init_
#define _test_init_

extern "C" {
  void FORTRAN(test_init)( int&  idim, 
			   int&  num_load_steps,
			   int&  use_config,
			   int&  search_type, 
			   int&  restart_flag,
			   int&  start_step,
			   int&  multiple_interaction_status,
			   Real* interaction_data,
			   int&  normal_smoothing_status,
			   Real* smoothing_data,
                           int&  compute_node_areas,
			   int&  num_as, 
			   int*  as_type,
			   Real* as_data, 
			   Real* search_data, 
			   int&  number_of_tables,
			   int*  table_ids,
			   int*  table_num_points,
			   Real* table_abscissas,
			   Real* table_ordinates,
			   int&  enforcement_type,
			   Real* enforcement_data,
			   Real* enforcement_data_vars,
			   int&  num_enf_models,
			   int*  enf_model_types,
			   int*  enf_model_ids,
  	                   int*  enf_model_sub_num,
  	                   int*  enf_model_sub_ids,
 	                   char* enf_model_sub_names, 
			   Real* enf_model_data,
			   int*  num_node_ignore, 
			   int*  node_ignore_list, 
			   int*  num_face_ignore, 
			   int*  face_ignore_list, 
			   int*  num_elem_ignore, 
			   int*  elem_ignore_list,
			   int&  update_topology,
			   int&  my_proc_id,
			   int&  num_procs,
                           int&  icomm );

  void FORTRAN(start_topology)(int & num_node_ignore,
                               int * node_ignore_list,
                               int & num_face_ignore,
                               int * face_ignore_list,
                               int & num_elem_ignore,
                               int * elem_ignore_list,
                               
                               int & number_node_blocks,
                               int * node_block_types,
                               int * mod_nodes_per_block,
                               int * mode_node_eids,
                               int * mod_node_ids,
                               
                               int & number_face_blocks,
                               int * face_block_types,
                               int * mod_faces_per_block,
                               int * mod_face_ids,
                               int * mod_face_conn,
                               
                               int & number_element_blocks,
                               int * element_block_types,
                               int * mod_elements_per_block,
                               int * mod_element_ids,
                               int * mod_element_conn,
                               
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
                               int & num_procs);

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
                             int & num_procs);

  void FORTRAN(get_topology_comm_plan)(int& num_comm_partners,
                                     int* comm_proc_ids,
                                     int* number_nodes_to_partner,
                                     int* comm_node,
                                     int& my_proc_id,
                                     int& num_procs);

  void FORTRAN(get_topology_variables)(double* position,
                                       double* nodal_vars,
                                       int&    num_nodal_vars,
                                       double* element_vars,
                                       int&    num_element_vars,
                                       double* attributes,
                                       int*    num_attributes,
                                       double* node_radius,
                                       double* shell_thickness,
                                       double* shell_offset);

  void FORTRAN(get_topology_host_ids)(int& number_of_nodes,
                                      int* node_eids,
                                      int* node_ids,
                                      int& number_of_faces,
                                      int* face_ids,
                                      int& number_of_elems,
                                      int* elem_ids);
                                    
  void FORTRAN(get_topology_dlb)(int& num_node_exports,
                                 int* node_export_id_list,
                                 int* node_export_pids,
		                 int& num_face_exports,
                                 int* face_export_id_list,
                                 int* face_export_pids,
		                 int& num_elem_exports,
                                 int* elem_export_id_list,
                                 int* elem_export_pids);
                               
  void FORTRAN(test_close)(int& exodus_id);

  void FORTRAN(create_plot_file)(int& icomm, int& my_proc_id,
				 int& plot_step, int& exodus_out);
  
  void get_config( int icomm, int my_proc_id, char current_dir[1024],
		   char base_name[1024], char chg_root_dir[1024],
		   char cfg_sub_dir[1024], int& cfg_num_procs, 
		   int& cfg_num_raids,
		   int& cfg_raid_offset, int& cfg_zero );
}
#endif
