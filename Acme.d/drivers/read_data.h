// $Id: read_data.h,v 2002.10 2004/03/11 16:07:27 mwglass Exp $

#ifndef _read_data_h_
#define _read_data_h_

#include "Contact_Defines.h"
#include "params.h"

extern "C" {
  void FORTRAN(read_data) ( int & num_load_steps,
			    int & use_config,
			    int & search_type, 
			    int & restart_flag,
			    int & start_step,
			    int & multiple_interaction_status,
			    Real* interaction_data,
			    int & normal_smoothing_status,
			    Real* smoothing_data,
			    int & compute_node_areas,
			    int & number_of_analytical_surfaces, 
			    int * as_type, 
			    Real* as_data,
			    int& num_face_blocks,
			    int* face_block_types,
			    Real* shell_block_offset,
			    Real* search_data,
			    int&  number_of_tables,
			    int*  table_ids,
			    int*  table_num_points,
			    Real* table_abscissas,
			    Real* table_ordinates,
			    char* filename,
			    int& enforcement_type,
			    Real* enforcement_data,
			    Real* enforcement_data_vars,
			    int& num_enf_models,
			    int* enf_model_types,
			    int* enf_model_ids, 
  	                    int* enf_model_sub_num,
  	                    int* enf_model_sub_ids,
 	                    char* enf_model_sub_names, 
			    Real* enf_model_data,
			    int& update_topology, 
			    int* face_ignore_list, 
			    int* num_face_ignore);
}

#endif

