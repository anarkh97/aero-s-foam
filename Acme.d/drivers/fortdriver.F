C $Id: fortdriver.F,v 2002.96 2004/06/18 17:47:04 mwglass Exp $
      program test_driver

      implicit real*8 (a-h,o-z)

#include "params.h"
      parameter (NSIZSD = 3)
      
      
      integer*4 num_node_deaths_per_block(MAX_BLOCKS)
      integer*4 node_deaths_global_ids(2*MAX_NODES)
      integer*4 num_face_deaths_per_block(MAX_BLOCKS)
      integer*4 face_deaths_global_ids(2*MAX_ELEMENTS)
      integer*4 num_element_deaths_per_block(MAX_BLOCKS)
      integer*4 element_deaths_global_ids(2*MAX_ELEMENTS)
      
      integer*4 num_node_births_per_block(MAX_BLOCKS)
      integer*4 nbirths_exodus_ids(MAX_NODES)
      integer*4 nbirths_global_ids(2*MAX_NODES)
      integer*4 number_face_births_per_block(MAX_BLOCKS)
      integer*4 fbirths_global_ids(2*MAX_ELEMENTS)
      integer*4 number_element_births_per_block(MAX_BLOCKS)
      integer*4 ebirths_global_ids(2*MAX_ELEMENTS)
      integer*4 node_export_id_list(2*MAX_NODES)
      integer*4 node_export_pid(MAX_NODES)
      integer*4 face_export_id_list(2*MAX_ELEMENTS)
      integer*4 face_export_pid(2*MAX_ELEMENTS)
      integer*4 element_export_id_list(2*MAX_ELEMENTS)
      integer*4 element_export_pid(2*MAX_ELEMENTS)
      
      
      integer*4 num_elem_ignore(MAX_STEPS)
      integer*4 elem_ignore_list((MAX_STEPS+1)*MAX_ELEMENTS)
      integer*4 num_face_ignore(MAX_STEPS)
      integer*4 face_ignore_list((MAX_STEPS+1)*MAX_ELEMENTS)
      integer*4 num_node_ignore(MAX_STEPS)
      integer*4 node_ignore_list((MAX_STEPS+1)*MAX_NODES)
      integer*4 node_eids(MAX_NODES)
      integer*4 node_ids(2*MAX_NODES)
      integer*4 face_ids(2*MAX_ELEMENTS)
      integer*4 elem_ids(2*MAX_ELEMENTS)
      dimension ifconn(2*MAX_NODES_PER_ELEM*MAX_ELEMENTS)
      dimension ieconn(2*MAX_NODES_PER_ELEM*MAX_ELEMENTS)
      dimension node_block_ids(MAX_BLOCKS)
      dimension node_block_types(MAX_BLOCKS)
      dimension nodes_per_block(MAX_BLOCKS)
      dimension iface_block_types(MAX_BLOCKS)
      dimension ifaces_per_block(MAX_BLOCKS)
      dimension ielem_block_types(MAX_BLOCKS)
      dimension ielems_per_block(MAX_BLOCKS)
      dimension pos_init(3*MAX_NODES)
      dimension pos1(3*MAX_NODES),pos2(3*MAX_NODES)
      dimension vars_nodal(MAX_NODES*MAX_NODAL_VARS)
      dimension vars_element(MAX_ELEMENTS*MAX_ELEMENT_VARS)
      dimension attributes(MAX_ELEMENTS*MAX_ATTRIBUTES)
      dimension num_attributes(MAX_BLOCKS)
      dimension radius(MAX_NODES)
      dimension shell_thickness(MAX_ELEMENTS)
      dimension shell_offset(MAX_ELEMENTS)
      dimension as_data(8*MAX_ANALYTIC_SURFACES)
      integer*4 as_type(MAX_ANALYTIC_SURFACES)
      dimension node_blk_ids(MAX_INTERACTIONS)
      dimension node_indx(MAX_INTERACTIONS)
      dimension node_entity_keys(MAX_INTERACTIONS)
      dimension iface_blk_ids(MAX_INTERACTIONS)
      dimension iface_indx(MAX_INTERACTIONS)
      dimension iface_proc(MAX_INTERACTIONS)
      dimension data(7*MAX_INTERACTIONS)
      dimension search_data(NSIZSD*MAX_BLOCKS*MAX_BLOCKS)
      integer*4 number_of_tables
      integer*4 table_ids(MAX_TABLES)
      integer*4 table_num_points(MAX_TABLES)
      dimension table_abscissas(MAX_TABLES*MAX_TABLE_POINTS)
      dimension table_ordinates(MAX_TABLES*MAX_TABLE_POINTS)

      dimension icomm_proc_id(MAX_PROCS) 
      dimension number_nodes_to_partner(MAX_PROCS) 
      dimension icomm_nodes(MAX_NODES*MAX_PROCS) 

      dimension r_interaction_data(2)
      dimension smoothing_data(3)
      dimension restart_buffer(MAX_RESTART_SIZE)
      character*1 message(80)
      character*80 message2
      integer*4 CURRENT_CONFIGURATION,PREDICTED_CONFIGURATION
      parameter( CURRENT_CONFIGURATION=1,PREDICTED_CONFIGURATION=2 )
      integer*4 enforcement_type
      dimension enforcement_data(5*MAX_BLOCKS*MAX_BLOCKS)
      dimension enforcement_vars(3*MAX_NODES),enf_data_vars(6)
      dimension results(3*MAX_NODES)
      dimension enf_array_1(MAX_NODAL_VARS*MAX_NODES)
      dimension enf_array_2(MAX_NODAL_VARS*MAX_NODES)
      dimension enf_array_3(MAX_NODES)
      dimension enf_array_4(MAX_ELEMENT_VARS*MAX_ELEMENTS)
      dimension enf_array_5(MAX_ELEMENT_VARS*MAX_ELEMENTS)
      dimension enf_array_6(MAX_ELEMENTS)
c arrays for mpc test harness
      dimension impc_array_1(3*MAX_NODES)
      dimension impc_array_2(3*MAX_NODES)
      dimension impc_array_3(3*MAX_NODES)
      dimension impc_array_4(3*MAX_NODES)
      dimension impc_array_5(3*MAX_NODES)
      dimension impc_array_6(8*3*MAX_NODES)
      dimension impc_array_7(8*3*MAX_NODES)
      dimension rmpc_array_1(8*3*MAX_NODES)

      integer*4 num_enf_models
      integer*4 enf_model_ids(MAX_ENFORCEMENT_MODELS)
      integer*4 enf_model_types(MAX_ENFORCEMENT_MODELS)
      integer*4 enf_model_sub_num(MAX_ENFORCEMENT_MODELS)
      integer*4 enf_model_sub_ids
     *   (MAX_ENFORCEMENT_MODELS*MAX_ENF_MODEL_SUBS)
      character*64 enf_model_sub_names
     *   (MAX_ENFORCEMENT_MODELS*MAX_ENF_MODEL_SUBS)

#define HARDCODE_USER_SUBS
#ifdef HARDCODE_USER_SUBS
      external user_sprngweld_init_model
      external user_sprngweld_init_tstep
      external user_sprngweld_init_sdata
      external user_sprngweld_is_active
      external user_sprngweld_itype
      external user_sprngweld_limit_force
#else
      external user_initialize_model
      external user_initialize_time_step
      external user_init_node_state_data
      external user_is_active
      external user_interaction_type
      external user_limit_force
#endif


      dimension 
     *  enf_model_data(MAX_ENF_MODEL_DATA_VARS*MAX_ENFORCEMENT_MODELS)
      dimension r_emod_data(MAX_ENF_MODEL_DATA_VARS)
      integer*4 i_emod_data(MAX_ENF_MODEL_DATA_VARS)
#ifndef CONTACT_NO_MPI
#define MPI_COMPILE 1
#include "mpif.h"
      integer ierr
      call MPI_INIT(ierr)
      if ( ierr .ne. 0) then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
         write (*,*) "ERROR in starting MPI."
#endif
         stop
      endif
      icomm = MPI_COMM_WORLD
      call MPI_COMM_SIZE(icomm, num_procs, ierr)
      call MPI_COMM_RANK(icomm, my_proc_id, ierr)
      call MPI_BARRIER(icomm,ierr)
#else
#define MPI_COMPILE 0
      num_procs = 1
      my_proc_id = 0
      icomm = 0
#endif
      if( my_proc_id .eq. 0 )then
        do i=1,80
          message2(i:i) = ' '
        enddo
        call acme_version( message2 )
        write(*,*) message2(1:40)
        do i=1,80
          message2(i:i) = ' '
        enddo
        call acme_versiondate( message2 )
        write(*,*) message2(1:40)
      endif
      
      num_enf_type_1_vars = 3

      do i=1,MAX_BLOCKS
        node_block_ids(i) = i
      end do
      do i=1,MAX_BLOCKS
        node_block_types(i) = 1
      end do
      number_nodes = 0
      nodes_per_block(1) = -1
      do i=1,MAX_STEPS
         num_node_ignore(i) = 0
         num_face_ignore(i) = 0
         num_elem_ignore(i) = 0
      enddo  

      do i=1,MAX_ENF_MODEL_DATA_VARS*MAX_ENFORCEMENT_MODELS
        enf_model_data(i) = 0.0
      enddo

      call acme_mpi_compatibility(MPI_COMPILE, ierror)
      if( ierror .ne. 0 ) then
        message2 = 'Incompatible MPI compilation'
        write(*,*) message2(1:40)
#ifndef CONTACT_NO_MPI
        call MPI_FINALIZE(ierr)
#endif
        stop
      endif
#ifndef HARDCODE_USER_SUBS
      call register_user_subs()  
#endif

      num_nodal_vars = MAX_NODAL_VARS
      num_element_vars = MAX_ELEMENT_VARS
      istart_step = 1
      num_comm_partners = 0
      itop_update = 0
      call test_init(
     &     idim, 
     &     num_load_steps,
     &     iuse_config, 
     &     isrch_type, 
     &     irestart_flag, 
     &     istart_step, 
     &     imi_status, 
     &     r_interaction_data,
     &     ins_status, 
     &     smoothing_data, 
     &     icompute_node_areas,
     &     num_as, 
     &     as_type, 
     &     as_data, 
     &     search_data, 
     &     number_of_tables,
     &     table_ids,
     &     table_num_points,
     &     table_abscissas,
     &     table_ordinates,
     &     enforcement_type, 
     &     enforcement_data, 
     &     enf_data_vars,
     &     num_enf_models, 
     &     enf_model_types, 
     &     enf_model_ids,
     &     enf_model_sub_num,
     &     enf_model_sub_ids,
     &     enf_model_sub_names,
     &     enf_model_data, 
     &     num_node_ignore,
     &     node_ignore_list,
     &     num_face_ignore,
     &     face_ignore_list,
     &     num_elem_ignore, 
     &     elem_ignore_list, 
     &     itop_update,
     &     my_proc_id, 
     &     num_procs,
     &     icomm)

C modify the topology for any ignored faces
      iflag_update = -1
      if (itop_update .ne. 0) iflag_update = 0
      call start_topology(
     &     num_node_ignore,
     &     node_ignore_list,
     &     num_face_ignore,
     &     face_ignore_list,
     &     num_elem_ignore,
     &     elem_ignore_list,
     
     &     num_node_blocks,
     &     node_block_types,
     &     nodes_per_block,
     &     node_eids,
     &     node_ids,
     
     &     num_face_blocks,
     &     iface_block_types,
     &     ifaces_per_block,
     &     face_ids,
     &     ifconn,
     
     &     num_elem_blocks,
     &     ielem_block_types,
     &     ielems_per_block,
     &     elem_ids,
     &     ieconn,
     
     &     pos_init, 
     &     vars_nodal,
     &     num_nodal_vars,
     &     vars_element,
     &     num_element_vars,
     &     attributes,
     &     num_attributes,
     &     radius,
     &     shell_thickness,
     &     shell_offset,
     
     &     num_comm_partners,
     &     icomm_proc_id,
     &     number_nodes_to_partner,
     &     icomm_nodes,
     
     &     my_proc_id,
     &     num_procs)
     
      num_nodes_tot = 0
      do i = 1, num_node_blocks
         num_nodes_tot = num_nodes_tot + nodes_per_block(i)
      enddo
      num_faces_tot = 0
      do i = 1, num_face_blocks
         num_faces_tot = num_faces_tot + ifaces_per_block(i)
      enddo
      num_elems_tot = 0
      do i = 1, num_elem_blocks
         num_elems_tot = num_elems_tot + ielems_per_block(i)
      enddo
      
      if( irestart_flag.ne.1 .and. irestart_flag.ne.3 )then
C Construct the faceted portion of the topology
        num_states = 1
        dummy = 0
        call build_search(idim,
     *                    num_states,
     *                    num_as,
     *                    num_node_blocks,
     *                    node_block_types,
     *                    nodes_per_block,
     *                    node_eids,
     *                    node_ids,
     *                    pos_init,
     *                    num_face_blocks,
     *                    iface_block_types,
     *                    ifaces_per_block,
     *                    face_ids,
     *                    ifconn,
     *                    shell_offset,
     *                    num_elem_blocks,
     *                    ielem_block_types,
     *                    ielems_per_block,
     *                    elem_ids,
     *                    ieconn,
     *                    num_comm_partners,
     *                    icomm_proc_id,
     *                    number_nodes_to_partner,
     *                    icomm_nodes,
     *                    icomm,
     *                    ierror)
        if( ierror .ne. 0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
          write(*,*) "Error in build_search: Error Code = ",ierror
#endif
	  call number_of_search_errors( num_errors )
          do ii=1,num_errors
             call get_search_error_message( ii,message )
             call print_message( message )
          end do
          call cleanup_search()
#ifndef CONTACT_NO_MPI
          call MPI_FINALIZE(ierr)
#endif
          stop
        endif

C Add any analytic surfaces
        do i=1,num_as
          call add_analytic_surface( as_type(i),as_data(8*(i-1)+1),
     *                               ierror )
          if( ierror.ne.0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
            write(*,*) "Error in Add_Analytic_Surface: Error Code = ",
     *         ierror
#endif
	    call number_of_search_errors( num_errors )
            do ii=1,num_errors
               call get_search_error_message( ii,message )
               call print_message( message )
            end do
            call cleanup_search()
#ifndef CONTACT_NO_MPI
            call MPI_FINALIZE(ierr)
#endif
            stop
          endif
        end do
C Set the multiple interactions option
        call set_search_option( 0,imi_status,r_interaction_data,ierror )
        if( ierror.ne.0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
          write(*,*) "Error in set_search_option: Error Code = ",
     *       ierror
#endif
	  call number_of_search_errors( num_errors )
          do ii=1,num_errors
             call get_search_error_message( ii,message )
             call print_message( message )
          end do
          call cleanup_search()
#ifndef CONTACT_NO_MPI
          call MPI_FINALIZE(ierr)
#endif
          stop
        endif
C Set the normal smoothing option
        call set_search_option( 1, ins_status, smoothing_data, ierror )
        if( ierror.ne.0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
          write(*,*) "Error in set_search_option: Error Code = ",
     *     ierror
#endif
	  call number_of_search_errors( num_errors )
          do ii=1,num_errors
             call get_search_error_message( ii,message )
             call print_message( message )
          end do
          call cleanup_search()
#ifndef CONTACT_NO_MPI
          call MPI_FINALIZE(ierr)
#endif
          stop
        endif
C Set the compute node areas option
        call set_search_option( 2, icompute_node_areas, smoothing_data,
     *                          ierror )
        if( ierror.ne.0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
          write(*,*) "Error in set_search_option: Error Code = ",
     *     ierror
#endif
	  call number_of_search_errors( num_errors )
          do ii=1,num_errors
             call get_search_error_message( ii,message )
             call print_message( message )
          end do
          call cleanup_search()
#ifndef CONTACT_NO_MPI
          call MPI_FINALIZE(ierr)
#endif
          stop
        endif
C Set the search data
        number_entity_keys = num_as +
     *                       num_node_blocks +
     *                       num_face_blocks +
     *                       num_elem_blocks
        if ((num_face_blocks+num_elem_blocks).gt.0) then
          number_entity_keys = number_entity_keys-1
        endif
        call check_search_data_size( 3,number_entity_keys,ierror )
        if( ierror.ne.0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
         write(*,*) "Error in check_search_data_size: Error Code = ",
     *       ierror
#endif
	  call number_of_search_errors( num_errors )
          do ii=1,num_errors
             call get_search_error_message( ii,message )
             call print_message( message )
          end do
          call cleanup_search()
#ifndef CONTACT_NO_MPI
          call MPI_FINALIZE(ierr)
#endif
          stop
        endif
        call set_search_data( search_data )

C add the tables
        itable_data_offset = 1
        do i=1,number_of_tables
          call add_table( table_ids(i), table_num_points(i),
     *      table_abscissas(itable_data_offset),
     *      table_ordinates(itable_data_offset),ierror )
          itable_data_offset = itable_data_offset + table_num_points(i)
        enddo

C initialize the enforcement if required
        if( enforcement_type .eq. 1 )then
          call build_td_enforcement( enforcement_data, ierror )
        else if( enforcement_type .eq. 2 )then
          call build_gap_removal( enforcement_data, ierror )
        else if( enforcement_type .eq. 3 )then
          call build_tied_kin( enforcement_data, ierror )
        else if( enforcement_type .eq. 4 )then
          call build_vol_tran( enforcement_data, ierror )
        else if( enforcement_type .eq. 5 )then
          call build_mpc_eqns( enforcement_data, ierror )
        end if
        
C add the TD_ENF_MODEL friction models if needed
        do i=1,num_enf_models
          if( enforcement_type .eq. 1 )then
            if( enf_model_types(i) .eq. 14 )then
              i_emod_data(1) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+1)
              i_emod_data(2) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+2)
              i_emod_data(3) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+3)              
              i_emod_data(4) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+4)              
              do ii=1,i_emod_data(1)
                   i_emod_data(ii+4) 
     &              = enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+ii+4)
              enddo
              do ii=1,i_emod_data(2)
                   r_emod_data(ii) 
     &              = enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+ii+4 
     &                               + i_emod_data(1))
              enddo
              call td_add_enf_model( enf_model_types(i), 
     &          enf_model_ids(i),i_emod_data(1), 
     &          r_emod_data(1), ierror )
              do ii=1,enf_model_sub_num(i)
                j = enf_model_sub_ids((i-1)*MAX_ENF_MODEL_SUBS+ii)
                jj = ((i-1)*MAX_ENF_MODEL_SUBS+(ii-1)+1)
                call str_len(enf_model_sub_names(jj),isize)
c               write(*,*) "getting : ", enf_model_sub_names(jj),
c    &                ", index: ", jj, ", str. len.", isize
#ifdef HARDCODE_USER_SUBS
c THIS ONLY WORKS FOR user_springweld
                if (j .eq. 0) then 
                  call user_fn_initialize_model(enf_model_ids(i),
     &              user_sprngweld_init_model)
                elseif (j .eq. 1) then 
                  call user_fn_initialize_time_step
     &                                         (enf_model_ids(i),
     &                user_sprngweld_init_tstep)
                elseif (j .eq. 2) then 
                  call user_fn_init_node_state_data
     &                                         (enf_model_ids(i),
     &                user_sprngweld_init_sdata)
                elseif (j .eq. 3) then 
                  call user_fn_interaction_type(enf_model_ids(i),
     &                user_sprngweld_itype)
                elseif (j .eq. 4) then 
                  call user_fn_active          (enf_model_ids(i),
     &                user_sprngweld_is_active)
                elseif (j .eq. 5) then 
                  call user_fn_limit_force     (enf_model_ids(i),
     &                user_sprngweld_limit_force)
#else
                if (j .eq. 0) then 
                  call get_usersub( user_initialize_model,
     &              enf_model_sub_names(jj),isize)
                  call user_fn_initialize_model(enf_model_ids(i),
     &              user_initialize_model)
                elseif (j .eq. 1) then 
                  call get_usersub( user_initialize_time_step,
     &              enf_model_sub_names(jj),isize)
                  call user_fn_initialize_time_step
     &                                         (enf_model_ids(i),
     &                user_initialize_time_step)
                elseif (j .eq. 2) then 
                  call get_usersub( user_init_node_state_data,
     &              enf_model_sub_names(jj),isize)
                  call user_fn_init_node_state_data
     &                                         (enf_model_ids(i),
     &                user_init_node_state_data)
                elseif (j .eq. 3) then 
                  call get_usersub( user_interaction_type,
     &              enf_model_sub_names(jj),isize)
                  call user_fn_interaction_type(enf_model_ids(i),
     &                user_interaction_type)
                elseif (j .eq. 4) then 
                  call get_usersub( user_is_active,
     &              enf_model_sub_names(jj),isize)
                  call user_fn_active          (enf_model_ids(i),
     &                user_is_active)
                elseif (j .eq. 5) then 
                  call get_usersub( user_limit_force,
     &              enf_model_sub_names(jj),isize)
                  call user_fn_limit_force     (enf_model_ids(i),
     &                user_limit_force)
#endif
                else
                  write(*,*) "unrecognized enf model sub id : ", j
                endif
              enddo
            else
              if((enf_model_types(i) .eq. 4) 
     &       .or.(enf_model_types(i) .eq. 13))then
                r_emod_data(1) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+1)
                r_emod_data(2) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+2)
                i_emod_data(1) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+3)
                i_emod_data(2) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+4)
              else if( enf_model_types(i) .eq. 7 )then
                i_emod_data(1) = 
     &          enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+1)
                i_emod_data(2) = 
     &          enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+2)
                r_emod_data(1) =  
     &          enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+3)
                i_emod_data(3) = 
     &          enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+4)
                i_emod_data(4) = 
     &          enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+5)
              else if( enf_model_types(i) .eq. 8 )then
                i_emod_data(1) = 
     &          enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+1)
              else if( enf_model_types(i) .eq. 9 )then
                r_emod_data(1) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+1)
                r_emod_data(2) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+2)
                i_emod_data(1) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+3)              
              else if( enf_model_types(i) .eq. 10 )then
                i_emod_data(1) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+1)
                i_emod_data(2) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+2)
                r_emod_data(1) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+3)              
              else if( enf_model_types(i) .eq. 11 )then
                i_emod_data(1) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+1)
                i_emod_data(2) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+2)
                i_emod_data(3) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+3)              
                i_emod_data(4) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+4)
                i_emod_data(5) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+5)
                r_emod_data(1) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+6)              
                r_emod_data(2) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+7)              
                r_emod_data(3) = 
     &             enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+8)              
              else
                do iii=1,MAX_ENF_MODEL_DATA_VARS
                  r_emod_data(iii) = 
     &            enf_model_data(MAX_ENF_MODEL_DATA_VARS*(i-1)+iii)
                  i_emod_data(iii) = 0
                enddo
              endif
              call td_add_enf_model( enf_model_types(i), 
     &          enf_model_ids(i),i_emod_data(1), 
     &          r_emod_data(1), ierror )
            endif
          endif
        end do  
      else
        call read_restart_buffer_size( ibuf_size,icomm,
     *    my_proc_id,num_procs )
        if( ibuf_size.gt.MAX_RESTART_SIZE )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
          write(*,*) 'need a larger restart buffer'
#endif
#ifndef CONTACT_NO_MPI
          call MPI_FINALIZE(ierr)
#endif
          stop
        endif
        call read_restart_data( restart_buffer,icomm,
     *                          my_proc_id,num_procs )
	call build_search_restart( restart_buffer, 
     *                             node_ids,
     *                             face_ids,
     *                             elem_ids,
     *                             icomm, ierror )
        if( enforcement_type .gt. 0 )then
          call read_enf_restart_buffer_size( ibuf_size,icomm,
     *                                       my_proc_id,num_procs )
          if( ibuf_size.gt.MAX_RESTART_SIZE )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
            write(*,*) 'need a larger restart buffer'
#endif
#ifndef CONTACT_NO_MPI
            call MPI_FINALIZE(ierr)
#endif
            stop
          endif
          call read_enf_restart_data( restart_buffer,icomm,
     *                                my_proc_id,num_procs )
          if( enforcement_type .eq. 1 )then
            call build_td_enf_restart( restart_buffer, ierror )
          else if( enforcement_type .eq. 2 )then
            call build_gap_removal_restart( restart_buffer, ierror )
          else if( enforcement_type .eq. 3 )then
            call build_tied_kin_restart( restart_buffer, ierror )
          else if( enforcement_type .eq. 4 )then
            call build_vol_tran_restart( restart_buffer, ierror )
          else if( enforcement_type .eq. 5)then
            call build_mpc_eqns_restart( restart_buffer, ierror )
          endif
        end if
      end if
      
C find number of load steps  We can determine
C the number of load steps from the number of nodal variables
      num_pos_vars = num_nodal_vars
      if ( enforcement_type .eq. 1 ) then
         num_pos_vars = num_pos_vars - num_enf_type_1_vars
      endif
      if ( mod(num_pos_vars,idim) .ne. 0) then
         write (*,*) "FATAL ERROR: number of nodal vars in mesh "
         write (*,*) "  must be a multiple of the spatial dimension,"
         write (*,*) "  plus enforcement variables if enforcement"
         write (*,*) "  has been requested."
#ifndef CONTACT_NO_MPI
         call MPI_FINALIZE(ierr)
#endif
         stop
      endif
      if (isrch_type .eq. 0) then
         num_steps = num_pos_vars/idim+1
      else 
         num_steps = num_pos_vars/idim
      endif
      if ( num_steps .lt. num_load_steps) then
         write (*,*) "FATAL ERROR: number of load steps specified "
         write (*,*) "  in input file is greater than number of load"
         write (*,*) "  steps defined in mesh file."
#ifndef CONTACT_NO_MPI
         call MPI_FINALIZE(ierr)
#endif
         stop
      endif

      if ( istart_step .gt. num_load_steps ) then
         write (*,*) "FATAL ERROR: start step of ",istart_step," is"
	 write (*,*) "  larger than the number of steps specified, ", 
     &        num_steps 
#ifndef CONTACT_NO_MPI
         call MPI_FINALIZE(ierr)
#endif
         stop
      endif
      
C copy enforcement data and initial position
      if (enforcement_type .eq. 1) then
         ioffset = num_nodes_tot*num_pos_vars
         do i = 1, num_nodes_tot*num_enf_type_1_vars
            enforcement_vars(i) = vars_nodal(ioffset + i)
         enddo
      endif
      
C loop over the number of load steps we plan to run
      do icur_step = istart_step, num_load_steps
        write (*,*) ">> Executing step ", icur_step
        
C     if a new decomposition has been computed, then read the
C     new mesh data
        if (icur_step .gt. 1) then
          if (itop_update .eq. 2 .or. itop_update .eq. 3 ) then
            ilast_step = icur_step - 1
#ifdef DOTHIS
            call new_dlb (icomm, 
     &        idim, 
     &        num_node_blocks,
     &        node_block_types,
     &        nodes_per_block, 
     &        node_eids,
     &        node_ids, 
     &        vars_nodal, 
     &        num_nodal_vars,
     &        num_face_blocks, 
     &        iface_block_types,
     &        shell_block_offsets, 
     &        ifaces_per_block, 
     &        face_ids, 
     &        ifconn,
     &        pos_init, 
     &        attributes,
     &        num_attributes,
     &        shell_offset,
     &        num_comm_partners, 
     &        icomm_proc_id, 
     &        number_nodes_to_partner, 
     &        icomm_nodes, 
     &        isrch_type, 
     &        enforcement_type, 
     &        irestart_flag, 
     &        ilast_step,
     &        my_proc_id, 
     &        num_procs)
#endif
          endif
    
C modify the topology for any ignored faces
          inode_ignore_off = 0
          iface_ignore_off = 0
          ielem_ignore_off = 0
          do i = 1, icur_step - 1
            inode_ignore_off = inode_ignore_off + 
     &        num_node_ignore(i)
            iface_ignore_off = iface_ignore_off + 
     &        num_face_ignore(i)
            ielem_ignore_off = ielem_ignore_off + 
     &        num_elem_ignore(i)
          enddo
          if (itop_update .ne. 0) then
            call get_topology_birthdeath
     &           ( num_node_ignore(icur_step),
     &             node_ignore_list(inode_ignore_off+1),
     &             num_face_ignore(icur_step),
     &             face_ignore_list(iface_ignore_off+1),
     &             num_elem_ignore(icur_step),
     &             elem_ignore_list(ielem_ignore_off+1),
     &             num_node_deaths_per_block, 
     &	           node_deaths_global_ids,
     &             num_face_deaths_per_block, 
     &	           face_deaths_global_ids,
     &             num_element_deaths_per_block, 
     &	           element_deaths_global_ids,
     &             num_node_births_per_block, 
     &	           nbirths_exodus_ids,
     &	           nbirths_global_ids,
     &	           number_face_births_per_block, 
     &	           fbirths_global_ids,
     &	           ifconn,
     &	           number_element_births_per_block, 
     &	           ebirths_global_ids,
     &	           ieconn,
     &             my_proc_id,
     &             num_procs )
            call get_topology_comm_plan
     &           ( num_comm_partners,
     &             icomm_proc_id,
     &             number_nodes_to_partner,
     &             icomm_nodes,
     &             my_proc_id,
     &             num_procs )
            call get_topology_host_ids
     &           ( num_nodes_tot,
     &             node_eids,
     &             node_ids,
     &             num_faces_tot,
     &             face_ids,
     &             num_elems_tot,
     &             elem_ids )
            call get_topology_variables
     &           ( pos_init,
     &             vars_nodal,
     &             num_nodal_vars,
     &	           vars_element,
     &	           num_element_vars,
     &	           attributes,
     &	           num_attributes,
     &	           radius,
     &	           shell_thickness,
     &             shell_offset )
          endif
        endif

C   Set up the current configuration
        if( icur_step .gt. 1 .and. iuse_config .ge. 1 )then
          if( iuse_config .eq. 1 )then
            idx_var = (icur_step-2)*idim*num_nodes_tot
            do i = 1, num_nodes_tot
              idx_pos = idim*(i-1)
              pos1(idx_pos+1) = pos_init(idx_pos+1) + 
     &          vars_nodal(idx_var+i)
              pos1(idx_pos+2) = pos_init(idx_pos+2) + 
     &          vars_nodal(idx_var+num_nodes_tot+i)
              if (idim .eq. 3) then
                pos1(idx_pos+3) = pos_init(idx_pos+3) + 
     &            vars_nodal(idx_var+2*num_nodes_tot+i)
              endif
            enddo
          else if( iuse_config .eq. 2 )then
            idx_var = (2*(icur_step-1)-1)*idim*num_nodes_tot
            do i = 1,num_nodes_tot
              idx_pos = idim*(i-1)
              pos1(idx_pos+1) = vars_nodal(idx_var+i)
              pos1(idx_pos+2) = vars_nodal(idx_var+num_nodes_tot+i)
              pos1(idx_pos+3) = vars_nodal(idx_var+2*num_nodes_tot+i)
            end do
          end if
        else
          do i = 1, num_nodes_tot*idim
            pos1(i) = pos_init(i)
          enddo
        endif
    
C     now call the routine to update the toplogy if the topology is different
        if( icur_step .gt. 1 )then
          if (itop_update .ne. 0) then
            num_node_exports=0
            num_face_exports=0
            num_element_exports=0
            call update_search
     &           ( num_node_deaths_per_block, 
     &		   node_deaths_global_ids,
     &             num_face_deaths_per_block, 
     &		   face_deaths_global_ids,
     &             num_element_deaths_per_block, 
     &		   element_deaths_global_ids,
                   
     &             num_node_births_per_block, 
     &		   nbirths_exodus_ids,
     &		   nbirths_global_ids,
     &		   number_face_births_per_block, 
     &		   fbirths_global_ids,
     &		   ifconn,
     &		   number_element_births_per_block, 
     &		   ebirths_global_ids,
     &		   ieconn,
                   
     &		   num_node_exports,
     &             node_export_id_list,
     &             node_export_pid,
     &		   num_face_exports,
     &             face_export_id_list,
     &             face_export_pid,
     &  	   num_element_exports,
     &             element_export_id_list,
     &             element_export_pid,
                   
     &             num_nodes_tot,
     &             node_ids,
     &             num_faces_tot,
     &             face_ids,
     &             num_elems_tot,
     &             elem_ids,
     &         
     &		   num_comm_partners,
     &		   icomm_proc_id,
     &		   number_nodes_to_partner,
     &		   icomm_nodes,
                   
     &		   ierr )
          endif
        endif

C   Supply the current configuration
         index = 1
         do i=1,num_node_blocks
            call set_node_block_configuration(CURRENT_CONFIGURATION,
     *           node_block_ids(i),
     *           pos1(index),ierror )
            if( ierror .ne. 0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
               write(*,*) 
     *              "Error in Set_NodeBlk_Configuration: Error Code = ",
     *              ierror
#endif
               call number_of_search_errors( num_errors )
               do ii=1,num_errors
                  call get_search_error_message( ii,message )
                  call print_message( message )
               end do
               call cleanup_search()
#ifndef CONTACT_NO_MPI
               call MPI_FINALIZE(ierr)
#endif
               stop
            endif
            index = index + idim*nodes_per_block(i)
         end do
         
         index_so = 1
         do i = 1,num_face_blocks
           if( iface_block_types(i) .eq. 5 .or. 
     &         iface_block_types(i) .eq. 6 ) then
             call set_face_block_attributes( 1, i, 
     &             shell_thickness(index_so), ierror )
             call set_face_block_attributes( 2, i,
     &             shell_offset(index_so), ierror )
             index_so = index_so + ifaces_per_block(i)
           endif
         end do
	
         index_so = 1
         do i = 1,num_node_blocks
           if( node_block_types(i) .eq. 2  ) then
             call set_node_block_attributes( 2, i, 
     &             radius(index_so), ierror )
             index_so = index_so + nodes_per_block(i)
           endif
         end do

C   Do a restart using the new (SIERRA-friendly) restart capability
         if (irestart_flag.eq.4 .or. irestart_flag.eq.6) then
            call read_restart_buffer_size( ibuf_size,icomm,
     *           my_proc_id,num_procs )
            call read_restart_data( restart_buffer,icomm,
     *           my_proc_id,num_procs )
            num_gen_rs_vars  = restart_buffer(1)
            num_node_rs_vars = restart_buffer(2)
            num_edge_rs_vars = restart_buffer(3)
            num_face_rs_vars = restart_buffer(4)
            num_elem_rs_vars = restart_buffer(5)
            num_gen_re_vars  = restart_buffer(6)
            num_node_re_vars = restart_buffer(7)
            num_edge_re_vars = restart_buffer(8)
            num_face_re_vars = restart_buffer(9)
            num_elem_re_vars = restart_buffer(10)
            ibase = 11
            call search_implant_general_rsvars(restart_buffer(ibase), 
     *           ierr)
            ibase = ibase + num_gen_rs_vars
            do i = 1, num_node_rs_vars
               call search_implant_node_rsvars(i,
     *              restart_buffer((i-1)*num_nodes_tot+ibase), ierr)
            enddo
            ibase = ibase + num_node_rs_vars * num_nodes_tot
            do i = 1, num_face_rs_vars
               call search_implant_face_rsvars(i,
     *              restart_buffer((i-1)*num_faces_tot+ibase), 
     *              ierr)
            enddo
            ibase = ibase + num_face_rs_vars * num_faces_tot
            do i = 1, num_elem_rs_vars
               call search_implant_elem_rsvars(i,
     *              restart_buffer((i-1)*num_elems_tot+ibase), 
     *              ierr)
            enddo
            ibase = ibase + num_elem_rs_vars * num_elems_tot
            call search_complete_restart(ierr)

            if( enforcement_type .eq. 1) then
               call td_implant_general_rsvars(restart_buffer(ibase),
     *              ierr)
               ibase = ibase + num_gen_re_vars
               do i = 1, num_node_re_vars
                  call td_implant_node_rsvars(i,
     *                 restart_buffer((i-1)*num_nodes_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_node_re_vars * num_nodes_tot
               do i = 1, num_face_re_vars
                  call td_implant_face_rsvars(i,
     *                 restart_buffer((i-1)*num_faces_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_face_re_vars * num_faces_tot
               do i = 1, num_elem_re_vars
                  call td_implant_element_rsvars(i,
     *                 restart_buffer((i-1)*num_elems_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_elem_re_vars * num_elems_tot
            else if( enforcement_type .eq. 2) then
               call gap_implant_general_rsvars(restart_buffer(ibase), 
     *              ierr)
               ibase = ibase + num_gen_re_vars
               do i = 1, num_node_re_vars
                  call gap_implant_node_rsvars(i,
     *                 restart_buffer((i-1)*num_nodes_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_node_re_vars * num_nodes_tot
               do i = 1, num_face_re_vars
                  call gap_implant_face_rsvars(i,
     *                 restart_buffer((i-1)*num_faces_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_face_re_vars * num_faces_tot
               do i = 1, num_elem_re_vars
                  call gap_implant_element_rsvars(i,
     *                 restart_buffer((i-1)*num_elems_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_elem_re_vars * num_elems_tot
            else if( enforcement_type .eq. 3) then
               call tied_kin_implant_general_rsvars(
     *              restart_buffer(ibase), ierr)
               ibase = ibase + num_gen_re_vars
               do i = 1, num_node_re_vars
                  call tied_kin_implant_node_rsvars(i,
     *                 restart_buffer((i-1)*num_nodes_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_node_re_vars * num_nodes_tot
               do i = 1, num_face_re_vars
                  call tied_kin_implant_face_rsvars(i,
     *                 restart_buffer((i-1)*num_faces_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_face_re_vars * num_faces_tot
               do i = 1, num_elem_re_vars
                  call tied_kin_implant_element_rsvars(i,
     *                 restart_buffer((i-1)*num_elems_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_elem_re_vars * num_elems_tot
            else if( enforcement_type .eq. 4) then
               call vol_tran_implant_general_rsvars(
     *              restart_buffer(ibase), ierr)
               ibase = ibase + num_gen_re_vars
               do i = 1, num_node_re_vars
                  call vol_tran_implant_node_rsvars(i,
     *                 restart_buffer((i-1)*num_nodes_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_node_re_vars * num_nodes_tot
               do i = 1, num_face_re_vars
                  call vol_tran_implant_face_rsvars(i,
     *                 restart_buffer((i-1)*num_faces_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_face_re_vars * num_faces_tot
               do i = 1, num_elem_re_vars
                  call vol_tran_implant_element_rsvars(i,
     *                 restart_buffer((i-1)*num_elems_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_elem_re_vars * num_elems_tot
            else if( enforcement_type .eq. 5) then
               call mpc_eqns_implant_general_rsvars(
     *              restart_buffer(ibase), ierr)
               ibase = ibase + num_gen_re_vars
               do i = 1, num_node_re_vars
                  call mpc_eqns_implant_node_rsvars(i,
     *                 restart_buffer((i-1)*num_nodes_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_node_re_vars * num_nodes_tot
               do i = 1, num_face_re_vars
                  call mpc_eqns_implant_face_rsvars(i,
     *                 restart_buffer((i-1)*num_faces_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_face_re_vars * num_faces_tot
               do i = 1, num_elem_re_vars
                  call mpc_eqns_implant_element_rsvars(i,
     *                 restart_buffer((i-1)*num_elems_tot+ibase),
     *                 ierr)
               enddo
               ibase = ibase + num_elem_re_vars * num_elems_tot
            endif
         endif
C    Perform the Search
         ierror = 0
         if( isrch_type .eq. 0 )then
            call static_search_1_configuration( ierror )
            if( ierror.ne.0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
               write(*,*) "Error in Static_Search: Error Code = ",
     *              ierror
#endif
               call number_of_search_errors( num_errors )
               do ii=1,num_errors
                  call get_search_error_message( ii,message )
                  call print_message( message )
               end do
               call cleanup_search()  
#ifndef CONTACT_NO_MPI
               call MPI_FINALIZE(ierr)
#endif
               stop
            endif
         else

C   Set up the predicted coordinates
            if( iuse_config .ne. 2 )then
              idx_var = (icur_step-1)*idim*num_nodes_tot
              do i = 1, num_nodes_tot
                idx_pos = idim*(i-1)
                pos2(idx_pos+1) = pos_init(idx_pos+1) + 
     &              vars_nodal(idx_var+i)
                pos2(idx_pos+2) = pos_init(idx_pos+2) + 
     &              vars_nodal(idx_var+num_nodes_tot+i)
                if (idim .eq. 3) then
                  pos2(idx_pos+3) = pos_init(idx_pos+3) + 
     &                 vars_nodal(idx_var+2*num_nodes_tot+i)
                endif
              enddo
            else
              idx_var = (icur_step-1)*2*idim*num_nodes_tot
              do i = 1, num_nodes_tot
                pos2(idim*(i-1)+1) = pos1(idim*(i-1)+1) +
     &                 vars_nodal(idx_var+i)
                pos2(idim*(i-1)+2) = pos1(idim*(i-1)+2) +
     &                 vars_nodal(idx_var+num_nodes_tot+i)
                pos2(idim*(i-1)+3) = pos1(idim*(i-1)+3) +
     &                 vars_nodal(idx_var+2*num_nodes_tot+i)
              enddo
            endif

C   Supply the predicted configuration
            index = 1
            do i=1,num_node_blocks
               call set_node_block_configuration(
     *              PREDICTED_CONFIGURATION,
     *              node_block_ids(i),
     *              pos2(index),ierror)
               if( ierror .ne. 0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
                  write(*,*) 
     *                 "Error in Set_NodeBlk_Configuration: Error ",
     *                 "Code = ", ierror
#endif
                  call number_of_search_errors( num_errors )
                  do ii=1,num_errors
                     call get_search_error_message( ii,message )
                     call print_message( message )
                  end do
                  call cleanup_search()
#ifndef CONTACT_NO_MPI
                  call MPI_FINALIZE(ierr)
#endif
                  stop
               endif
               index = index + idim*nodes_per_block(i)
            end do
            if( isrch_type .eq. 1 )then
               call static_search_2_configuration( ierror )
               if( ierror.ne.0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
                  write(*,*) "Error in Static_Search_2: Error Code = ",
     &                 ierror
#endif
                  call number_of_search_errors( num_errors )
                  do ii=1,num_errors
                     call get_search_error_message( ii,message )
                     call print_message( message )
                  end do
                  call cleanup_search()
#ifndef CONTACT_NO_MPI
                  call MPI_FINALIZE(ierr)
#endif
                  stop
               endif
            else if( isrch_type .eq. 2 )then
               call dynamic_search_2_configuration( ierror )
               if( ierror.ne.0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
                  write(*,*) "Error in Dynamic_Search: Error Code = ",
     &                 ierror
#endif
                  call number_of_search_errors( num_errors )
                  do ii=1,num_errors
                     call get_search_error_message( ii,message )
                     call print_message( message )
                  end do
                  call cleanup_search()
#ifndef CONTACT_NO_MPI
                  call MPI_FINALIZE(ierr)
#endif
                  stop
               endif
            else if (isrch_type .eq. 3 )then
               do ii = 1, num_nodes_tot
                 enf_array_1(ii) = enforcement_vars(ii)
               enddo
               call dynamic_search_aug_2_config( 
     *              enf_array_1, enf_data_vars(1),
     *              enf_data_vars(2), ierror )
               if( ierror.ne.0 )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
                  write(*,*) "Error in Dynamic_Search: Error Code = ",
     &                 ierror
#endif
                  call number_of_search_errors( num_errors )
                  do ii=1,num_errors
                     call get_search_error_message( ii,message )
                     call print_message( message )
                  end do
                  call cleanup_search()
#ifndef CONTACT_NO_MPI
                  call MPI_FINALIZE(ierr)
#endif
                  stop
              endif
           endif
         endif
         
         call display()
         
C    Get the interactions
         call size_nodeface_interactions( num_interactions,idata_size )
         call get_nodeface_interactions( node_blk_ids,node_indx,
     *        node_entity_keys,
     *        iface_blk_ids,iface_indx,
     *        iface_proc,data)
         call output_nodeface_interactions( num_interactions, 
     *        node_blk_ids,node_indx,iface_blk_ids,iface_indx,
     *        iface_proc,data,idata_size)
         
         call size_nodesurface_interactions( num_interactions,
     *        idata_size )
         call get_nodesurface_interactions( node_blk_ids,node_indx,
     *        iface_blk_ids, data )
         call output_nodesurface_interactions( num_interactions, 
     *        node_blk_ids,node_indx,iface_blk_ids, data )
         
C    Perform the enforcement if required
         if( enforcement_type .eq. 1 )then
            do i=1,num_nodes_tot
              enf_array_1(i) = enforcement_vars(i)
              enf_array_2(i) = enforcement_vars(num_nodes_tot + i)
              enf_array_3(i) = enforcement_vars(num_nodes_tot*2 + i)
            enddo
            iterations = enf_data_vars(3)
            call set_td_iterations( iterations, ierror )
            if( enf_data_vars(4) .eq. 1.0 )then
              n1 = enf_data_vars(5)
              n2 = enf_data_vars(6)
              call set_td_enf_symm_nodes( n1, n2, ierror )
            endif
            call compute_td_contact_force( enf_data_vars(1), 
     *           enf_data_vars(2),
     *           enf_array_1, enf_array_2, enf_array_3,
     *           results, ierror )
            if( ierror .gt. 0 )then
               call number_of_td_errors(ne)
               do ii=1,ne
                  call get_td_error_message( ii, message )
                  call print_message( message )
               end do
            endif
         else if( enforcement_type .eq. 2 )then
            max_iter = enf_data_vars(1)
            triv_gap = enf_data_vars(2)
            call compute_gap_removal( max_iter, triv_gap, results, 
     *                                ierror )
            if( ierror .gt. 0 )then
              call number_of_gap_errors(ne)
              do ii=1,ne
                call get_gap_error_message( ii, message )
                call print_message( message )
              end do
            endif
         else if( enforcement_type .eq. 3 )then
           call compute_tied_kin_position( results, ierror )
           if( ierror .gt. 0 )then
             call number_of_tied_errors(ne)
             do ii=1,ne
               call get_tied_kin_error_message( ii, message )
               call print_message( message )
             end do
           endif
         else if( enforcement_type .eq. 4 )then
c          set number of element and node vars to transfer to all in genesis file
           num_node_vars = num_nodal_vars
           num_elem_vars = num_element_vars
c          pack up nodal data
           loop_size = 1
           do ii=1, num_nodes_tot
             do jj=1, num_node_vars
               enf_array_1(loop_size) = 
     &           vars_nodal((jj-1)*num_nodes_tot+ii)
               loop_size = loop_size + 1
             enddo
           enddo
c          pack up element data
           loop_size = 1
           do ii=1, num_elems_tot
             do jj=1, num_elem_vars 
               enf_array_4(loop_size) = 
     &           vars_element((jj-1)*num_elems_tot+ii)
               loop_size = loop_size + 1
             enddo
           enddo
           call compute_vol_tran( num_node_vars,
     &       num_elem_vars,
     &       enf_array_1,
     &       enf_array_4,
     &       enf_array_2,
     &       enf_array_5,
     &       enf_array_6,
     &       error )
           if( ierror .gt. 0 )then
             call number_of_voltrans_errors(ne)
             do ii=1,ne
               call get_vol_tran_error_message( ii, message )
               call print_message( message )
             end do
           endif
         else if( enforcement_type .eq. 5 )then
c ContactMPCs enforcement
           impc_num_s = 1
           call compute_mpc_eqns( impc_num_s, ierr )
           call get_num_mpceqns( nmpceq, ierr )
           call get_mpc_eqns(
     &       nmpceq,
     &       impc_array_1,
     &       impc_array_2,
     &       impc_array_3,
     &       impc_array_4,
     &       impc_array_5,
     &       impc_array_6,
     &       impc_array_7,
     &       rmpc_array_1,
     &       ierr )

         end if
#ifndef CONTACT_NO_EXODUS_OUTPUT
         time = icur_step
         call create_plot_file( icomm, my_proc_id, 
     *	                        icur_step, iout_exodus )
         call contact_exodus_output( iout_exodus, time, ierr )
	 call test_close( iout_exodus )
#endif

C end of load step loop
      enddo    
               
C do restart file if requested      
      if( irestart_flag.eq.2 .or. irestart_flag.eq.3 )then
         call search_restart_size( irestart_size )
         if( irestart_size.gt.MAX_RESTART_SIZE )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
            write(*,*) "Need a larger buffer for restart"
#endif
#ifndef CONTACT_NO_MPI
            call MPI_FINALIZE(ierr)
#endif
            stop
         end if
         call search_extract_restart_data( restart_buffer, ierror )
         call write_restart(irestart_size,restart_buffer,
     *        icomm,my_proc_id,num_procs)
         if( enforcement_type.eq.1 )then
           call td_enf_restart_size( irestart_size )
           if( irestart_size.gt.MAX_RESTART_SIZE )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
              write(*,*) "Need a larger buffer for restart"
#endif
#ifndef CONTACT_NO_MPI
              call MPI_FINALIZE(ierr)
#endif
              stop
           end if
           call td_enf_extract_restart( restart_buffer, ierror )
           call write_enf_restart( irestart_size,restart_buffer,
     *                             icomm,my_proc_id,num_procs )
         else if ( enforcement_type.eq.2 ) then
           call gap_removal_restart_size( irestart_size )
           if( irestart_size.gt.MAX_RESTART_SIZE )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
              write(*,*) "Need a larger buffer for restart"
#endif
#ifndef CONTACT_NO_MPI
              call MPI_FINALIZE(ierr)
#endif
              stop
           end if
           call gap_removal_extract_restart( restart_buffer, ierror )
           call write_enf_restart( irestart_size,restart_buffer,
     *                             icomm,my_proc_id,num_procs )
         else if ( enforcement_type.eq.3 ) then
           call tied_kin_restart_size( irestart_size )
           if( irestart_size.gt.MAX_RESTART_SIZE )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
              write(*,*) "Need a larger buffer for restart"
#endif
#ifndef CONTACT_NO_MPI
              call MPI_FINALIZE(ierr)
#endif
              stop
           end if
           call tied_kin_extract_restart( restart_buffer, ierror )
           call write_enf_restart( irestart_size,restart_buffer,
     *                             icomm,my_proc_id,num_procs )
         else if ( enforcement_type.eq.4 ) then
           call vol_tran_restart_size( irestart_size )
           if( irestart_size.gt.MAX_RESTART_SIZE )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
              write(*,*) "Need a larger buffer for restart"
#endif
#ifndef CONTACT_NO_MPI
              call MPI_FINALIZE(ierr)
#endif
              stop
           end if
           call vol_tran_extract_restart( restart_buffer, ierror )
           call write_enf_restart( irestart_size,restart_buffer,
     *                             icomm,my_proc_id,num_procs )
         else if ( enforcement_type.eq.5 ) then
           call mpc_eqns_restart_size( irestart_size )
           if( irestart_size.gt.MAX_RESTART_SIZE )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
             write(*,*) "Need a larger buffer for restart"
#endif
#ifndef CONTACT_NO_MPI
             call MPI_FINALIZE(ierr)
#endif
             stop
           end if
           call mpc_eqns_extract_restart( restart_buffer, ierror )
           call write_enf_restart( irestart_size,restart_buffer,
     *                             icomm,my_proc_id,num_procs )
         end if
      else if (irestart_flag.eq.5 .or. irestart_flag.eq.6) then
         call search_num_general_rsvars (num_gen_rs_vars)
         call search_num_node_rsvars (num_node_rs_vars)
         call search_num_edge_rsvars (num_edge_rs_vars)
         call search_num_face_rsvars (num_face_rs_vars)
         call search_num_elem_rsvars (num_elem_rs_vars)
         num_gen_re_vars  = 0
         num_node_re_vars = 0
         num_edge_re_vars = 0
         num_face_re_vars = 0
         num_elem_re_vars = 0
         if ( enforcement_type.eq.1) then
            call td_num_general_rsvars(num_gen_re_vars)
            call td_num_node_rsvars(num_node_re_vars)
            call td_num_edge_rsvars(num_edge_re_vars)
            call td_num_face_rsvars(num_face_re_vars)
            call td_num_element_rsvars(num_elem_re_vars)
         else if ( enforcement_type.eq.2 ) then
            call gap_num_general_rsvars(num_gen_re_vars)
            call gap_num_node_rsvars(num_node_re_vars)
            call gap_num_edge_rsvars(num_edge_re_vars)
            call gap_num_face_rsvars(num_face_re_vars)
            call gap_num_element_rsvars(num_elem_re_vars)
         else if ( enforcement_type.eq.3 ) then
            call tied_kin_num_general_rsvars(num_gen_re_vars)
            call tied_kin_num_node_rsvars(num_node_re_vars)
            call tied_kin_num_edge_rsvars(num_edge_re_vars)
            call tied_kin_num_face_rsvars(num_face_re_vars)
            call tied_kin_num_element_rsvars(num_elem_re_vars)
         else if ( enforcement_type.eq.4 ) then
            call vol_tran_num_general_rsvars(num_gen_re_vars)
            call vol_tran_num_node_rsvars(num_node_re_vars)
            call vol_tran_num_edge_rsvars(num_edge_re_vars)
            call vol_tran_num_face_rsvars(num_face_re_vars)
            call vol_tran_num_element_rsvars(num_elem_re_vars)
         else if ( enforcement_type.eq.5 ) then
            call mpc_eqns_num_general_rsvars(num_gen_re_vars)
            call mpc_eqns_num_node_rsvars(num_node_re_vars)
            call mpc_eqns_num_edge_rsvars(num_edge_re_vars)
            call mpc_eqns_num_face_rsvars(num_face_re_vars)
            call mpc_eqns_num_element_rsvars(num_elem_re_vars)
         endif
         irestart_size = 10
     &        + num_gen_rs_vars
     &        + num_node_rs_vars * num_nodes_tot
     &        + num_face_rs_vars * num_faces_tot
     &        + num_elem_rs_vars * num_elems_tot
     &        + num_gen_re_vars
     &        + num_node_re_vars * num_nodes_tot
     &        + num_face_re_vars * num_faces_tot
     &        + num_elem_re_vars * num_elems_tot
         if( irestart_size.gt.MAX_RESTART_SIZE )then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
            write(*,*) 'need a larger restart buffer'
#endif
#ifndef CONTACT_NO_MPI
            call MPI_FINALIZE(ierr)
#endif
            stop
         endif
         
         restart_buffer(1) = num_gen_rs_vars
         restart_buffer(2) = num_node_rs_vars
         restart_buffer(3) = num_edge_rs_vars
         restart_buffer(4) = num_face_rs_vars
         restart_buffer(5) = num_elem_rs_vars
         restart_buffer(6) = num_gen_re_vars
         restart_buffer(7) = num_node_re_vars
         restart_buffer(8) = num_edge_re_vars
         restart_buffer(9) = num_face_re_vars
         restart_buffer(10) = num_elem_re_vars

         ibase = 10+1
         call search_extract_general_rsvars(restart_buffer(ibase), 
     *        ierr)
         ibase = ibase + num_gen_rs_vars
         do i = 1, num_node_rs_vars
            call search_extract_node_rsvars(i,
     *           restart_buffer((i-1)*num_nodes_tot+ibase), ierr)
         enddo
         ibase = ibase + num_node_rs_vars * num_nodes_tot
         do i = 1, num_face_rs_vars
            call search_extract_face_rsvars(i,
     *           restart_buffer((i-1)*num_faces_tot+ibase), 
     *           ierr)
         enddo
         ibase = ibase + num_face_rs_vars * num_faces_tot
         do i = 1, num_elem_rs_vars
            call search_extract_elem_rsvars(i,
     *           restart_buffer((i-1)*num_elems_tot+ibase), 
     *           ierr)
         enddo
         ibase = ibase + num_elem_rs_vars * num_elems_tot
         
         if( enforcement_type .eq. 1) then
            call td_extract_general_rsvars(restart_buffer(ibase), 
     *           ierr)
            ibase = ibase + num_gen_re_vars
            do i = 1, num_node_re_vars
               call td_extract_node_rsvars(i,
     *              restart_buffer((i-1)*num_nodes_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_node_re_vars * num_nodes_tot
            do i = 1, num_face_re_vars
               call td_extract_face_rsvars(i,
     *              restart_buffer((i-1)*num_faces_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_face_re_vars * num_faces_tot
            do i = 1, num_elem_re_vars
               call td_extract_element_rsvars(i,
     *              restart_buffer((i-1)*num_elems_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_elem_re_vars * num_elems_tot
         else if( enforcement_type .eq. 2) then
            call gap_extract_general_rsvars(restart_buffer(ibase), 
     *           ierr)
            ibase = ibase + num_gen_re_vars
            do i = 1, num_node_re_vars
               call gap_extract_node_rsvars(i,
     *              restart_buffer((i-1)*num_nodes_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_node_re_vars * num_nodes_tot
            do i = 1, num_face_re_vars
               call gap_extract_face_rsvars(i,
     *              restart_buffer((i-1)*num_faces_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_face_re_vars * num_faces_tot
            do i = 1, num_elem_re_vars
               call gap_extract_element_rsvars(i,
     *              restart_buffer((i-1)*num_elems_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_elem_re_vars * num_elems_tot
         else if( enforcement_type .eq. 3) then
            call tied_kin_extract_general_rsvars(
     *           restart_buffer(ibase),ierr)
            ibase = ibase + num_gen_re_vars
            do i = 1, num_node_re_vars
               call tied_kin_extract_node_rsvars(i,
     *              restart_buffer((i-1)*num_nodes_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_node_re_vars * num_nodes_tot
            do i = 1, num_face_re_vars
               call tied_kin_extract_face_rsvars(i,
     *              restart_buffer((i-1)*num_faces_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_face_re_vars * num_faces_tot
            do i = 1, num_elem_re_vars
               call tied_kin_extract_element_rsvars(i,
     *              restart_buffer((i-1)*num_elems_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_elem_re_vars * num_elems_tot
         else if( enforcement_type .eq. 4) then
            call vol_tran_extract_general_rsvars(
     *           restart_buffer(ibase), ierr)
            ibase = ibase + num_gen_re_vars
            do i = 1, num_node_re_vars
               call vol_tran_extract_node_rsvars(i,
     *              restart_buffer((i-1)*num_nodes_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_node_re_vars * num_nodes_tot
            do i = 1, num_face_re_vars
               call vol_tran_extract_face_rsvars(i,
     *              restart_buffer((i-1)*num_faces_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_face_re_vars * num_faces_tot
            do i = 1, num_elem_re_vars
               call vol_tran_extract_element_rsvars(i,
     *              restart_buffer((i-1)*num_elems_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_elem_re_vars * num_elems_tot
         else if( enforcement_type .eq. 5) then
            call mpc_eqns_extract_general_rsvars(
     *           restart_buffer(ibase), ierr)
            ibase = ibase + num_gen_re_vars
            do i = 1, num_node_re_vars
               call mpc_eqns_extract_node_rsvars(i,
     *              restart_buffer((i-1)*num_nodes_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_node_re_vars * num_nodes_tot
            do i = 1, num_face_re_vars
               call mpc_eqns_extract_face_rsvars(i,
     *              restart_buffer((i-1)*num_faces_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_face_re_vars * num_faces_tot
            do i = 1, num_elem_re_vars
               call mpc_eqns_extract_element_rsvars(i,
     *              restart_buffer((i-1)*num_elems_tot+ibase),
     *              ierr)
            enddo
            ibase = ibase + num_elem_re_vars * num_elems_tot
         endif
         call write_restart(irestart_size,restart_buffer,
     *        icomm,my_proc_id,num_procs)
      endif
      
      if( enforcement_type .eq. 1 )then
         call cleanup_td_enforcement()
      else if( enforcement_type .eq. 2 )then
         call cleanup_gap_removal()
      else if( enforcement_type .eq. 3 )then
         call cleanup_tied_kinematics()
      else if( enforcement_type .eq. 4 )then
         call cleanup_vol_tran()
      else if( enforcement_type .eq. 5 )then
         call cleanup_mpc_eqns()
      end if
      call cleanup_search()
      
#ifndef CONTACT_NO_MPI
      call MPI_FINALIZE(ierr)
      if ( ierr .ne. 0) then
#if CONTACT_DEBUG_PRINT_LEVEL>=1
         write (*,*) "ERROR in finishing MPI."
#endif
         stop
      endif
#endif

      end

