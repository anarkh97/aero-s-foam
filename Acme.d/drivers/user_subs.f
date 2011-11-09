c $Id: user_subs.f,v 2002.8 2004/06/18 17:47:04 mwglass Exp $

      subroutine register_user_subs()

      implicit none
      integer isize

      external user_adhesion_init_model
      external user_adhesion_itype
      external user_adhesion_is_active
      external user_adhesion_limit_force
      
      external user_sprngweld_init_model
      external user_sprngweld_init_tstep
      external user_sprngweld_init_sdata
      external user_sprngweld_itype
      external user_sprngweld_is_active
      external user_sprngweld_limit_force
      
       
c    SUN limit 31 char:123456789a123456789b123456789c
      isize =     len('user_adhesion_init_model')
      call reg_usersub(user_adhesion_init_model, 
     *                'user_adhesion_init_model',isize)
      isize =     len('user_adhesion_itype')
      call reg_usersub(user_adhesion_itype, 
     *                'user_adhesion_itype',isize)
      isize =     len('user_adhesion_is_active')
      call reg_usersub(user_adhesion_is_active, 
     *                'user_adhesion_is_active',isize)
      isize =     len('user_adhesion_limit_force')
      call reg_usersub(user_adhesion_limit_force, 
     *                'user_adhesion_limit_force',isize)
     
      isize =     len('user_sprngweld_init_model')
      call reg_usersub(user_sprngweld_init_model, 
     *                'user_sprngweld_init_model',isize)
      isize =     len('user_sprngweld_init_tstep')
      call reg_usersub(user_sprngweld_init_tstep, 
     *                'user_sprngweld_init_tstep',isize)
      isize =     len('user_sprngweld_init_sdata')
      call reg_usersub(user_sprngweld_init_sdata, 
     *                'user_sprngweld_init_sdata',isize)
      isize =     len('user_sprngweld_itype')
      call reg_usersub(user_sprngweld_itype, 
     *                'user_sprngweld_itype',isize)
      isize =     len('user_sprngweld_is_active')
      call reg_usersub(user_sprngweld_is_active, 
     *                'user_sprngweld_is_active',isize)
      isize =     len('user_sprngweld_limit_force')
      call reg_usersub(user_sprngweld_limit_force, 
     *                'user_sprngweld_limit_force',isize)
      
      return
      end
     
      subroutine user_adhesion_init_model(enf, id, rdata, 
     *                                          idata, istat)
      integer istat, id, idata(*)
      real*8  rdata(*), value
      istat = 1
      call userquery_table_last_abscissa(enf, idata(1), value)
      rdata(1) = value
      return
      end
      
      subroutine user_adhesion_itype(enf, zni, id, rdata, 
     *                              idata, itype, index, istat)
      integer istat, id, itype, index, idata(*)
      real*8  rdata(*)
      istat = 3
      return
      end
  
      subroutine user_adhesion_is_active(enf, zni, id, rdata, idata, 
     *                                   itype, index, gap, istat)
      integer istat, id, itype, index, idata(*)
      real*8  gap, rdata(*)
      if (gap.lt.rdata(1)) then
        istat = 1
      else
        istat = 0
      endif  
      return
      end

      subroutine user_adhesion_limit_force(enf, zni, id, rdata, idata,
     *                                     itype, index, gap, rel_disp, 
     *                                     slip, normal, dt, area, 
     *                                     force, istat)
      integer istat, id, itype, index, idata(*)
      real*8 gap, dt, area, g, f_n, f_np, value
      real*8 rel_disp(*), slip(*), normal(*), force(*), rdata(*)
      
      if (itype.eq.0) then
        g   = 0.0
        f_n = 0.0
        do i=1,3
          g   = g+rel_disp(i)*normal(i)
          f_n = f_n+force(i)*normal(i)
        enddo
        if (f_n.ge.0.0) then
c         compression
          istat = 0
        else 
c         tension
          if (g.lt.rdata(1)) then
            call userquery_table_interpolate(enf, idata(1), g, value)
            f_np = -value
            do i=1,3
              force(i) = normal(i)*f_np*area
            enddo
          else
            do i=1,3
              force(i) = 0.0
            enddo
          endif
          istat = 1
        endif
      else
        force(1) = 0.0
        force(2) = 0.0
        force(3) = 0.0
        istat = 1
      endif
      return
      end
      

     
      subroutine user_sprngweld_init_model(enf, id, rdata, 
     *                                     idata, istat)
      integer istat, id, idata(*)
      real*8  rdata(*), value
      call userquery_table_last_abscissa(enf, idata(2), value)
      rdata(2) = value
      call userquery_table_last_abscissa(enf, idata(3), value)
      rdata(3) = value
      istat = 1
      return
      end
      
      subroutine user_sprngweld_init_tstep(enf, id, rdata, 
     *                                         idata, istat)
      integer idata(*), id, istat, index
      real*8  rdata(*), value
      call userquery_number_of_nodes(enf, nnodes)
      value = 0.0
      index = 4
      do i=1,nnodes
        call userset_node_state_data(enf, id, i, index, value)
      enddo
      istat = 1
      return
      end
      
      subroutine user_sprngweld_init_sdata(enf, id, rdata, idata,
     *                                          sdata, istat)
      integer idata(*), id, istat
      real*8  rdata(*), sdata(*)
      sdata(1) = idata(1)
      sdata(2) = 0.0
      sdata(3) = 0.0
      sdata(4) = 0.0
      istat = 1
      return
      end
      
      subroutine user_sprngweld_itype(enf, zni, id, rdata, 
     *                           idata, itype, index, istat)
      integer istat, id, itype, index, idata(*), indx
      real*8  rdata(*), value
      indx = 1
      call userquery_node_state_data(enf, id, index, indx, value)
      if (value.gt.0.0) then
        istat = 4
      else 
        istat = -1
      endif
      return
      end
  
      subroutine user_sprngweld_is_active(enf, zni, id, rdata, idata, 
     *                                   itype, index, gap, istat)
      integer istat, id, itype, index, idata(*), indx
      real*8  gap, rdata(*), value
      indx = 1
      call userquery_node_state_data(enf, id, index, indx, value)
      if (value.gt.0.0) then
        istat = 1
      else 
        istat = -1
      endif
      return
      end

      subroutine user_sprngweld_limit_force(enf, zni, id, rdata, idata,
     *                                      itype, node_index, 
     *                                      gap, rel_disp, 
     *                                      slip, normal, dt, area, 
     *                                      force, istat)
      integer istat, id, itype, node_index, index, idata(*)
      real*8 gap, dt, area, f_n, value, s_i
      real*8 strength, already_decremented, s_n, s_t(3)
      real*8 mag_f_t, mag_s_t, norm_ratio, tang_ratio
      real*8 failure_crit, scale, f_n_fail, f_t_fail
      real*8 rel_disp(*), slip(*), normal(*), force(*), rdata(*)
      
      s_i = idata(1)
      if (itype.eq.0) then
        index = 1
        call userquery_node_state_data(enf, id, node_index, 
     *                                 index, strength)
        if (strength.gt.0.0) then
          s_n = 0.0
          do i=1,3
            s_n = s_n+rel_disp(i)*normal(i)
          enddo
          do i=1,3
            s_t(i) = rel_disp(i) - s_n*normal(i)
          enddo
          mag_s_t = 0.0
          do i=1,3
            mag_s_t = mag_s_t+s_t(i)*s_t(i)
          enddo
          mag_s_t = sqrt(mag_s_t)
          index = 4
          call userquery_node_state_data(enf, id, node_index, index, 
     *                                    already_decremented)
          if (strength.eq.s_i) then
            f_n = 0.0
            do i=1,3
              f_n = f_n+force(i)*normal(i)
            enddo
            if (s_n.gt.0.0) then
              call userquery_table_interpolate(enf, idata(2), 
     *                                         s_n, value)
              f_n = -value
            endif
            norm_ratio = abs(max(s_n,0.0d0)/rdata(2))
            call userquery_table_interpolate(enf, idata(3), 
     *                                       mag_s_t, value)
            mag_f_t      = -value
            tang_ratio   = abs(mag_s_t/rdata(3))
            failure_crit = norm_ratio**rdata(1) + 
     *                     tang_ratio**rdata(1)
            if (failure_crit.ge.1.0) then
              if( already_decremented.eq.0.0 ) then
                strength = strength - 1.0
                already_decremented = 1.0
                index = 1
                call userset_node_state_data(enf, id, node_index, 
     *                                       index, strength)
                index = 4
                call userset_node_state_data(enf, id, node_index, 
     *                                       index, already_decremented)
              endif
              index = 2
              call userset_node_state_data(enf, id, node_index, 
     *                                     index, f_n)
              index = 3
              call userset_node_state_data(enf, id, node_index, 
     *                                     index, mag_f_t)
            endif
            if (mag_s_t.gt.1.0e-10) then
              do i=1,3
                force(i) = f_n*normal(i) + (mag_f_t/mag_s_t)*s_t(i)
              enddo
            else
              do i=1,3
                force(i) = f_n*normal(i)
              enddo
            endif
          else
            if (already_decremented.eq.0.0) then
              strength = strength - 1.0
              already_decremented = 1.0
              index = 1
              call userset_node_state_data(enf, id, node_index, 
     *                                     index, strength)
              index = 4
              call userset_node_state_data(enf, id, node_index, 
     *                                     index, already_decremented)
            endif
            scale = (strength+1)/s_i
            index = 2
            call userquery_node_state_data(enf, id, node_index, 
     *                                      index, f_n_fail)
            if (s_n.gt.0.0) then
              f_n = scale*f_n_fail
            else
              f_n = 0.0
              do i=1,3
                f_n = f_n+force(i)*normal(i)
              enddo
            endif
            if (mag_s_t.gt.1.0e-10) then 
              index = 3
              call userquery_node_state_data(enf, id, node_index, 
     *                                        index, f_t_fail)
              do i=1,3
                force(i) = f_n*normal(i) + 
     *                     (scale*f_t_fail/mag_s_t)*s_t(i)
              enddo 
            else
              do i=1,3
                force(i) = f_n*normal(i)
              enddo
            endif
            if (strength.le.0.0) then
              call userset_nfi_failure_tied(enf, id, zni)
            endif
          endif
        endif
      else
        force(1) = 0.0
        force(2) = 0.0
        force(3) = 0.0
      endif
      istat = 1
      return
      end

      
