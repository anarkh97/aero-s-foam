c-----------------------------------------------------------------------
      subroutine ehsfid(name1,name2)
      return
      end
c-----------------------------------------------------------------------
      subroutine ehufDO (i1,d1_w, d2_v, d1_p,g_ehfid,i2)
      write(*,*) 'Warning: ADIFOR Exception ehufDO called. Nothing done'
      return
      end
c-----------------------------------------------------------------------
      subroutine ehbfDO (i1,etop, ebot, d3_v, d1_p, d2_p,g_ehfid,name)
      character*(*)  name 
      write(*,*) 'Warning: ADIFOR Exception ehbfDO called.'
      return
      end
c-----------------------------------------------------------------------
      subroutine ehufdv(i,a1,a2,db,name)
      character*(*)  name 
      integer    i
      real*8     a1,a2,db
c      write(*,*) 
c      write(*,*) 'Warning: ADIFOR Exception ehufdv called.'
c      write(*,*) 'Derivative set to zero'
c      write(*,*) 
      db = 0.0d0
      return
      end
c-----------------------------------------------------------------------
      subroutine ehbfdv(i1,upvms, lwvms, d3_v, d1_p, d2_p,g_compvms,i2)   
      write(*,*) 'Warning: ADIFOR Exception ehbdv called. Nothing done'
      return
      end
