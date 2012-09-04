      subroutine nlpsal (ialg,
     &                   x,xl,xu,f,g,df,dg,
     &                   u,wa,dfloc,xupp,xlow,dsave,xsave,
     &                   nbd,iwa,iactive,isave,
     &                   n,m,me,
     &                   acc,rpen,penfac,uinit,accloc,bscal,relbo,
     &                   nbfgs,lenwa,itmax,maxloc,
     &                   ifail,iprint,fname,ifnsize)
c
      implicit real*8 (a-h,o-z)
      character fname*(*),cdum*1
c
c     open output file 
c
      IOUT= 45
      OPEN (IOUT,FILE=FNAME(1:IFNSIZE))
    1 READ (IOUT,1000,END=2) CDUM
      GOTO 1
    2 CONTINUE
c      
 1000 format(a1) 
c
c     call optimizer
c
      floc = 0.0d0
c
      call lbfgs(x,xl,xu,f,g,df,dg,
     &           floc,u,wa,dfloc,xupp,xlow,dsave,xsave,
     &           nbd,iwa,iactive,isave,
     &           n,m,me,nbfgs,lenwa,
     &           acc,rpen,penfac,uinit,accloc,bscal,relbo,
     &           itmax,maxloc,iout,ialg)
c
      return
      end

c=======================================================================

      subroutine lbfgs (x,xl,xu,f,g,df,dg,
     &                  floc,u,wa,dfloc,xupp,xlow,dsave,xsave,
     &                  nbd,iwa,iactive,isave,
     &                  n,m,me,nbfgs,lenwa,
     &                  acc,rpen,penfac,uinit,accloc,bscal,relbo,
     &                  itmax,maxloc,iout,ialg)
c 
c     Declare the variables needed by the code.
c     A description of all these variables is given at the end of 
c     the driver.

      implicit real*8 (a-h,o-z)

c     global arrays

      dimension nbd(n),u(m),iactive(m)

      dimension x(n),xl(n),xu(n),dfloc(n),xupp(n),xlow(n),xsave(n)

c     local variables and ararays

      character*60     task, csave
      logical          lsave(4)

c     nbd(i)=0 if x(i) is unbounded,
c            1 if x(i) has only a lower bound,
c            2 if x(i) has both lower and upper bounds, 
c            3 if x(i) has only an upper bound.
c
c     assume that variables has upper and lower bounds

      do i=1,n
        nbd(i) = 2
      enddo

c     initialize lagrange multipliers

      if (uinit.lt.0) then
        open(43,file='sal.data')
        do j=1,m
          read(43,*) u(j)
        enddo
        close(43)       
      else
        do j=1,m
          u(j)  = uinit
        enddo
      endif

c     initialize set of active constraints

      do j=1,m
        iactive(j)  = 1
      enddo

c     set convergence criteria of local minimizer extremly high

      factr  = 1.0d+1
      pgtol  = 1.0d-12

c     no output of local minimizer 

      iprint = -1

c     initialize lower and upper bounds

      do i=1,n
        xlow(i)  = x(i)
        xupp(i)  = x(i)
        xsave(i) = x(i)
      enddo

c     adjust lower and upper bounds

      do i=1,n
        delmax  = relbo*xu(i)-xl(i)
        xlow(i) = max(xl(i),x(i)-delmax)
        xupp(i) = min(xu(i),x(i)+delmax)
      enddo

c     start the iteration by initializing task.
 
      task = 'START'

      niter = 1
c
      nfunc = 0
      ngrad = 0
      inipr = 0

c     ------- Printing initial configuration -----

      write(iout,*)
      write(iout,*) 'NLPSAL - Algorithm' 
      write(iout,*) '========================' 
      write(iout,*)
      CALL FLUSH(IOUT)

c     ------- The beginning of the loop ----------

  111 continue
      
c     Printing initial configuration

      if (inipr.eq.1) then
        call lbfgsprint(n,m,x,u,f,floc,g,df,dg,0,iout)
      endif

      inipr=inipr+1

c     This is the call to the L-BFGS-B code.
 
      call setulb(n,nbfgs,x,xlow,xupp,nbd,floc,dfloc,factr,pgtol,wa,iwa,
     +            task,iprint,csave,lsave,isave,dsave)


      if (task(1:2) .eq. 'FG') then

c        we build a uncontraint minimization problem

         call lbfgsfunc(m,me,n,f,g,df,dg,x,xl,xu,xlow,xupp,xsave,
     &                  u,floc,dfloc,iactive,
     &                  acc,maxloc,rpen,penfac,accloc,bscal,
     &                  relbo,niter,task,nfunc,ngrad,iout,ialg)     

c        stop if constraint problem has converged

         if (task(1:4).eq.'STOP') goto 999
      
c        Go back to the minimization routine.

         goto 111

      elseif (task(1:5) .eq. 'NEW_X') then
 
c         print intermediate design

          call lbfgsprint(n,m,x,u,f,floc,g,df,dg,niter,iout)

c         check for maximum iterations and maximum total function evaluations

          if (niter.gt.itmax .or. nfunc.gt.itmax*5) goto 999

          niter = niter + 1

c        Go back to the minimization routine.

         goto 111

      else

c        check for violation of constraints otherwise stop

         task = 'CHECK'

         call lbfgsfunc(m,me,n,f,g,df,dg,x,xl,xu,xlow,xupp,xsave,
     &                  u,floc,dfloc,iactive,
     &                  acc,maxloc,rpen,penfac,accloc,bscal,
     &                  relbo,niter,task,nfunc,ngrad,iout,ialg)          

c        stop if constraint problem has converged

         if (task(1:4).eq.'STOP') goto 999
      
c        Go back to the minimization routine.

         goto 111

      endif

c     ---------- The end of the loop -------------

 999  continue

c     Final output

      write(iout,*)
      write(iout,*) 'FINAL CONVERGENCE OUTPUT' 
      write(iout,*) '========================' 
      write(iout,*)
     
      if (niter.eq.itmax)
     &  write(iout,*) 'WARNING: maximum number of iterations exceeded'
      if (nfunc.eq.itmax*5)
     &  write(iout,*) 
     &  'WARNING: maximum number of function evaluations exceeded'

      call lbfgsprint(n,m,x,u,f,floc,g,df,dg,niter,iout)

      write(iout,*)
      write(iout,*) 'NUMBER OF ITERATIONS    :',niter
      write(iout,*) 'NUMBER OF FUNCTION CALLS:',nfunc
      write(iout,*) 'NUMBER OF GRADIENT CALLS:',ngrad
      write(iout,*)
      call flush(iout)

      return
      end

c=======================================================================

      subroutine lbfgsfunc(m,me,n,f,g,df,dg,x,xl,xu,xlow,xupp,xsave,
     &                     u,floc,dfloc,iactive,
     &                     acc,maxloc,rpen,penfac,accloc,bscal,
     &                     relbo,niter,task,nfunc,ngrad,iout,ialg)

      implicit real*8 (a-h,o-z)
c
      parameter(zero=0.0d0)
      parameter(climit=1.0d-18)
c
      character*(*) task
c
      dimension x(n), xl(n), xu(n),  xlow(n), xupp(n), xsave(n)
      dimension g(n), df(n), dg(m,n)
      dimension u(m), dfloc(n) ,f(*)
c
      dimension iactive(m)
c
      save fold,flast,nold,nloc,iskip,dlcini
c
      if (nfunc.eq.0) then
        nold  = 0
        nloc  = 0
        iskip = 0
        flast = 0
      endif
c
      write(iout,*) 'entering lbfgsfunc with ',task(1:5)
      write(iout,*) 'fold ',fold
      write(iout,*) 'flast ',flast
c
c     initialize control parameters
c
      inew   = 0
      icheck = 0
c
      if (task(1:5).eq.'CHECK') icheck=1
c
c     local evaluation counter
c
      if (nloc.eq.maxloc) icheck=1
      nloc  = nloc + 1
c
c     check whether to evaluate functions & gradients
c     (when restart the first evaluation can be skipped)
c
      if (iskip.gt.0) goto 1

c       reset variables if last optimization procedure was not
c       successfull
c
   11 continue
c
      if (icheck.gt.0) then

        if (floc.gt.flast) then
          write(iout,*) 'maximum number of iteration step reached'
          write(iout,*) 'not progress made - reset variables'
          write(iout,*)
          write(iout,*) 'previous solution:',fold
          write(iout,*) 'new solution	  :',flast

          do i=1,n
            x(i) = xsave(i)
          enddo
        endif

      endif
c
c     evaluate functions
c
      call func(ialg,niter)
      nfunc = nfunc + 1
c
      write(iout,*) 'objective',f(1)
      do i=1,m
        write(iout,*) 'constraint ',i,' :',g(i)
      enddo
      call flush(iout)
c
c     determine active constraints
c
      do j=me+1,m
        iactive(j) = 1
        val = max(zero,u(j)-rpen*g(j))
        if (abs(val).lt.climit) iactive(j) = 0
      enddo
c
c     evaluate gradients
c
      call grad(ialg,iactive)
      ngrad = ngrad + 1

c      write(iout,*)
c      do i=1,n
c        write(iout,*) 'gradient of objective wrt. var',i,' :',df(i)
c      enddo
c      do j=1,m
c        if(iactive(j).gt.0) then 
c          do i=1,n
c            write(iout,*) 'gradient of constraint ',j,' wrt. var',
c     *                     i,' :',dg(j,i)
c          enddo
c        else
c          write(iout,*) 'no gradients, constraint ',j,' not active'
c        endif
c      enddo
c      write(iout,*)
c      call flush(iout)
c
    1 continue
c
c     build augmented Lagrangian
c
      floc = f(1) 
c
      do i=1,me
        floc  = floc + u(i)*g(i) + rpen*g(i)*g(i)/2.0d0
      enddo
      do i=me+1,m
        if(iactive(i).ne.0) then
          val  = max(0.0,u(i)-rpen*g(i))
          floc = floc + (val*val-u(i)*u(i))/rpen/2.0d0
        endif
      enddo
c
c     build gradients of augmented Lagrangian
c
      do i=1,n
        sum = df(i)
        do j=1,me
          sum = sum + u(j)*dg(j,i) + rpen*g(j)*dg(j,i)
        enddo
        do j=me+1,m
          if(iactive(j).ne.0) then
            val  = max(zero,u(j)-rpen*g(j))
            sum  = sum - dg(j,i)*val
          endif
        enddo
        dfloc(i) = sum
      enddo
c
c     scale augmented Lagrangian and its gradients
c
      dlloc = zero
c
      floc = bscal * floc
c
      do i=1,n
        dfloc(i) = bscal * dfloc(i)
        if (x(i).ne.xlow(i) .and. x(i).ne.xupp(i))
     &    dlloc = dlloc + dfloc(i)*dfloc(i)
      enddo
c
      dlloc=sqrt(dlloc)
c
      if (nloc.eq.1) dlcini=dlloc
c
      write(iout,'(a,2i5,2e20.10,i3)') 'IN LINE SEARCH - AUGMENTED:',
     &   nloc,nfunc,FLOC,f(1),iskip
      write(iout,'(a,2e20.10)') 'IN LINE SEARCH - CONVERGENCE:',
     &   abs((floc-fold)/floc),dlloc/dlcini

c
      call flush(iout)
c
c     initialize flast
c
      if (nfunc.eq.1) flast = floc
c
c     save last best solution
c
      if (floc.lt.flast.or.nloc.eq.0) then       
        write(iout,*) 'save solution for internal restart'
        flast=floc
        write(iout,*) 'new best solution :',flast

        do i=1,n
          xsave(i) = x(i)
        enddo
      endif
c
c     reset skip flag
c
      if (iskip.gt.0 .and. inew.eq.0) iskip = 0
c
c     check convergence of unconstraint problem 
c     and update lagrange multipliers
c
      if ( (niter.ne.nold .and. inew.eq.0) .or. icheck.eq.1) then

        delf=abs((floc-fold)/floc)
        delf=dlloc/dlcini

        if( delf.lt.accloc .or. icheck.eq.1) then
c
          if (icheck.eq.0 .and. floc.gt.flast) then
            icheck=1
            goto 11
          endif
c
          write(iout,*) 'Update of Lagrange multipliers'
          write(iout,*) '=============================='
          write(iout,*)
c
          cv=zero
          dl=zero
          du=zero
          db=zero

c         evaluate constraint violation

          do i=1,me
            cv   =  cv + abs(g(i))
          enddo

          do i=me+1,m
            cv   = cv + max(zero,-g(i))
          enddo

c         update lagrange multiplier

c          if (cv.eq.0.0) then
            do i=1,me
              val  =  rpen*g(i)
              du   =  du + (u(i)-val)*(u(i)-val)
              u(i) =  u(i) + val
            enddo
c
            do i=me+1,m
              if(iactive(i).ne.0) then
                val  = max(zero,u(i)-rpen*g(i))
                du   = du + (u(i)-val)*(u(i)-val)
                u(i) = val
              endif
            enddo

            du = sqrt(du)

c          endif

c         increase penalty factor

c          if (cv.gt.0.0) 
          rpen = rpen * penfac

c         adjust lower and upper bounds

          do i=1,n
            delmax  = relbo*xu(i)-xl(i)
            valo    = max(xl(i),x(i)-delmax)
            valu    = min(xu(i),x(i)+delmax)
            db      = db + (xlow(i) - valo)**2 + (valu - xupp(i))**2
            xlow(i) = valo
            xupp(i) = valu
          enddo

c         determine norm of gradient of Lagrangian wrt x

          do i=1,n
            sum = df(i)
            do j=1,me
              sum = sum + dg(j,i)*u(j)
            enddo
            do j=me+1,m
              if(iactive(j).ne.0) then
                sum = sum - dg(j,i)*u(j)
              endif
            enddo
            dl = dl + sum*sum
          enddo

c         print statistics of inner minimization

          write(iout,*)
          write(iout,*) 'ITERATION  IN LOCAL OPTIMIZER  :',nloc
          write(iout,*) 'CONVERGENCE IN LOCAL OPTIMIZER :',delf
          write(iout,*) 'SUM OF CONSTRAIN VIOLATION     :',cv
          write(iout,*) 'NORM OF GRADIENT OF LAGRANGIAN :',sqrt(dl)
          write(iout,*) 'NORM OF DLAMBDA                :',du
          write(iout,*) 'NORM OF DBOUNDS                :',sqrt(db)
          write(iout,*) 'PENALTY FACTOR                 :',rpen
          write(iout,*)
          call flush(iout)

          istop = 0

          if (du  .lt.acc)  istop=istop+1
          if (delf.lt.acc)  istop=istop+1
          if (db  .lt.acc)  istop=istop+1

          if (istop.eq.3) then
            task = 'STOP'
            return
          endif

          if (cv.gt.acc .or. dl.gt.acc .or. db.gt.acc) then
            task   = 'START'
            inew   = 1
            icheck = 0
            nloc   = 0
            iskip  = 1
            goto  1
          else
            task = 'STOP'
            return
          endif

        endif

      endif
c
      if (niter.ne.nold .or. inew.eq.1) then
        fold=floc
        nold=niter
      endif

      if (inew.eq.1) flast=floc
c
      return
      end
 
c=======================================================================

      subroutine lbfgsprint(n,m,x,u,f,floc,g,df,dg,niter,iout)

      implicit real*8 (a-h,o-z)

      dimension x(n)
      dimension u(m),g(m),df(n),dg(m,n)

      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,'(A,I5)') 'ITERATION ',NITER
      WRITE(IOUT,*)
      WRITE(IOUT,'(A)') 'VARIABLES:'
      DO IX=1,N
        WRITE (IOUT,'(A,I5,E20.10)') 'VARIABLE ',IX,X(IX)
      ENDDO
      WRITE(IOUT,*)
      WRITE(IOUT,'(A,E20.10)') 'OBJECTIVE          :',F
      WRITE(IOUT,'(A,E20.10)') 'AUGMENTED OBJECTIVE:',FLOC
c      DO IX=1,N
c        WRITE (IOUT,'(A,I5,E20.10)') 'GRADIENT VARIABLE ',IX,DF(IX)
c      ENDDO
      WRITE(IOUT,*)
      DO JX=1,M
        WRITE(IOUT,'(A,I5,2E20.10)') 'CONSTRAINT:',JX,G(JX),U(JX)
c        DO IX=1,N
c          WRITE (IOUT,'(A,I5,E20.10)') 'GRADIENT VARIABLE ',IX,DG(JX,IX)
c        ENDDO
      ENDDO
      WRITE(IOUT,*)
      CALL FLUSH(IOUT)

      return
      end
