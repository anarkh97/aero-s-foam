! ============================================
! eigen values and vectors of symmetric matrix
! ============================================
!
! ref:
! ---
! Harwell Subroutine Library: ea06 / ea07
! http://www.cse.clrc.ac.uk/nag/hsl/contents.shtml
!
! note:
! ----
! call ea06cd(a,value,vector,m,ia,iv,w) : compute eigen values and vectors
! call ea07cd(a,value,m,ia,w) : compute eigen values only
!
! =========================================================================================================



* *******************************************************************
* copyright (c) 1970 hyprotech uk
* all rights reserved.
*
* none of the comments in this copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* this package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* all use is subject to licence. for full details of an hsl archive
* licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* please note that for an hsl archive licence:
*
* 1. the package must not be copied for use by any other person.
*    supply of any part of the library by the licensee to a third party
*    shall be subject to prior written agreement between aea
*    hyprotech uk limited and the licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. all information on the package is provided to the licensee on the
*    understanding that the details thereof are confidential.
* 3. all publications issued by the licensee that include results obtained
*    with the help of one or more of the packages shall acknowledge the
*    use of the packages. the licensee will notify the numerical analysis
*    group at rutherford appleton laboratory of any such publication.
* 4. the packages may be modified by or on behalf of the licensee
*    for such use in research applications but at no time shall such
*    packages or modifications thereof become the property of the
*    licensee. the licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. neither cclrc nor hyprotech uk limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of packages by the licensee.
*
* *******************************************************************
*
*######date 26 oct 1992
c       toolpack tool decs employed.
c       variables declared but not used have been removed.
c       parameter statement added.
c       m22 removed.
c       modified for obsolescent features (feb 1997)
c
      subroutine ea06cd(a,value,vector,m,ia,iv,w)
      double precision one,zero
      parameter (one=1.0d0,zero=0.0d0)
      integer ia,iv,m
      double precision a(ia,m),value(m),vector(iv,m),w(5*m)
      double precision aiw,pp
      integer i,ida,j,l,m1,m2
      external ea06dd,ea08cd,mc04bd
      do 20 i = 1,m
        do 10 j = 1,m
          vector(j,i) = zero
   10   continue
        vector(i,i) = one
   20 continue
      m1 = m + 1
      m2 = m1 + m
      w(1) = a(1,1)
      if(m-2.lt.0) go to 50
      if(m-2.gt.0) go to 40
      w(2) = a(2,2)
      w(4) = a(2,1)
      go to 50
   40 continue
      call mc04bd(a,w,w(m1),m,ia,w(m2))
   50 continue
      call ea08cd(w,w(m1),value,vector,m,iv,w(m2))
      if (m.le.2) go to 110
      call ea06dd(vector,iv,m)
      do 100 i = m - 2,1,-1
        if (w(m1+i).ne.0) then
          aiw = one/ (a(i,i+1)*w(m1+i))
          do 60 ida = i + 1,m
            w(ida) = a(i,ida)
   60     continue
          do 90 l = 1,m
            pp = zero
            do 70 ida = i + 1,m
              pp = pp + w(ida)*vector(l,ida)
   70       continue
            pp = pp*aiw
            do 80 ida = i + 1,m
              vector(l,ida) = vector(l,ida) + pp*w(ida)
   80       continue
   90     continue
        end if
  100 continue
      call ea06dd(vector,iv,m)
  110 return
      end
      subroutine ea06dd(v,iv,m)
      integer iv,m
      double precision v(iv,m)
      double precision x
      integer i,j
      do 20 i = 1,m - 1
        do 10 j = i + 1,m
          x = v(i,j)
          v(i,j) = v(j,i)
          v(j,i) = x
   10   continue
   20 continue
      return
      end





* *******************************************************************
* copyright (c) 1970 hyprotech uk
* all rights reserved.
*
* none of the comments in this copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* this package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* all use is subject to licence. for full details of an hsl archive
* licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* please note that for an hsl archive licence:
*
* 1. the package must not be copied for use by any other person.
*    supply of any part of the library by the licensee to a third party
*    shall be subject to prior written agreement between aea
*    hyprotech uk limited and the licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. all information on the package is provided to the licensee on the
*    understanding that the details thereof are confidential.
* 3. all publications issued by the licensee that include results obtained
*    with the help of one or more of the packages shall acknowledge the
*    use of the packages. the licensee will notify the numerical analysis
*    group at rutherford appleton laboratory of any such publication.
* 4. the packages may be modified by or on behalf of the licensee
*    for such use in research applications but at no time shall such
*    packages or modifications thereof become the property of the
*    licensee. the licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. neither cclrc nor hyprotech uk limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of packages by the licensee.
* *******************************************************************
*
*######date 17 february 1994
c       toolpack tool decs employed.
c       modified for obsolescent features (feb 1997)
c
      subroutine ea07cd(a,value,m,ia,w)
      integer ia,m
      double precision a(ia,m),value(m),w(*)
      integer m2
      external ea09cd,mc04bd
      w(1) = a(1,1)
      m2 = m + m + 1
      if(m-2.lt.0) go to 60
      if(m-2.gt.0) go to 15
      w(2) = a(2,2)
      w(4) = a(2,1)
      go to 60
   15 call mc04bd(a,w(1),w(m+1),m,ia,w(m2))
   60 call ea09cd(w(1),w(m+1),value,m,w(m2))
      return
      end





* *******************************************************************
* copyright (c) 1970 hyprotech uk
* all rights reserved.
*
* none of the comments in this copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* this package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* all use is subject to licence. for full details of an hsl archive
* licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* please note that for an hsl archive licence:
*
* 1. the package must not be copied for use by any other person.
*    supply of any part of the library by the licensee to a third party
*    shall be subject to prior written agreement between aea
*    hyprotech uk limited and the licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. all information on the package is provided to the licensee on the
*    understanding that the details thereof are confidential.
* 3. all publications issued by the licensee that include results obtained
*    with the help of one or more of the packages shall acknowledge the
*    use of the packages. the licensee will notify the numerical analysis
*    group at rutherford appleton laboratory of any such publication.
* 4. the packages may be modified by or on behalf of the licensee
*    for such use in research applications but at no time shall such
*    packages or modifications thereof become the property of the
*    licensee. the licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. neither cclrc nor hyprotech uk limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of packages by the licensee.
*
* *******************************************************************
*
*######date 26 oct 1992
c       toolpack tool decs employed.
c       variables declared but not used have been removed.
c       parameter statement added and a34 set.
c       b(1) set.
c
      subroutine ea08cd(a,b,value,vec,m,iv,w)
      double precision one,zero
      parameter (one=1.0d0,zero=0.0d0)
      integer iv,m
      double precision a(m),b(m),value(m),vec(iv,m),w(2*m)
      double precision a11,a12,a13,a21,a22,a23,a33,a34,bb,cc,co,eps,
     +                 root,s,si,sim,sml,v1,xax,xx
      integer i,ii,iter,j,ji,k,mi,mn,mn2,n1,n2,n2m1,n3
      double precision fd05ad
      external fd05ad
      external ea09cd
      intrinsic abs,float,min,sign,sqrt
      a34 = zero
      eps = fd05ad(1)*10.0d0
      sml = eps*float(m)
      call ea09cd(a,b,w(m+1),m,w)
      b(1) = zero
      do 20 i = 1,m
        value(i) = a(i)
        w(i) = b(i)
        do 10 j = 1,m
          vec(j,i) = zero
   10   continue
        vec(i,i) = one
   20 continue
      k = 0
      if (m.eq.1) return
      do 90 n3 = 2,m
        n2 = m + 2 - n3
        mn2 = m + n2
        root = w(mn2)
        do 80 iter = 1,20
          bb = (value(n2)-value(n2-1))*0.5d0
          cc = w(n2)*w(n2)
          a22 = value(n2)
          if (cc.ne.0.0d0) a22 = a22 + cc/
     +                           (bb+sign(1.0d0,bb)*sqrt(bb*bb+cc))
          do 30 i = 1,n2
            mi = m + i
            if (abs(root-a22).le.abs(w(mi)-a22)) go to 30
            root = w(mi)
            mn = m + n2
            w(mi) = w(mn)
            w(mn) = root
   30     continue
          do 40 ii = 2,n2
            n1 = 2 + n2 - ii
            if (abs(w(n1)).le. (abs(value(n1-1))+abs(value(n1)))*
     +          sml) go to 50
   40     continue
          n1 = 1
   50     continue
          if (n2.eq.n1) go to 90
          n2m1 = n2 - 1
          if (iter.ge.3) root = a22
          k = k + 1
          a22 = value(n1)
          a12 = a22 - root
          a23 = w(n1+1)
          a13 = a23
          do 70 i = n1,n2m1
            a33 = value(i+1)
            if (i.ne.n2m1) a34 = w(i+2)
            s = sign(sqrt(a12*a12+a13*a13),a12)
            si = a13/s
            co = a12/s
            sim = -si
            do 60 ji = 1,min(m,i+k)
              v1 = vec(ji,i)
              vec(ji,i) = co*vec(ji,i) + si*vec(ji,i+1)
              vec(ji,i+1) = co*vec(ji,i+1) + sim*v1
   60       continue
            if (i.ne.n1) w(i) = s
            a11 = co*a22 + si*a23
            a12 = co*a23 + si*a33
            a13 = si*a34
            a21 = co*a23 - si*a22
            a22 = co*a33 - si*a23
            a23 = co*a34
            value(i) = a11*co + a12*si
            a12 = -a11*si + a12*co
            w(i+1) = a12
            a22 = a22*co - a21*si
   70     continue
          value(n2) = a22
   80   continue
        write (6,fmt=9000)
        stop
   90 continue
      do 110 j = 1,m
        xx = vec(1,j)**2
        xax = xx*a(1)
        do 100 i = 2,m
          xx = xx + vec(i,j)**2
          xax = xax + vec(i,j)* (2.0d0*b(i)*vec(i-1,j)+a(i)*vec(i,j))
  100   continue
        value(j) = xax/xx
  110 continue
 9000 format ('1cycle detected in subroutine ea08 -stopping now')
      return
      end





* *******************************************************************
* copyright (c) 1970 hyprotech uk.
* all rights reserved.
*
* none of the comments in this copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* this package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* all use is subject to licence. for full details of an hsl archive
* licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* please note that for an hsl archive licence:
*
* 1. the package must not be copied for use by any other person.
*    supply of any part of the library by the licensee to a third party
*    shall be subject to prior written agreement between aea
*    hyprotech uk limited and the licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. all information on the package is provided to the licensee on the
*    understanding that the details thereof are confidential.
* 3. all publications issued by the licensee that include results obtained
*    with the help of one or more of the packages shall acknowledge the
*    use of the packages. the licensee will notify the numerical analysis
*    group at rutherford appleton laboratory of any such publication.
* 4. the packages may be modified by or on behalf of the licensee
*    for such use in research applications but at no time shall such
*    packages or modifications thereof become the property of the
*    licensee. the licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. neither cclrc nor hyprotech uk limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of packages by the licensee.
* *******************************************************************
*
*######date   29 sept. 1992
c      13/10/92. modified to parameterize zero and set a34.
c      09/12/92. toolpack tool decs employed.
c      label 75 removed.
c      modified for obsolescent features (feb 1997)
c
      subroutine ea09cd(a,b,value,m,off)
!      real zero,one
      double precision zero,one
      parameter (zero=0.0d0,one=1.0d0)
      integer m
      double precision a(m),b(m),off(m),value(m)
      double precision a11,a12,a13,a21,a22,a23,a33,a34,bb,cc,co,eps,
     +                 root,s,sbb,si,sml,xmax
      integer i,ii,iter,maxit,n1,n2,n2m1,n3
      double precision fd05ad
      external fd05ad
      intrinsic dabs,dmax1,dsign,dsqrt,float
      eps = fd05ad(1)*10.0d0
      sml = eps*float(m)
      value(1) = a(1)
      if (m.eq.1) return
      do 10 i = 2,m
        value(i) = a(i)
        off(i) = b(i)
   10 continue
      maxit = 10*m
      do 90 iter = 1,maxit
        do 45 n3 = 2,m
          n2 = m + 2 - n3
          do 30 ii = 2,n2
            n1 = 2 + n2 - ii
            if (dabs(off(n1)).le. (dabs(value(n1-1))+dabs(value(n1)))*
     +          sml) go to 40
   30     continue
          n1 = 1
   40     if (n2.ne.n1) go to 50
   45   continue
        return
   50   bb = (value(n2)-value(n2-1))*0.5d0
        cc = off(n2)*off(n2)
        sbb = 1.0d0
        if (bb.lt.zero) sbb = -1.0d0
        root = value(n2) + cc/ (bb+sbb*dsqrt(bb*bb+cc))
        n2m1 = n2 - 1
        a22 = value(n1)
        a12 = a22 - root
        a23 = off(n1+1)
        a13 = a23
        a34 = zero
        do 80 i = n1,n2m1
          a33 = value(i+1)
          if (i.ne.n2m1) a34 = off(i+2)
          xmax = dmax1(dabs(a12),dabs(a13))
          s = dsign(dsqrt((a12/xmax)**2+ (a13/xmax)**2)*xmax,a12)
          si = a13/s
          co = a12/s
          if (i.ne.n1) off(i) = s
          a11 = co*a22 + si*a23
          a12 = co*a23 + si*a33
          a13 = si*a34
          a21 = co*a23 - si*a22
          a22 = co*a33 - si*a23
          a23 = co*a34
          value(i) = a11*co + a12*si
          a12 = -a11*si + a12*co
          off(i+1) = a12
          a22 = a22*co - a21*si
   80   continue
        value(n2) = a22
   90 continue
      write (6,fmt=100)
  100 format ( '1looping detected in ea09-stopping now ')
      stop
      end





* *******************************************************************
* copyright (c) 1966 hyprotech uk
* all rights reserved.
*
* none of the comments in this copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* this package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* all use is subject to licence. for full details of an hsl archive
* licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* please note that for an hsl archive licence:
*
* 1. the package must not be copied for use by any other person.
*    supply of any part of the library by the licensee to a third party
*    shall be subject to prior written agreement between aea
*    hyprotech uk limited and the licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. all information on the package is provided to the licensee on the
*    understanding that the details thereof are confidential.
* 3. all publications issued by the licensee that include results obtained
*    with the help of one or more of the packages shall acknowledge the
*    use of the packages. the licensee will notify the numerical analysis
*    group at rutherford appleton laboratory of any such publication.
* 4. the packages may be modified by or on behalf of the licensee
*    for such use in research applications but at no time shall such
*    packages or modifications thereof become the property of the
*    licensee. the licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. neither cclrc nor hyprotech uk limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of packages by the licensee.
* *******************************************************************
*
*######date 17 february 1994
c       toolpack tool decs employed.
c       data made parameter.
c
      subroutine mc04bd(a,alpha,beta,m,ia,q)
      integer ia,m
      double precision a(ia,m),alpha(m),beta(m),q(m)
      double precision aik,beta1,h,hrec,one,p5,pp,pp1,qk,qw1,qw2,zero
      integer i,j,j1,k,ki,kj,l
      intrinsic dsqrt
      parameter (zero=0.0d0,one=1.0d0,p5=0.5d0)
      alpha(1) = a(1,1)
      do 21 j = 2,m
        j1 = j - 1
        do 22 i = 1,j1
          a(i,j) = a(j,i)
   22   continue
        alpha(j) = a(j,j)
   21 continue
      do 1 i = 1,m - 2
        pp = zero
        do 2 j = i + 1,m
          pp = pp + a(i,j)**2
          beta(j) = a(i,j)
    2   continue
        pp1 = dsqrt(pp)
        if (a(i,i+1).ge.zero) then
          beta1 = -pp1
        else
          beta1 = pp1
        end if
        if (pp.gt.zero) then
          h = pp - beta1*a(i,i+1)
          a(i,i+1) = a(i,i+1) - beta1
          beta(i+1) = a(i,i+1)
          hrec = one/h
          do 207 k = i + 1,m
            qk = zero
            do 208 l = i + 1,k
              qk = qk + beta(l)*a(l,k)
  208       continue
            q(k) = qk
  207     continue
          do 7 k = i + 2,m
            aik = beta(k)
            do 8 l = i + 1,k - 1
              q(l) = q(l) + aik*a(l,k)
    8       continue
    7     continue
          pp = zero
          do 118 k = i + 1,m
            q(k) = hrec*q(k)
            pp = pp + beta(k)*q(k)
  118     continue
          pp = -pp*p5*hrec
          do 11 ki = i + 1,m
            q(ki) = q(ki) + pp*beta(ki)
            qw1 = -q(ki)
            qw2 = -beta(ki)
            do 12 kj = i + 1,ki
              a(kj,ki) = a(kj,ki) + qw1*beta(kj) + qw2*q(kj)
   12       continue
   11     continue
        end if
        beta(i+1) = beta1
    1 continue
      do 23 i = 2,m
        h = alpha(i)
        alpha(i) = a(i,i)
        a(i,i) = h
   23 continue
      beta(m) = a(m-1,m)
      return
      end





* *******************************************************************
* copyright (c) 1988 aea technology
*######date 21 jan 1993
c       toolpack tool decs employed.
c       save statement added.
c 1/10/98 dc(3) not initialized to avoid sun f90 failure
c 16 october 2001: stop and write statements removed.
* *******************************************************************
      double precision function fd05ad(inum)
c----------------------------------------------------------------
c  real constants for: ieee double precision (8-byte arithmetic)
c
c  obtained from h.s.l. subroutine ze02am.
c  nick gould and sid marlow, harwell laboratory, april 1988.
c----------------------------------------------------------------
c     .. scalar arguments ..
      integer inum
c     ..
c     .. local arrays ..
      double precision dc(5)
c     ..
c     .. save statement ..
      save dc
c     ..
c     .. data statements ..
c
c  dc(1) the smallest positive number: 1.0 + dc(1) > 1.0.
c  dc(2) the smallest positive number: 1.0 - dc(2) < 1.0.
c  dc(3) the smallest nonzero +ve real number.
c  dc(4) the smallest full precision +ve real number.
c  dc(5) the largest finite +ve real number.
c
      data dc(1)/2.2204460492504d-16/
      data dc(2)/1.1102230246253d-16/
c     data dc(3)/4.9406564584126d-324/
      data dc(4)/2.2250738585073d-308/
      data dc(5)/1.7976931348622d+308/
c     ..
c     .. executable statements ..

      if ( inum .le. 0 ) then
         fd05ad = dc( 1 )
      else if ( inum .ge. 6 ) then
         fd05ad = dc( 5 )
      else if ( inum .eq. 3 ) then
         fd05ad = dc(4)/2.0d0**52
      else
         fd05ad = dc( inum )
      endif
      return
      end