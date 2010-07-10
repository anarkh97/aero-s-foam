! ==================================
! geometry: tripack
! ==================================
!
! ref:
! ---
! Algorithm 751: TRIPACK: a constrained two-dimensional Delaunay triangulation package
! R. J. Renka, ACM Transactions on Mathematical Software (TOMS) archive, Vol. 22, No. 1, 1996
!
! =========================================================================================================



      subroutine addcst (ncc,lcc,n,x,y, lwk,iwk,list,lptr,
     .                   lend, ier)
      integer ncc, lcc(*), n, lwk, iwk(lwk), list(*),
     .        lptr(*), lend(n), ier
      real(8) x(n), y(n)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   11/12/94
c
c   this subroutine provides for creation of a constrained
c delaunay triangulation which, in some sense, covers an
c arbitrary connected region r rather than the convex hull
c of the nodes.  this is achieved simply by forcing the
c presence of certain adjacencies (triangulation arcs) cor-
c responding to constraint curves.  the union of triangles
c coincides with the convex hull of the nodes, but triangles
c in r can be distinguished from those outside of r.  the
c only modification required to generalize the definition of
c the delaunay triangulation is replacement of property 5
c (refer to trmesh) by the following:
c
c  5')  if a node is contained in the interior of the cir-
c       cumcircle of a triangle, then every interior point
c       of the triangle is separated from the node by a
c       constraint arc.
c
c   in order to be explicit, we make the following defini-
c tions.  a constraint region is the open interior of a
c simple closed positively oriented polygonal curve defined
c by an ordered sequence of three or more distinct nodes
c (constraint nodes) p(1),p(2),...,p(k), such that p(i) is
c adjacent to p(i+1) for i = 1,...,k with p(k+1) = p(1).
c thus, the constraint region is on the left (and may have
c nonfinite area) as the sequence of constraint nodes is
c traversed in the specified order.  the constraint regions
c must not contain nodes and must not overlap.  the region
c r is the convex hull of the nodes with constraint regions
c excluded.
c
c   note that the terms boundary node and boundary arc are
c reserved for nodes and arcs on the boundary of the convex
c hull of the nodes.
c
c   the algorithm is as follows:  given a triangulation
c which includes one or more sets of constraint nodes, the
c corresponding adjacencies (constraint arcs) are forced to
c be present (subroutine edge).  any additional new arcs
c required are chosen to be locally optimal (satisfy the
c modified circumcircle property).
c
c
c on input:
c
c       ncc = number of constraint curves (constraint re-
c             gions).  ncc .ge. 0.
c
c       lcc = array of length ncc (or dummy array of length
c             1 if ncc = 0) containing the index (for x, y,
c             and lend) of the first node of constraint i in
c             lcc(i) for i = 1 to ncc.  thus, constraint i
c             contains k = lcc(i+1) - lcc(i) nodes, k .ge.
c             3, stored in (x,y) locations lcc(i), ...,
c             lcc(i+1)-1, where lcc(ncc+1) = n+1.
c
c       n = number of nodes in the triangulation, including
c           constraint nodes.  n .ge. 3.
c
c       x,y = arrays of length n containing the coordinates
c             of the nodes with non-constraint nodes in the
c             first lcc(1)-1 locations, followed by ncc se-
c             quences of constraint nodes.  only one of
c             these sequences may be specified in clockwise
c             order to represent an exterior constraint
c             curve (a constraint region with nonfinite
c             area).
c
c the above parameters are not altered by this routine.
c
c       lwk = length of iwk.  this must be at least 2*ni
c             where ni is the maximum number of arcs which
c             intersect a constraint arc to be added.  ni
c             is bounded by n-3.
c
c       iwk = integer work array of length lwk (used by
c             subroutine edge to add constraint arcs).
c
c       list,lptr,lend = data structure defining the trian-
c                        gulation.  refer to subroutine
c                        trmesh.
c
c on output:
c
c       lwk = required length of iwk unless ier = 1 or ier =
c             3.  in the case of ier = 1, lwk is not altered
c             from its input value.
c
c       iwk = array containing the endpoint indexes of the
c             new arcs which were swapped in by the last
c             call to subroutine edge.
c
c       list,lptr,lend = triangulation data structure with
c                        all constraint arcs present unless
c                        ier .ne. 0.  these arrays are not
c                        altered if ier = 1.
c
c       ier = error indicator:
c             ier = 0 if no errors were encountered.
c             ier = 1 if ncc, n, or an lcc entry is outside
c                     its valid range, or lwk .lt. 0 on
c                     input.
c             ier = 2 if more space is required in iwk.
c             ier = 3 if the triangulation data structure is
c                     invalid, or failure (in edge or optim)
c                     was caused by collinear nodes on the
c                     convex hull boundary.  an error mes-
c                     sage is written to logical unit 6 in
c                     this case.
c             ier = 4 if intersecting constraint arcs were
c                     encountered.
c             ier = 5 if a constraint region contains a
c                     node.
c
c modules required by addcst:  edge, left, lstptr, optim,
c                                swap, swptst
c
c intrinsic functions called by addcst:  abs, max
c
c***********************************************************
c
      integer i, ifrst, ilast, k, kbak, kfor, kn, lccip1,
     .        lp, lpb, lpf, lpl, lw, lwd2, n1, n2
      lwd2 = lwk/2
c
c  test for errors in input parameters.
c
      ier = 1
      if (ncc .lt. 0  .or.  lwk .lt. 0) return
      if (ncc .eq. 0) then
        if (n .lt. 3) return
        lwk = 0
        go to 9
      else
        lccip1 = n+1
        do 1 i = ncc,1,-1
          if (lccip1 - lcc(i) .lt. 3) return
          lccip1 = lcc(i)
    1     continue
        if (lccip1 .lt. 1) return
      endif
c
c  force the presence of constraint arcs.  the outer loop is
c    on constraints in reverse order.  ifrst and ilast are
c    the first and last nodes of constraint i.
c
      lwk = 0
      ifrst = n+1
      do 3 i = ncc,1,-1
        ilast = ifrst - 1
        ifrst = lcc(i)
c
c   inner loop on constraint arcs n1-n2 in constraint i.
c
        n1 = ilast
        do 2 n2 = ifrst,ilast
          lw = lwd2
          call edge (n1,n2,x,y, lw,iwk,list,lptr,lend, ier)
          lwk = max(lwk,2*lw)
          if (ier .eq. 4) ier = 3
          if (ier .ne. 0) return
          n1 = n2
    2     continue
    3   continue
c
c test for errors.  the outer loop is on constraint i with
c   first and last nodes ifrst and ilast, and the inner loop
c   is on constraint nodes k with (kbak,k,kfor) a subse-
c   quence of constraint i.
c
      ier = 4
      ifrst = n+1
      do 8 i = ncc,1,-1
        ilast = ifrst - 1
        ifrst = lcc(i)
        kbak = ilast
        do 7 k = ifrst,ilast
          kfor = k + 1
          if (k .eq. ilast) kfor = ifrst
c
c   find the list pointers lpf and lpb of kfor and kbak as
c     neighbors of k.
c
          lpf = 0
          lpb = 0
          lpl = lend(k)
          lp = lpl
c
    4     lp = lptr(lp)
            kn = abs(list(lp))
            if (kn .eq. kfor) lpf = lp
            if (kn .eq. kbak) lpb = lp
            if (lp .ne. lpl) go to 4
c
c   a pair of intersecting constraint arcs was encountered
c     if and only if a constraint arc is missing (introduc-
c     tion of the second caused the first to be swapped out).
c
          if (lpf .eq. 0  .or.  lpb .eq. 0) return
c
c   loop on neighbors kn of node k which follow kfor and
c     precede kbak.  the constraint region contains no nodes
c     if and only if all such nodes kn are in constraint i.
c
          lp = lpf
    5     lp = lptr(lp)
            if (lp .eq. lpb) go to 6
            kn = abs(list(lp))
            if (kn .lt. ifrst  .or.  kn .gt. ilast) go to 10
            go to 5
c
c   bottom of loop.
c
    6     kbak = k
    7     continue
    8   continue
c
c no errors encountered.
c
    9 ier = 0
      return
c
c a constraint region contains a node.
c
   10 ier = 5
      return
      end
      subroutine addnod (k,xk,yk,ist,ncc, lcc,n,x,y,list,
     .                   lptr,lend,lnew, ier)
      integer k, ist, ncc, lcc(*), n, list(*), lptr(*),
     .        lend(*), lnew, ier
      real(8)    xk, yk, x(*), y(*)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   06/27/98
c
c   given a triangulation of n nodes in the plane created by
c subroutine trmesh or trmshr, this subroutine updates the
c data structure with the addition of a new node in position
c k.  if node k is inserted into x and y (k .le. n) rather
c than appended (k = n+1), then a corresponding insertion
c must be performed in any additional arrays associated
c with the nodes.  for example, an array of data values z
c must be shifted down to open up position k for the new
c value:  set z(i+1) to z(i) for i = n,n-1,...,k.  for
c optimal efficiency, new nodes should be appended whenever
c possible.  insertion is necessary, however, to add a non-
c constraint node when constraints are present (refer to
c subroutine addcst).
c
c   note that a constraint node cannot be added by this
c routine.  in order to insert a constraint node, it is
c necessary to add the node with no constraints present
c (call this routine with ncc = 0), update lcc by increment-
c ing the appropriate entries, and then create (or restore)
c the constraints by a call to addcst.
c
c   the algorithm consists of the following steps:  node k
c is located relative to the triangulation (trfind), its
c index is added to the data structure (intadd or bdyadd),
c and a sequence of swaps (swptst and swap) are applied to
c the arcs opposite k so that all arcs incident on node k
c and opposite node k (excluding constraint arcs) are local-
c ly optimal (satisfy the circumcircle test).  thus, if a
c (constrained) delaunay triangulation is input, a (con-
c strained) delaunay triangulation will result.  all indexes
c are incremented as necessary for an insertion.
c
c
c on input:
c
c       k = nodal index (index for x, y, and lend) of the
c           new node to be added.  1 .le. k .le. lcc(1).
c           (k .le. n+1 if ncc=0).
c
c       xk,yk = cartesian coordinates of the new node (to be
c               stored in x(k) and y(k)).  the node must not
c               lie in a constraint region.
c
c       ist = index of a node at which trfind begins the
c             search.  search time depends on the proximity
c             of this node to node k.  1 .le. ist .le. n.
c
c       ncc = number of constraint curves.  ncc .ge. 0.
c
c the above parameters are not altered by this routine.
c
c       lcc = list of constraint curve starting indexes (or
c             dummy array of length 1 if ncc = 0).  refer to
c             subroutine addcst.
c
c       n = number of nodes in the triangulation before k is
c           added.  n .ge. 3.  note that n will be incre-
c           mented following the addition of node k.
c
c       x,y = arrays of length at least n+1 containing the
c             cartesian coordinates of the nodes in the
c             first n positions with non-constraint nodes
c             in the first lcc(1)-1 locations if ncc > 0.
c
c       list,lptr,lend,lnew = data structure associated with
c                             the triangulation of nodes 1
c                             to n.  the arrays must have
c                             sufficient length for n+1
c                             nodes.  refer to trmesh.
c
c on output:
c
c       lcc = list of constraint curve starting indexes in-
c             cremented by 1 to reflect the insertion of k
c             unless ncc = 0 or (ier .ne. 0 and ier .ne.
c             -4).
c
c       n = number of nodes in the triangulation including k
c           unless ier .ne. 0 and ier .ne. -4.  note that
c           all comments refer to the input value of n.
c
c       x,y = arrays updated with the insertion of xk and yk
c             in the k-th positions (node i+1 was node i be-
c             fore the insertion for i = k to n if k .le. n)
c             unless ier .ne. 0 and ier .ne. -4.
c
c       list,lptr,lend,lnew = data structure updated with
c                             the addition of node k unless
c                             ier .ne. 0 and ier .ne. -4.
c
c       ier = error indicator:
c             ier =  0 if no errors were encountered.
c             ier = -1 if k, ist, ncc, n, or an lcc entry is
c                      outside its valid range on input.
c             ier = -2 if all nodes (including k) are col-
c                      linear.
c             ier =  l if nodes l and k coincide for some l.
c             ier = -3 if k lies in a constraint region.
c             ier = -4 if an error flag is returned by swap
c                      implying that the triangulation
c                      (geometry) was bad on input.
c
c             the errors conditions are tested in the order
c             specified.
c
c modules required by addnod:  bdyadd, crtri, indxcc,
c                                insert, intadd, jrand,
c                                left, lstptr, swap,
c                                swptst, trfind
c
c intrinsic function called by addnod:  abs
c
c***********************************************************
c
      integer indxcc, lstptr
      integer i, i1, i2, i3, ibk, io1, io2, in1, kk, l,
     .        lccip1, lp, lpf, lpo1, nm1
      logical crtri, swptst
      kk = k
c
c test for an invalid input parameter.
c
      if (kk .lt. 1  .or.  ist .lt. 1  .or.  ist .gt. n
     .    .or.  ncc .lt. 0  .or.  n .lt. 3) go to 7
      lccip1 = n+1
      do 1 i = ncc,1,-1
        if (lccip1-lcc(i) .lt. 3) go to 7
        lccip1 = lcc(i)
    1   continue
      if (kk .gt. lccip1) go to 7
c
c find a triangle (i1,i2,i3) containing k or the rightmost
c   (i1) and leftmost (i2) visible boundary nodes as viewed
c   from node k.
c
      call trfind (ist,xk,yk,n,x,y,list,lptr,lend, i1,i2,i3)
c
c test for collinear nodes, duplicate nodes, and k lying in
c   a constraint region.
c
      if (i1 .eq. 0) go to 8
      if (i3 .ne. 0) then
        l = i1
        if (xk .eq. x(l)  .and.  yk .eq. y(l)) go to 9
        l = i2
        if (xk .eq. x(l)  .and.  yk .eq. y(l)) go to 9
        l = i3
        if (xk .eq. x(l)  .and.  yk .eq. y(l)) go to 9
        if (ncc .gt. 0  .and.  crtri(ncc,lcc,i1,i2,i3) )
     .    go to 10
      else
c
c   k is outside the convex hull of the nodes and lies in a
c     constraint region iff an exterior constraint curve is
c     present.
c
        if (ncc .gt. 0  .and.  indxcc(ncc,lcc,n,list,lend)
     .      .ne. 0) go to 10
      endif
c
c no errors encountered.
c
      ier = 0
      nm1 = n
      n = n + 1
      if (kk .lt. n) then
c
c open a slot for k in x, y, and lend, and increment all
c   nodal indexes which are greater than or equal to k.
c   note that list, lptr, and lnew are not yet updated with
c   either the neighbors of k or the edges terminating on k.
c
        do 2 ibk = nm1,kk,-1
          x(ibk+1) = x(ibk)
          y(ibk+1) = y(ibk)
          lend(ibk+1) = lend(ibk)
    2     continue
        do 3 i = 1,ncc
          lcc(i) = lcc(i) + 1
    3     continue
        l = lnew - 1
        do 4 i = 1,l
          if (list(i) .ge. kk) list(i) = list(i) + 1
          if (list(i) .le. -kk) list(i) = list(i) - 1
    4     continue
        if (i1 .ge. kk) i1 = i1 + 1
        if (i2 .ge. kk) i2 = i2 + 1
        if (i3 .ge. kk) i3 = i3 + 1
      endif
c
c insert k into x and y, and update list, lptr, lend, and
c   lnew with the arcs containing node k.
c
      x(kk) = xk
      y(kk) = yk
      if (i3 .eq. 0) then
        call bdyadd (kk,i1,i2, list,lptr,lend,lnew )
      else
        call intadd (kk,i1,i2,i3, list,lptr,lend,lnew )
      endif
c
c initialize variables for optimization of the triangula-
c   tion.
c
      lp = lend(kk)
      lpf = lptr(lp)
      io2 = list(lpf)
      lpo1 = lptr(lpf)
      io1 = abs(list(lpo1))
c
c begin loop:  find the node opposite k.
c
    5 lp = lstptr(lend(io1),io2,list,lptr)
        if (list(lp) .lt. 0) go to 6
        lp = lptr(lp)
        in1 = abs(list(lp))
        if ( crtri(ncc,lcc,io1,io2,in1) ) go to 6
c
c swap test:  if a swap occurs, two new arcs are
c             opposite k and must be tested.
c
        if ( .not. swptst(in1,kk,io1,io2,x,y) ) go to 6
        call swap (in1,kk,io1,io2, list,lptr,lend, lpo1)
        if (lpo1 .eq. 0) go to 11
        io1 = in1
        go to 5
c
c no swap occurred.  test for termination and reset
c   io2 and io1.
c
    6   if (lpo1 .eq. lpf  .or.  list(lpo1) .lt. 0) return
        io2 = io1
        lpo1 = lptr(lpo1)
        io1 = abs(list(lpo1))
        go to 5
c
c a parameter is outside its valid range on input.
c
    7 ier = -1
      return
c
c all nodes are collinear.
c
    8 ier = -2
      return
c
c nodes l and k coincide.
c
    9 ier = l
      return
c
c node k lies in a constraint region.
c
   10 ier = -3
      return
c
c zero pointer returned by swap.
c
   11 ier = -4
      return
      end
      real(8) function areap (x,y,nb,nodes)
      integer nb, nodes(nb)
      real(8)    x(*), y(*)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   09/21/90
c
c   given a sequence of nb points in the plane, this func-
c tion computes the signed area bounded by the closed poly-
c gonal curve which passes through the points in the
c specified order.  each simple closed curve is positively
c oriented (bounds positive area) if and only if the points
c are specified in counterclockwise order.  the last point
c of the curve is taken to be the first point specified, and
c this point should therefore not be specified twice.
c
c   the area of a triangulation may be computed by calling
c areap with values of nb and nodes determined by subroutine
c bnodes.
c
c
c on input:
c
c       x,y = arrays of length n containing the cartesian
c             coordinates of a set of points in the plane
c             for some n .ge. nb.
c
c       nb = length of nodes.
c
c       nodes = array of length nb containing the ordered
c               sequence of nodal indexes (in the range
c               1 to n) which define the polygonal curve.
c
c input parameters are not altered by this function.
c
c on output:
c
c       areap = signed area bounded by the polygonal curve,
c              or zero if nb < 3.
c
c modules required by areap:  none
c
c***********************************************************
c
      integer i, nd1, nd2, nnb
      real(8)    a
c
c local parameters:
c
c a =       partial sum of signed (and doubled) trapezoid
c             areas
c i =       do-loop and nodes index
c nd1,nd2 = elements of nodes
c nnb =     local copy of nb
c
      nnb = nb
      a = 0.
      if (nnb .lt. 3) go to 2
      nd2 = nodes(nnb)
c
c loop on line segments nodes(i-1) -> nodes(i), where
c   nodes(0) = nodes(nb), adding twice the signed trapezoid
c   areas (integrals of the linear interpolants) to a.
c
      do 1 i = 1,nnb
        nd1 = nd2
        nd2 = nodes(i)
        a = a + (x(nd2)-x(nd1))*(y(nd1)+y(nd2))
    1   continue
c
c a contains twice the negative signed area of the region.
c
    2 areap = -a/2.
      return
      end
      subroutine bdyadd (kk,i1,i2, list,lptr,lend,lnew )
      integer kk, i1, i2, list(*), lptr(*), lend(*), lnew
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   02/22/91
c
c   this subroutine adds a boundary node to a triangulation
c of a set of points in the plane.  the data structure is
c updated with the insertion of node kk, but no optimization
c is performed.
c
c
c on input:
c
c       kk = index of a node to be connected to the sequence
c            of all visible boundary nodes.  kk .ge. 1 and
c            kk must not be equal to i1 or i2.
c
c       i1 = first (rightmost as viewed from kk) boundary
c            node in the triangulation which is visible from
c            node kk (the line segment kk-i1 intersects no
c            arcs.
c
c       i2 = last (leftmost) boundary node which is visible
c            from node kk.  i1 and i2 may be determined by
c            subroutine trfind.
c
c the above parameters are not altered by this routine.
c
c       list,lptr,lend,lnew = triangulation data structure
c                             created by trmesh or trmshr.
c                             nodes i1 and i2 must be in-
c                             cluded in the triangulation.
c
c on output:
c
c       list,lptr,lend,lnew = data structure updated with
c                             the addition of node kk.  node
c                             kk is connected to i1, i2, and
c                             all boundary nodes in between.
c
c module required by bdyadd:  insert
c
c***********************************************************
c
      integer k, lp, lsav, n1, n2, next, nsav
      k = kk
      n1 = i1
      n2 = i2
c
c add k as the last neighbor of n1.
c
      lp = lend(n1)
      lsav = lptr(lp)
      lptr(lp) = lnew
      list(lnew) = -k
      lptr(lnew) = lsav
      lend(n1) = lnew
      lnew = lnew + 1
      next = -list(lp)
      list(lp) = next
      nsav = next
c
c loop on the remaining boundary nodes between n1 and n2,
c   adding k as the first neighbor.
c
    1 lp = lend(next)
        call insert (k,lp,list,lptr,lnew)
        if (next .eq. n2) go to 2
        next = -list(lp)
        list(lp) = next
        go to 1
c
c add the boundary nodes between n1 and n2 as neighbors
c   of node k.
c
    2 lsav = lnew
      list(lnew) = n1
      lptr(lnew) = lnew + 1
      lnew = lnew + 1
      next = nsav
c
    3 if (next .eq. n2) go to 4
        list(lnew) = next
        lptr(lnew) = lnew + 1
        lnew = lnew + 1
        lp = lend(next)
        next = list(lp)
        go to 3
c
    4 list(lnew) = -n2
      lptr(lnew) = lsav
      lend(k) = lnew
      lnew = lnew + 1
      return
      end
      subroutine bnodes (n,list,lptr,lend, nodes,nb,na,nt)
      integer n, list(*), lptr(*), lend(n), nodes(*), nb,
     .        na, nt
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   09/01/88
c
c   given a triangulation of n points in the plane, this
c subroutine returns an array containing the indexes, in
c counterclockwise order, of the nodes on the boundary of
c the convex hull of the set of points.
c
c
c on input:
c
c       n = number of nodes in the triangulation.  n .ge. 3.
c
c       list,lptr,lend = data structure defining the trian-
c                        gulation.  refer to subroutine
c                        trmesh.
c
c the above parameters are not altered by this routine.
c
c       nodes = integer array of length at least nb
c               (nb .le. n).
c
c on output:
c
c       nodes = ordered sequence of boundary node indexes
c               in the range 1 to n.
c
c       nb = number of boundary nodes.
c
c       na,nt = number of arcs and triangles, respectively,
c               in the triangulation.
c
c modules required by bnodes:  none
c
c***********************************************************
c
      integer k, lp, n0, nst
c
c set nst to the first boundary node encountered.
c
      nst = 1
    1 lp = lend(nst)
        if (list(lp) .lt. 0) go to 2
        nst = nst + 1
        go to 1
c
c initialization.
c
    2 nodes(1) = nst
      k = 1
      n0 = nst
c
c traverse the boundary in counterclockwise order.
c
    3 lp = lend(n0)
        lp = lptr(lp)
        n0 = list(lp)
        if (n0 .eq. nst) go to 4
        k = k + 1
        nodes(k) = n0
        go to 3
c
c termination.
c
    4 nb = k
      nt = 2*n - nb - 2
      na = nt + n - 1
      return
      end
      subroutine circum (x1,y1,x2,y2,x3,y3,ratio, xc,yc,cr,
     .                   sa,ar)
      logical ratio
      real(8)    x1, y1, x2, y2, x3, y3, xc, yc, cr, sa, ar
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   12/10/96
c
c   given three vertices defining a triangle, this subrou-
c tine returns the circumcenter, circumradius, signed
c triangle area, and, optionally, the aspect ratio of the
c triangle.
c
c
c on input:
c
c       x1,...,y3 = cartesian coordinates of the vertices.
c
c       ratio = logical variable with value true if and only
c               if the aspect ratio is to be computed.
c
c input parameters are not altered by this routine.
c
c on output:
c
c       xc,yc = cartesian coordinates of the circumcenter
c               (center of the circle defined by the three
c               points) unless sa = 0, in which xc and yc
c               are not altered.
c
c       cr = circumradius (radius of the circle defined by
c            the three points) unless sa = 0 (infinite
c            radius), in which case cr is not altered.
c
c       sa = signed triangle area with positive value if
c            and only if the vertices are specified in
c            counterclockwise order:  (x3,y3) is strictly
c            to the left of the directed line from (x1,y1)
c            toward (x2,y2).
c
c       ar = aspect ratio r/cr, where r is the radius of the
c            inscribed circle, unless ratio = false, in
c            which case ar is not altered.  ar is in the
c            range 0 to .5, with value 0 iff sa = 0 and
c            value .5 iff the vertices define an equilateral
c            triangle.
c
c modules required by circum:  none
c
c intrinsic functions called by circum:  abs, sqrt
c
c***********************************************************
c
      integer i
      real(8)    ds(3), fx, fy, u(3), v(3)
c
c set u(k) and v(k) to the x and y components, respectively,
c   of the directed edge opposite vertex k.
c
      u(1) = x3 - x2
      u(2) = x1 - x3
      u(3) = x2 - x1
      v(1) = y3 - y2
      v(2) = y1 - y3
      v(3) = y2 - y1
c
c set sa to the signed triangle area.
c
      sa = (u(1)*v(2) - u(2)*v(1))/2.
      if (sa .eq. 0.) then
        if (ratio) ar = 0.
        return
      endif
c
c set ds(k) to the squared distance from the origin to
c   vertex k.
c
      ds(1) = x1*x1 + y1*y1
      ds(2) = x2*x2 + y2*y2
      ds(3) = x3*x3 + y3*y3
c
c compute factors of xc and yc.
c
      fx = 0.
      fy = 0.
      do 1 i = 1,3
        fx = fx - ds(i)*v(i)
        fy = fy + ds(i)*u(i)
    1   continue
      xc = fx/(4.*sa)
      yc = fy/(4.*sa)
      cr = sqrt( (xc-x1)**2 + (yc-y1)**2 )
      if (.not. ratio) return
c
c compute the squared edge lengths and aspect ratio.
c
      do 2 i = 1,3
        ds(i) = u(i)*u(i) + v(i)*v(i)
    2   continue
      ar = 2.*abs(sa)/
     .     ( (sqrt(ds(1)) + sqrt(ds(2)) + sqrt(ds(3)))*cr )
      return
      end
      logical function crtri (ncc,lcc,i1,i2,i3)
      integer ncc, lcc(*), i1, i2, i3
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   08/14/91
c
c   this function returns true if and only if triangle (i1,
c i2,i3) lies in a constraint region.
c
c
c on input:
c
c       ncc,lcc = constraint data structure.  refer to sub-
c                 routine addcst.
c
c       i1,i2,i3 = nodal indexes of the counterclockwise-
c                  ordered vertices of a triangle.
c
c input parameters are altered by this function.
c
c       crtri = true iff (i1,i2,i3) is a constraint region
c               triangle.
c
c note that input parameters are not tested for validity.
c
c modules required by crtri:  none
c
c intrinsic functions called by crtri:  max, min
c
c***********************************************************
c
      integer i, imax, imin
      imax = max(i1,i2,i3)
c
c   find the index i of the constraint containing imax.
c
      i = ncc + 1
    1 i = i - 1
        if (i .le. 0) go to 2
        if (imax .lt. lcc(i)) go to 1
      imin = min(i1,i2,i3)
c
c p lies in a constraint region iff i1, i2, and i3 are nodes
c   of the same constraint (imin >= lcc(i)), and (imin,imax)
c   is (i1,i3), (i2,i1), or (i3,i2).
c
      crtri = imin .ge. lcc(i)  .and.  ((imin .eq. i1 .and.
     .        imax .eq. i3)  .or.  (imin .eq. i2  .and.
     .        imax .eq. i1)  .or.  (imin .eq. i3  .and.
     .        imax .eq. i2))
      return
c
c ncc .le. 0 or all vertices are non-constraint nodes.
c
    2 crtri = .false.
      return
      end
      subroutine delarc (n,io1,io2, list,lptr,lend,
     .                   lnew, ier)
      integer n, io1, io2, list(*), lptr(*), lend(n), lnew,
     .        ier
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   11/12/94
c
c   this subroutine deletes a boundary arc from a triangula-
c tion.  it may be used to remove a null triangle from the
c convex hull boundary.  note, however, that if the union of
c triangles is rendered nonconvex, subroutines delnod, edge,
c and trfind may fail.  thus, subroutines addcst, addnod,
c delnod, edge, and nearnd should not be called following
c an arc deletion.
c
c
c on input:
c
c       n = number of nodes in the triangulation.  n .ge. 4.
c
c       io1,io2 = indexes (in the range 1 to n) of a pair of
c                 adjacent boundary nodes defining the arc
c                 to be removed.
c
c the above parameters are not altered by this routine.
c
c       list,lptr,lend,lnew = triangulation data structure
c                             created by trmesh or trmshr.
c
c on output:
c
c       list,lptr,lend,lnew = data structure updated with
c                             the removal of arc io1-io2
c                             unless ier > 0.
c
c       ier = error indicator:
c             ier = 0 if no errors were encountered.
c             ier = 1 if n, io1, or io2 is outside its valid
c                     range, or io1 = io2.
c             ier = 2 if io1-io2 is not a boundary arc.
c             ier = 3 if the node opposite io1-io2 is al-
c                     ready a boundary node, and thus io1
c                     or io2 has only two neighbors or a
c                     deletion would result in two triangu-
c                     lations sharing a single node.
c             ier = 4 if one of the nodes is a neighbor of
c                     the other, but not vice versa, imply-
c                     ing an invalid triangulation data
c                     structure.
c
c modules required by delarc:  delnb, lstptr
c
c intrinsic function called by delarc:  abs
c
c***********************************************************
c
      integer lstptr
      integer lp, lph, lpl, n1, n2, n3
      n1 = io1
      n2 = io2
c
c test for errors, and set n1->n2 to the directed boundary
c   edge associated with io1-io2:  (n1,n2,n3) is a triangle
c   for some n3.
c
      if (n .lt. 4  .or.  n1 .lt. 1  .or.  n1 .gt. n  .or.
     .    n2 .lt. 1  .or.  n2 .gt. n  .or.  n1 .eq. n2) then
        ier = 1
        return
      endif
c
      lpl = lend(n2)
      if (-list(lpl) .ne. n1) then
        n1 = n2
        n2 = io1
        lpl = lend(n2)
        if (-list(lpl) .ne. n1) then
          ier = 2
          return
        endif
      endif
c
c set n3 to the node opposite n1->n2 (the second neighbor
c   of n1), and test for error 3 (n3 already a boundary
c   node).
c
      lpl = lend(n1)
      lp = lptr(lpl)
      lp = lptr(lp)
      n3 = abs(list(lp))
      lpl = lend(n3)
      if (list(lpl) .le. 0) then
        ier = 3
        return
      endif
c
c delete n2 as a neighbor of n1, making n3 the first
c   neighbor, and test for error 4 (n2 not a neighbor
c   of n1).  note that previously computed pointers may
c   no longer be valid following the call to delnb.
c
      call delnb (n1,n2,n, list,lptr,lend,lnew, lph)
      if (lph .lt. 0) then
        ier = 4
        return
      endif
c
c delete n1 as a neighbor of n2, making n3 the new last
c   neighbor.
c
      call delnb (n2,n1,n, list,lptr,lend,lnew, lph)
c
c make n3 a boundary node with first neighbor n2 and last
c   neighbor n1.
c
      lp = lstptr(lend(n3),n1,list,lptr)
      lend(n3) = lp
      list(lp) = -n1
c
c no errors encountered.
c
      ier = 0
      return
      end
      subroutine delnb (n0,nb,n, list,lptr,lend,lnew, lph)
      integer n0, nb, n, list(*), lptr(*), lend(n), lnew,
     .        lph
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   07/30/98
c
c   this subroutine deletes a neighbor nb from the adjacency
c list of node n0 (but n0 is not deleted from the adjacency
c list of nb) and, if nb is a boundary node, makes n0 a
c boundary node.  for pointer (list index) lph to nb as a
c neighbor of n0, the empty list,lptr location lph is filled
c in with the values at lnew-1, pointer lnew-1 (in lptr and
c possibly in lend) is changed to lph, and lnew is decremen-
c ted.  this requires a search of lend and lptr entailing an
c expected operation count of o(n).
c
c
c on input:
c
c       n0,nb = indexes, in the range 1 to n, of a pair of
c               nodes such that nb is a neighbor of n0.
c               (n0 need not be a neighbor of nb.)
c
c       n = number of nodes in the triangulation.  n .ge. 3.
c
c the above parameters are not altered by this routine.
c
c       list,lptr,lend,lnew = data structure defining the
c                             triangulation.
c
c on output:
c
c       list,lptr,lend,lnew = data structure updated with
c                             the removal of nb from the ad-
c                             jacency list of n0 unless
c                             lph < 0.
c
c       lph = list pointer to the hole (nb as a neighbor of
c             n0) filled in by the values at lnew-1 or error
c             indicator:
c             lph > 0 if no errors were encountered.
c             lph = -1 if n0, nb, or n is outside its valid
c                      range.
c             lph = -2 if nb is not a neighbor of n0.
c
c modules required by delnb:  none
c
c intrinsic function called by delnb:  abs
c
c***********************************************************
c
      integer i, lnw, lp, lpb, lpl, lpp, nn
c
c local parameters:
c
c i =   do-loop index
c lnw = lnew-1 (output value of lnew)
c lp =  list pointer of the last neighbor of nb
c lpb = pointer to nb as a neighbor of n0
c lpl = pointer to the last neighbor of n0
c lpp = pointer to the neighbor of n0 that precedes nb
c nn =  local copy of n
c
      nn = n
c
c test for error 1.
c
      if (n0 .lt. 1  .or.  n0 .gt. nn  .or.  nb .lt. 1  .or.
     .    nb .gt. nn  .or.  nn .lt. 3) then
        lph = -1
        return
      endif
c
c   find pointers to neighbors of n0:
c
c     lpl points to the last neighbor,
c     lpp points to the neighbor np preceding nb, and
c     lpb points to nb.
c
      lpl = lend(n0)
      lpp = lpl
      lpb = lptr(lpp)
    1 if (list(lpb) .eq. nb) go to 2
        lpp = lpb
        lpb = lptr(lpp)
        if (lpb .ne. lpl) go to 1
c
c   test for error 2 (nb not found).
c
      if (abs(list(lpb)) .ne. nb) then
        lph = -2
        return
      endif
c
c   nb is the last neighbor of n0.  make np the new last
c     neighbor and, if nb is a boundary node, then make n0
c     a boundary node.
c
      lend(n0) = lpp
      lp = lend(nb)
      if (list(lp) .lt. 0) list(lpp) = -list(lpp)
      go to 3
c
c   nb is not the last neighbor of n0.  if nb is a boundary
c     node and n0 is not, then make n0 a boundary node with
c     last neighbor np.
c
    2 lp = lend(nb)
      if (list(lp) .lt. 0  .and.  list(lpl) .gt. 0) then
        lend(n0) = lpp
        list(lpp) = -list(lpp)
      endif
c
c   update lptr so that the neighbor following nb now fol-
c     lows np, and fill in the hole at location lpb.
c
    3 lptr(lpp) = lptr(lpb)
      lnw = lnew-1
      list(lpb) = list(lnw)
      lptr(lpb) = lptr(lnw)
      do 4 i = nn,1,-1
        if (lend(i) .eq. lnw) then
          lend(i) = lpb
          go to 5
        endif
    4   continue
c
    5 do 6 i = 1,lnw-1
        if (lptr(i) .eq. lnw) then
          lptr(i) = lpb
        endif
    6   continue
c
c no errors encountered.
c
      lnew = lnw
      lph = lpb
      return
      end
      subroutine delnod (k,ncc, lcc,n,x,y,list,lptr,lend,
     .                   lnew,lwk,iwk, ier)
      integer k, ncc, lcc(*), n, list(*), lptr(*),
     .        lend(*), lnew, lwk, iwk(2,*), ier
      real(8)    x(*), y(*)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   06/28/98
c
c   this subroutine deletes node k (along with all arcs
c incident on node k) from a triangulation of n nodes in the
c plane, and inserts arcs as necessary to produce a triangu-
c lation of the remaining n-1 nodes.  if a delaunay triangu-
c lation is input, a delaunay triangulation will result, and
c thus, delnod reverses the effect of a call to subroutine
c addnod.
c
c   note that a constraint node cannot be deleted by this
c routine.  in order to delete a constraint node, it is
c necessary to call this routine with ncc = 0, decrement the
c appropriate lcc entries (lcc(i) such that lcc(i) > k), and
c then create (or restore) the constraints by a call to sub-
c routine addcst.
c
c
c on input:
c
c       k = index (for x and y) of the node to be deleted.
c           1 .le. k .lt. lcc(1).  (k .le. n if ncc=0).
c
c       ncc = number of constraint curves.  ncc .ge. 0.
c
c the above parameters are not altered by this routine.
c
c       lcc = list of constraint curve starting indexes (or
c             dummy array of length 1 if ncc = 0).  refer to
c             subroutine addcst.
c
c       n = number of nodes in the triangulation on input.
c           n .ge. 4.  note that n will be decremented
c           following the deletion.
c
c       x,y = arrays of length n containing the coordinates
c             of the nodes with non-constraint nodes in the
c             first lcc(1)-1 locations if ncc > 0.
c
c       list,lptr,lend,lnew = data structure defining the
c                             triangulation.  refer to sub-
c                             routine trmesh.
c
c       lwk = number of columns reserved for iwk.  lwk must
c             be at least nnb-3, where nnb is the number of
c             neighbors of node k, including an extra
c             pseudo-node if k is a boundary node.
c
c       iwk = integer work array dimensioned 2 by lwk (or
c             array of length .ge. 2*lwk).
c
c on output:
c
c       lcc = list of constraint curve starting indexes de-
c             cremented by 1 to reflect the deletion of k
c             unless ncc = 0 or 1 .le. ier .le. 4.
c
c       n = new number of nodes (input value minus one) un-
c           less 1 .le. ier .le. 4.
c
c       x,y = updated arrays of length n-1 containing nodal
c             coordinates (with elements k+1,...,n shifted
c             up a position and thus overwriting element k)
c             unless 1 .le. ier .le. 4.  (n here denotes the
c             input value.)
c
c       list,lptr,lend,lnew = updated triangulation data
c                             structure reflecting the dele-
c                             tion unless ier .ne. 0.  note
c                             that the data structure may
c                             have been altered if ier .ge.
c                             3.
c
c       lwk = number of iwk columns required unless ier = 1
c             or ier = 3.
c
c       iwk = indexes of the endpoints of the new arcs added
c             unless lwk = 0 or 1 .le. ier .le. 4.  (arcs
c             are associated with columns, or pairs of
c             adjacent elements if iwk is declared as a
c             singly-subscripted array.)
c
c       ier = error indicator:
c             ier = 0 if no errors were encountered.
c             ier = 1 if k, ncc, n, or an lcc entry is out-
c                     side its valid range or lwk < 0 on
c                     input.
c             ier = 2 if more space is required in iwk.
c                     refer to lwk.
c             ier = 3 if the triangulation data structure is
c                     invalid on input.
c             ier = 4 if k is an interior node with 4 or
c                     more neighbors, and the number of
c                     neighbors could not be reduced to 3
c                     by swaps.  this could be caused by
c                     floating point errors with collinear
c                     nodes or by an invalid data structure.
c             ier = 5 if an error flag was returned by
c                     optim.  an error message is written
c                     to the standard output unit in this
c                     event.
c
c   note that the deletion may result in all remaining nodes
c being collinear.  this situation is not flagged.
c
c modules required by delnod:  delnb, left, lstptr, nbcnt,
c                                optim, swap, swptst
c
c intrinsic function called by delnod:  abs
c
c***********************************************************
c
      integer lstptr, nbcnt
      logical left
      integer i, ierr, iwl, j, lccip1, lnw, lp, lp21, lpf,
     .        lph, lpl, lpl2, lpn, lwkl, n1, n2, nfrst, nit,
     .        nl, nn, nnb, nr
      logical bdry
      real(8)    x1, x2, xl, xr, y1, y2, yl, yr
c
c set n1 to k and nnb to the number of neighbors of n1 (plus
c   one if n1 is a boundary node), and test for errors.  lpf
c   and lpl are list indexes of the first and last neighbors
c   of n1, iwl is the number of iwk columns containing arcs,
c   and bdry is true iff n1 is a boundary node.
c
      n1 = k
      nn = n
      if (ncc .lt. 0  .or.  n1 .lt. 1  .or.  nn .lt. 4  .or.
     .    lwk .lt. 0) go to 21
      lccip1 = nn+1
      do 1 i = ncc,1,-1
        if (lccip1-lcc(i) .lt. 3) go to 21
        lccip1 = lcc(i)
    1   continue
      if (n1 .ge. lccip1) go to 21
      lpl = lend(n1)
      lpf = lptr(lpl)
      nnb = nbcnt(lpl,lptr)
      bdry = list(lpl) .lt. 0
      if (bdry) nnb = nnb + 1
      if (nnb .lt. 3) go to 23
      lwkl = lwk
      lwk = nnb - 3
      if (lwkl .lt. lwk) go to 22
      iwl = 0
      if (nnb .eq. 3) go to 5
c
c initialize for loop on arcs n1-n2 for neighbors n2 of n1,
c   beginning with the second neighbor.  nr and nl are the
c   neighbors preceding and following n2, respectively, and
c   lp indexes nl.  the loop is exited when all possible
c   swaps have been applied to arcs incident on n1.  if n1
c   is interior, the number of neighbors will be reduced
c   to 3.
c
      x1 = x(n1)
      y1 = y(n1)
      nfrst = list(lpf)
      nr = nfrst
      xr = x(nr)
      yr = y(nr)
      lp = lptr(lpf)
      n2 = list(lp)
      x2 = x(n2)
      y2 = y(n2)
      lp = lptr(lp)
c
c top of loop:  set nl to the neighbor following n2.
c
    2 nl = abs(list(lp))
      if (nl .eq. nfrst  .and.  bdry) go to 5
      xl = x(nl)
      yl = y(nl)
c
c   test for a convex quadrilateral.  to avoid an incorrect
c     test caused by collinearity, use the fact that if n1
c     is a boundary node, then n1 left nr->nl and if n2 is
c     a boundary node, then n2 left nl->nr.
c
      lpl2 = lend(n2)
      if ( (bdry  .or.  left(xr,yr,xl,yl,x1,y1))  .and.
     .     (list(lpl2) .lt. 0  .or.
     .      left(xl,yl,xr,yr,x2,y2)) ) go to 3
c
c   nonconvex quadrilateral -- no swap is possible.
c
      nr = n2
      xr = x2
      yr = y2
      go to 4
c
c   the quadrilateral defined by adjacent triangles
c     (n1,n2,nl) and (n2,n1,nr) is convex.  swap in
c     nl-nr and store it in iwk.  indexes larger than n1
c     must be decremented since n1 will be deleted from
c     x and y.
c
    3 call swap (nl,nr,n1,n2, list,lptr,lend, lp21)
      iwl = iwl + 1
      if (nl .le. n1) then
        iwk(1,iwl) = nl
      else
        iwk(1,iwl) = nl - 1
      endif
      if (nr .le. n1) then
        iwk(2,iwl) = nr
      else
        iwk(2,iwl) = nr - 1
      endif
c
c   recompute the list indexes lpl,lp and decrement nnb.
c
      lpl = lend(n1)
      nnb = nnb - 1
      if (nnb .eq. 3) go to 5
      lp = lstptr(lpl,nl,list,lptr)
      if (nr .eq. nfrst) go to 4
c
c   nr is not the first neighbor of n1.
c     back up and test n1-nr for a swap again:  set n2 to
c     nr and nr to the previous neighbor of n1 -- the
c     neighbor of nr which follows n1.  lp21 points to nl
c     as a neighbor of nr.
c
      n2 = nr
      x2 = xr
      y2 = yr
      lp21 = lptr(lp21)
      lp21 = lptr(lp21)
      nr = abs(list(lp21))
      xr = x(nr)
      yr = y(nr)
      go to 2
c
c   bottom of loop -- test for invalid termination.
c
    4 if (n2 .eq. nfrst) go to 24
      n2 = nl
      x2 = xl
      y2 = yl
      lp = lptr(lp)
      go to 2
c
c delete n1 from the adjacency list of n2 for all neighbors
c   n2 of n1.  lpl points to the last neighbor of n1.
c   lnew is stored in local variable lnw.
c
    5 lp = lpl
      lnw = lnew
c
c loop on neighbors n2 of n1, beginning with the first.
c
    6 lp = lptr(lp)
        n2 = abs(list(lp))
        call delnb (n2,n1,n, list,lptr,lend,lnw, lph)
        if (lph .lt. 0) go to 23
c
c   lp and lpl may require alteration.
c
        if (lpl .eq. lnw) lpl = lph
        if (lp .eq. lnw) lp = lph
        if (lp .ne. lpl) go to 6
c
c delete n1 from x, y, and lend, and remove its adjacency
c   list from list and lptr.  list entries (nodal indexes)
c   which are larger than n1 must be decremented.
c
      nn = nn - 1
      if (n1 .gt. nn) go to 9
      do 7 i = n1,nn
        x(i) = x(i+1)
        y(i) = y(i+1)
        lend(i) = lend(i+1)
    7   continue
c
      do 8 i = 1,lnw-1
        if (list(i) .gt. n1) list(i) = list(i) - 1
        if (list(i) .lt. -n1) list(i) = list(i) + 1
    8   continue
c
c   for lpn = first to last neighbors of n1, delete the
c     preceding neighbor (indexed by lp).
c
c   each empty list,lptr location lp is filled in with the
c     values at lnw-1, and lnw is decremented.  all pointers
c     (including those in lptr and lend) with value lnw-1
c     must be changed to lp.
c
c  lpl points to the last neighbor of n1.
c
    9 if (bdry) nnb = nnb - 1
      lpn = lpl
      do 13 j = 1,nnb
        lnw = lnw - 1
        lp = lpn
        lpn = lptr(lp)
        list(lp) = list(lnw)
        lptr(lp) = lptr(lnw)
        if (lptr(lpn) .eq. lnw) lptr(lpn) = lp
        if (lpn .eq. lnw) lpn = lp
        do 10 i = nn,1,-1
          if (lend(i) .eq. lnw) then
            lend(i) = lp
            go to 11
          endif
   10     continue
c
   11   do 12 i = lnw-1,1,-1
          if (lptr(i) .eq. lnw) lptr(i) = lp
   12     continue
   13   continue
c
c decrement lcc entries.
c
      do 14 i = 1,ncc
        lcc(i) = lcc(i) - 1
   14   continue
c
c update n and lnew, and optimize the patch of triangles
c   containing k (on input) by applying swaps to the arcs
c   in iwk.
c
      n = nn
      lnew = lnw
      if (iwl .gt. 0) then
        nit = 4*iwl
        call optim (x,y,iwl, list,lptr,lend,nit,iwk, ierr)
        if (ierr .ne. 0) go to 25
      endif
c
c successful termination.
c
      ier = 0
      return
c
c invalid input parameter.
c
   21 ier = 1
      return
c
c insufficient space reserved for iwk.
c
   22 ier = 2
      return
c
c invalid triangulation data structure.  nnb < 3 on input or
c   n2 is a neighbor of n1 but n1 is not a neighbor of n2.
c
   23 ier = 3
      return
c
c k is an interior node with 4 or more neighbors, but the
c   number of neighbors could not be reduced.
c
   24 ier = 4
      return
c
c error flag returned by optim.
c
   25 ier = 5
      write (*,100) nit, ierr
      return
  100 format (//5x,'*** error in optim:  nit = ',i4,
     .        ', ier = ',i1,' ***'/)
      end
      subroutine edge (in1,in2,x,y, lwk,iwk,list,lptr,
     .                 lend, ier)
      integer in1, in2, lwk, iwk(2,*), list(*), lptr(*),
     .        lend(*), ier
      real(8)    x(*), y(*)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   06/23/98
c
c   given a triangulation of n nodes and a pair of nodal
c indexes in1 and in2, this routine swaps arcs as necessary
c to force in1 and in2 to be adjacent.  only arcs which
c intersect in1-in2 are swapped out.  if a delaunay triangu-
c lation is input, the resulting triangulation is as close
c as possible to a delaunay triangulation in the sense that
c all arcs other than in1-in2 are locally optimal.
c
c   a sequence of calls to edge may be used to force the
c presence of a set of edges defining the boundary of a non-
c convex and/or multiply connected region (refer to subrou-
c tine addcst), or to introduce barriers into the triangula-
c tion.  note that subroutine getnp will not necessarily
c return closest nodes if the triangulation has been con-
c strained by a call to edge.  however, this is appropriate
c in some applications, such as triangle-based interpolation
c on a nonconvex domain.
c
c
c on input:
c
c       in1,in2 = indexes (of x and y) in the range 1 to n
c                 defining a pair of nodes to be connected
c                 by an arc.
c
c       x,y = arrays of length n containing the cartesian
c             coordinates of the nodes.
c
c the above parameters are not altered by this routine.
c
c       lwk = number of columns reserved for iwk.  this must
c             be at least ni -- the number of arcs which
c             intersect in1-in2.  (ni is bounded by n-3.)
c
c       iwk = integer work array of length at least 2*lwk.
c
c       list,lptr,lend = data structure defining the trian-
c                        gulation.  refer to subroutine
c                        trmesh.
c
c on output:
c
c       lwk = number of arcs which intersect in1-in2 (but
c             not more than the input value of lwk) unless
c             ier = 1 or ier = 3.  lwk = 0 if and only if
c             in1 and in2 were adjacent (or lwk=0) on input.
c
c       iwk = array containing the indexes of the endpoints
c             of the new arcs other than in1-in2 unless ier
c             .gt. 0 or lwk = 0.  new arcs to the left of
c             in2-in1 are stored in the first k-1 columns
c             (left portion of iwk), column k contains
c             zeros, and new arcs to the right of in2-in1
c             occupy columns k+1,...,lwk.  (k can be deter-
c             mined by searching iwk for the zeros.)
c
c       list,lptr,lend = data structure updated if necessary
c                        to reflect the presence of an arc
c                        connecting in1 and in2 unless ier
c                        .ne. 0.  the data structure has
c                        been altered if ier = 4.
c
c       ier = error indicator:
c             ier = 0 if no errors were encountered.
c             ier = 1 if in1 .lt. 1, in2 .lt. 1, in1 = in2,
c                     or lwk .lt. 0 on input.
c             ier = 2 if more space is required in iwk.
c             ier = 3 if in1 and in2 could not be connected
c                     due to either an invalid data struc-
c                     ture or collinear nodes (and floating
c                     point error).
c             ier = 4 if an error flag was returned by
c                     optim.
c
c   an error message is written to the standard output unit
c in the case of ier = 3 or ier = 4.
c
c modules required by edge:  left, lstptr, optim, swap,
c                              swptst
c
c intrinsic function called by edge:  abs
c
c***********************************************************
c
      logical left
      integer i, ierr, iwc, iwcp1, iwend, iwf, iwl, lft, lp,
     .        lpl, lp21, next, nit, nl, nr, n0, n1, n2,
     .        n1frst, n1lst
      real(8)    dx, dy, x0, y0, x1, y1, x2, y2
c
c local parameters:
c
c dx,dy =   components of arc n1-n2
c i =       do-loop index and column index for iwk
c ierr =    error flag returned by subroutine optim
c iwc =     iwk index between iwf and iwl -- nl->nr is
c             stored in iwk(1,iwc)->iwk(2,iwc)
c iwcp1 =   iwc + 1
c iwend =   input or output value of lwk
c iwf =     iwk (column) index of the first (leftmost) arc
c             which intersects in1->in2
c iwl =     iwk (column) index of the last (rightmost) are
c             which intersects in1->in2
c lft =     flag used to determine if a swap results in the
c             new arc intersecting in1-in2 -- lft = 0 iff
c             n0 = in1, lft = -1 implies n0 left in1->in2,
c             and lft = 1 implies n0 left in2->in1
c lp21 =    unused parameter returned by swap
c lp =      list pointer (index) for list and lptr
c lpl =     pointer to the last neighbor of in1 or nl
c n0 =      neighbor of n1 or node opposite nr->nl
c n1,n2 =   local copies of in1 and in2
c n1frst =  first neighbor of in1
c n1lst =   (signed) last neighbor of in1
c next =    node opposite nl->nr
c nit =     flag or number of iterations employed by optim
c nl,nr =   endpoints of an arc which intersects in1-in2
c             with nl left in1->in2
c x0,y0 =   coordinates of n0
c x1,y1 =   coordinates of in1
c x2,y2 =   coordinates of in2
c
c
c store in1, in2, and lwk in local variables and test for
c   errors.
c
      n1 = in1
      n2 = in2
      iwend = lwk
      if (n1 .lt. 1  .or.  n2 .lt. 1  .or.  n1 .eq. n2  .or.
     .    iwend .lt. 0) go to 31
c
c test for n2 as a neighbor of n1.  lpl points to the last
c   neighbor of n1.
c
      lpl = lend(n1)
      n0 = abs(list(lpl))
      lp = lpl
    1 if (n0 .eq. n2) go to 30
        lp = lptr(lp)
        n0 = list(lp)
        if (lp .ne. lpl) go to 1
c
c initialize parameters.
c
      iwl = 0
      nit = 0
c
c store the coordinates of n1 and n2.
c
    2 x1 = x(n1)
      y1 = y(n1)
      x2 = x(n2)
      y2 = y(n2)
c
c set nr and nl to adjacent neighbors of n1 such that
c   nr left n2->n1 and nl left n1->n2,
c   (nr forward n1->n2 or nl forward n1->n2), and
c   (nr forward n2->n1 or nl forward n2->n1).
c
c   initialization:  set n1frst and n1lst to the first and
c     (signed) last neighbors of n1, respectively, and
c     initialize nl to n1frst.
c
      lpl = lend(n1)
      n1lst = list(lpl)
      lp = lptr(lpl)
      n1frst = list(lp)
      nl = n1frst
      if (n1lst .lt. 0) go to 4
c
c   n1 is an interior node.  set nl to the first candidate
c     for nr (nl left n2->n1).
c
    3 if ( left(x2,y2,x1,y1,x(nl),y(nl)) ) go to 4
        lp = lptr(lp)
        nl = list(lp)
        if (nl .ne. n1frst) go to 3
c
c   all neighbors of n1 are strictly left of n1->n2.
c
      go to 5
c
c   nl = list(lp) left n2->n1.  set nr to nl and nl to the
c     following neighbor of n1.
c
    4 nr = nl
        lp = lptr(lp)
        nl = abs(list(lp))
        if ( left(x1,y1,x2,y2,x(nl),y(nl)) ) then
c
c   nl left n1->n2 and nr left n2->n1.  the forward tests
c     are employed to avoid an error associated with
c     collinear nodes.
c
          dx = x2-x1
          dy = y2-y1
          if ((dx*(x(nl)-x1)+dy*(y(nl)-y1) .ge. 0.  .or.
     .         dx*(x(nr)-x1)+dy*(y(nr)-y1) .ge. 0.)  .and.
     .        (dx*(x(nl)-x2)+dy*(y(nl)-y2) .le. 0.  .or.
     .         dx*(x(nr)-x2)+dy*(y(nr)-y2) .le. 0.)) go to 6
c
c   nl-nr does not intersect n1-n2.  however, there is
c     another candidate for the first arc if nl lies on
c     the line n1-n2.
c
          if ( .not. left(x2,y2,x1,y1,x(nl),y(nl)) ) go to 5
        endif
c
c   bottom of loop.
c
        if (nl .ne. n1frst) go to 4
c
c either the triangulation is invalid or n1-n2 lies on the
c   convex hull boundary and an edge nr->nl (opposite n1 and
c   intersecting n1-n2) was not found due to floating point
c   error.  try interchanging n1 and n2 -- nit > 0 iff this
c   has already been done.
c
    5 if (nit .gt. 0) go to 33
      nit = 1
      n1 = n2
      n2 = in1
      go to 2
c
c store the ordered sequence of intersecting edges nl->nr in
c   iwk(1,iwl)->iwk(2,iwl).
c
    6 iwl = iwl + 1
      if (iwl .gt. iwend) go to 32
      iwk(1,iwl) = nl
      iwk(2,iwl) = nr
c
c   set next to the neighbor of nl which follows nr.
c
      lpl = lend(nl)
      lp = lptr(lpl)
c
c   find nr as a neighbor of nl.  the search begins with
c     the first neighbor.
c
    7 if (list(lp) .eq. nr) go to 8
        lp = lptr(lp)
        if (lp .ne. lpl) go to 7
c
c   nr must be the last neighbor, and nl->nr cannot be a
c     boundary edge.
c
      if (list(lp) .ne. nr) go to 33
c
c   set next to the neighbor following nr, and test for
c     termination of the store loop.
c
    8 lp = lptr(lp)
      next = abs(list(lp))
      if (next .eq. n2) go to 9
c
c   set nl or nr to next.
c
      if ( left(x1,y1,x2,y2,x(next),y(next)) ) then
        nl = next
      else
        nr = next
      endif
      go to 6
c
c iwl is the number of arcs which intersect n1-n2.
c   store lwk.
c
    9 lwk = iwl
      iwend = iwl
c
c initialize for edge swapping loop -- all possible swaps
c   are applied (even if the new arc again intersects
c   n1-n2), arcs to the left of n1->n2 are stored in the
c   left portion of iwk, and arcs to the right are stored in
c   the right portion.  iwf and iwl index the first and last
c   intersecting arcs.
c
      iwf = 1
c
c top of loop -- set n0 to n1 and nl->nr to the first edge.
c   iwc points to the arc currently being processed.  lft
c   .le. 0 iff n0 left n1->n2.
c
   10 lft = 0
      n0 = n1
      x0 = x1
      y0 = y1
      nl = iwk(1,iwf)
      nr = iwk(2,iwf)
      iwc = iwf
c
c   set next to the node opposite nl->nr unless iwc is the
c     last arc.
c
   11 if (iwc .eq. iwl) go to 21
      iwcp1 = iwc + 1
      next = iwk(1,iwcp1)
      if (next .ne. nl) go to 16
      next = iwk(2,iwcp1)
c
c   next right n1->n2 and iwc .lt. iwl.  test for a possible
c     swap.
c
      if ( .not. left(x0,y0,x(nr),y(nr),x(next),y(next)) )
     .   go to 14
      if (lft .ge. 0) go to 12
      if ( .not. left(x(nl),y(nl),x0,y0,x(next),y(next)) )
     .   go to 14
c
c   replace nl->nr with n0->next.
c
      call swap (next,n0,nl,nr, list,lptr,lend, lp21)
      iwk(1,iwc) = n0
      iwk(2,iwc) = next
      go to 15
c
c   swap nl-nr for n0-next, shift columns iwc+1,...,iwl to
c     the left, and store n0-next in the right portion of
c     iwk.
c
   12 call swap (next,n0,nl,nr, list,lptr,lend, lp21)
      do 13 i = iwcp1,iwl
        iwk(1,i-1) = iwk(1,i)
        iwk(2,i-1) = iwk(2,i)
   13   continue
      iwk(1,iwl) = n0
      iwk(2,iwl) = next
      iwl = iwl - 1
      nr = next
      go to 11
c
c   a swap is not possible.  set n0 to nr.
c
   14 n0 = nr
      x0 = x(n0)
      y0 = y(n0)
      lft = 1
c
c   advance to the next arc.
c
   15 nr = next
      iwc = iwc + 1
      go to 11
c
c   next left n1->n2, next .ne. n2, and iwc .lt. iwl.
c     test for a possible swap.
c
   16 if ( .not. left(x(nl),y(nl),x0,y0,x(next),y(next)) )
     .   go to 19
      if (lft .le. 0) go to 17
      if ( .not. left(x0,y0,x(nr),y(nr),x(next),y(next)) )
     .   go to 19
c
c   replace nl->nr with next->n0.
c
      call swap (next,n0,nl,nr, list,lptr,lend, lp21)
      iwk(1,iwc) = next
      iwk(2,iwc) = n0
      go to 20
c
c   swap nl-nr for n0-next, shift columns iwf,...,iwc-1 to
c     the right, and store n0-next in the left portion of
c     iwk.
c
   17 call swap (next,n0,nl,nr, list,lptr,lend, lp21)
      do 18 i = iwc-1,iwf,-1
        iwk(1,i+1) = iwk(1,i)
        iwk(2,i+1) = iwk(2,i)
   18   continue
      iwk(1,iwf) = n0
      iwk(2,iwf) = next
      iwf = iwf + 1
      go to 20
c
c   a swap is not possible.  set n0 to nl.
c
   19 n0 = nl
      x0 = x(n0)
      y0 = y(n0)
      lft = -1
c
c   advance to the next arc.
c
   20 nl = next
      iwc = iwc + 1
      go to 11
c
c   n2 is opposite nl->nr (iwc = iwl).
c
   21 if (n0 .eq. n1) go to 24
      if (lft .lt. 0) go to 22
c
c   n0 right n1->n2.  test for a possible swap.
c
      if ( .not. left(x0,y0,x(nr),y(nr),x2,y2) ) go to 10
c
c   swap nl-nr for n0-n2 and store n0-n2 in the right
c     portion of iwk.
c
      call swap (n2,n0,nl,nr, list,lptr,lend, lp21)
      iwk(1,iwl) = n0
      iwk(2,iwl) = n2
      iwl = iwl - 1
      go to 10
c
c   n0 left n1->n2.  test for a possible swap.
c
   22 if ( .not. left(x(nl),y(nl),x0,y0,x2,y2) ) go to 10
c
c   swap nl-nr for n0-n2, shift columns iwf,...,iwl-1 to the
c     right, and store n0-n2 in the left portion of iwk.
c
      call swap (n2,n0,nl,nr, list,lptr,lend, lp21)
      i = iwl
   23 iwk(1,i) = iwk(1,i-1)
      iwk(2,i) = iwk(2,i-1)
      i = i - 1
      if (i .gt. iwf) go to 23
      iwk(1,iwf) = n0
      iwk(2,iwf) = n2
      iwf = iwf + 1
      go to 10
c
c iwf = iwc = iwl.  swap out the last arc for n1-n2 and
c   store zeros in iwk.
c
   24 call swap (n2,n1,nl,nr, list,lptr,lend, lp21)
      iwk(1,iwc) = 0
      iwk(2,iwc) = 0
c
c optimization procedure --
c
      if (iwc .gt. 1) then
c
c   optimize the set of new arcs to the left of in1->in2.
c
        nit = 3*(iwc-1)
        call optim (x,y,iwc-1, list,lptr,lend,nit,iwk, ierr)
        if (ierr .ne. 0) go to 34
      endif
      if (iwc .lt. iwend) then
c
c   optimize the set of new arcs to the right of in1->in2.
c
        nit = 3*(iwend-iwc)
        call optim (x,y,iwend-iwc, list,lptr,lend,nit,
     .              iwk(1,iwc+1), ierr)
        if (ierr .ne. 0) go to 34
      endif
c
c successful termination.
c
      ier = 0
      return
c
c in1 and in2 were adjacent on input.
c
   30 ier = 0
      return
c
c invalid input parameter.
c
   31 ier = 1
      return
c
c insufficient space reserved for iwk.
c
   32 ier = 2
      return
c
c invalid triangulation data structure or collinear nodes
c   on convex hull boundary.
c
   33 ier = 3
      write (*,130) in1, in2
  130 format (//5x,'*** error in edge:  invalid triangula',
     .        'tion or null triangles on boundary'/
     .        9x,'in1 =',i4,', in2=',i4/)
      return
c
c error flag returned by optim.
c
   34 ier = 4
      write (*,140) nit, ierr
  140 format (//5x,'*** error in optim:  nit = ',i4,
     .        ', ier = ',i1,' ***'/)
      return
      end
      subroutine getnp (ncc,lcc,n,x,y,list,lptr,lend,
     .                  l, npts,ds, ier)
      integer ncc, lcc(*), n, list(*), lptr(*), lend(n),
     .        l, npts(l), ier
      real(8)    x(n), y(n), ds(l)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   11/12/94
c
c   given a triangulation of n nodes and an array npts con-
c taining the indexes of l-1 nodes ordered by distance from
c npts(1), this subroutine sets npts(l) to the index of the
c next node in the sequence -- the node, other than npts(1),
c ...,npts(l-1), which is closest to npts(1).  thus, the
c ordered sequence of k closest nodes to n1 (including n1)
c may be determined by k-1 calls to getnp with npts(1) = n1
c and l = 2,3,...,k for k .ge. 2.  note that npts must in-
c clude constraint nodes as well as non-constraint nodes.
c thus, a sequence of k1 closest non-constraint nodes to n1
c must be obtained as a subset of the closest k2 nodes to n1
c for some k2 .ge. k1.
c
c   the terms closest and distance have special definitions
c when constraint nodes are present in the triangulation.
c nodes n1 and n2 are said to be visible from each other if
c and only if the line segment n1-n2 intersects no con-
c straint arc (except possibly itself) and is not an interi-
c or constraint arc (arc whose interior lies in a constraint
c region).  a path from n1 to n2 is an ordered sequence of
c nodes, with n1 first and n2 last, such that adjacent path
c elements are visible from each other.  the path length is
c the sum of the euclidean distances between adjacent path
c nodes.  finally, the distance from n1 to n2 is defined to
c be the length of the shortest path from n1 to n2.
c
c   the algorithm uses the property of a delaunay triangula-
c tion that the k-th closest node to n1 is a neighbor of one
c of the k-1 closest nodes to n1.  with the definition of
c distance used here, this property holds when constraints
c are present as long as non-constraint arcs are locally
c optimal.
c
c
c on input:
c
c       ncc = number of constraints.  ncc .ge. 0.
c
c       lcc = list of constraint curve starting indexes (or
c             dummy array of length 1 if ncc = 0).  refer to
c             subroutine addcst.
c
c       n = number of nodes in the triangulation.  n .ge. 3.
c
c       x,y = arrays of length n containing the coordinates
c             of the nodes with non-constraint nodes in the
c             first lcc(1)-1 locations if ncc > 0.
c
c       list,lptr,lend = triangulation data structure.  re-
c                        fer to subroutine trmesh.
c
c       l = number of nodes in the sequence on output.  2
c           .le. l .le. n.
c
c       npts = array of length .ge. l containing the indexes
c              of the l-1 closest nodes to npts(1) in the
c              first l-1 locations.
c
c       ds = array of length .ge. l containing the distance
c            (defined above) between npts(1) and npts(i) in
c            the i-th position for i = 1,...,l-1.  thus,
c            ds(1) = 0.
c
c input parameters other than npts(l) and ds(l) are not
c   altered by this routine.
c
c on output:
c
c       npts = array updated with the index of the l-th
c              closest node to npts(1) in position l unless
c              ier .ne. 0.
c
c       ds = array updated with the distance between npts(1)
c            and npts(l) in position l unless ier .ne. 0.
c
c       ier = error indicator:
c             ier =  0 if no errors were encountered.
c             ier = -1 if ncc, n, l, or an lcc entry is
c                      outside its valid range on input.
c             ier =  k if npts(k) is not a valid index in
c                      the range 1 to n.
c
c module required by getnp:  intsec
c
c intrinsic functions called by getnp:  abs, min, sqrt
c
c***********************************************************
c
      logical intsec
      integer i, ifrst, ilast, j, k, km1, lcc1, lm1, lp,
     .        lpcl, lpk, lpkl, n1, nc, nf1, nf2, nj, nk,
     .        nkbak, nkfor, nl, nn
      logical isw, vis, ncf, njf, skip, sksav, lft1, lft2,
     .        lft12
      real(8)    dc, dl, x1, xc, xj, xk, y1, yc, yj, yk
c
c store parameters in local variables and test for errors.
c   lcc1 indexes the first constraint node.
c
      ier = -1
      nn = n
      lcc1 = nn+1
      lm1 = l-1
      if (ncc .lt. 0  .or.  lm1 .lt. 1  .or.  lm1 .ge. nn)
     .   return
      if (ncc .eq. 0) then
        if (nn .lt. 3) return
      else
        do 1 i = ncc,1,-1
          if (lcc1 - lcc(i) .lt. 3) return
          lcc1 = lcc(i)
    1     continue
        if (lcc1 .lt. 1) return
      endif
c
c test for an invalid index in npts.
c
      do 2 k = 1,lm1
        nk = npts(k)
        if (nk .lt. 1  .or.  nk .gt. nn) then
          ier = k
          return
        endif
    2   continue
c
c store n1 = npts(1) and mark the elements of npts.
c
      n1 = npts(1)
      x1 = x(n1)
      y1 = y(n1)
      do 3 k = 1,lm1
        nk = npts(k)
        lend(nk) = -lend(nk)
    3   continue
c
c candidates nc for nl = npts(l) are the unmarked visible
c   neighbors of nodes nk in npts.  isw is an initialization
c   switch set to .true. when nl and its distance dl from n1
c   have been initialized with the first candidate encount-
c   ered.
c
      isw = .false.
      dl = 0.
c
c loop on marked nodes nk = npts(k).  lpkl indexes the last
c   neighbor of nk in list.
c
      do 16 k = 1,lm1
        km1 = k - 1
        nk = npts(k)
        xk = x(nk)
        yk = y(nk)
        lpkl = -lend(nk)
        nkfor = 0
        nkbak = 0
        vis = .true.
        if (nk .ge. lcc1) then
c
c   nk is a constraint node.  set nkfor and nkbak to the
c     constraint nodes which follow and precede nk.  ifrst
c     and ilast are set to the first and last nodes in the
c     constraint containing nk.
c
          ifrst = nn + 1
          do 4 i = ncc,1,-1
            ilast = ifrst - 1
            ifrst = lcc(i)
            if (nk .ge. ifrst) go to 5
    4       continue
c
    5     if (nk .lt. ilast) then
            nkfor = nk + 1
          else
            nkfor = ifrst
          endif
          if (nk .gt. ifrst) then
            nkbak = nk - 1
          else
            nkbak = ilast
          endif
c
c   initialize vis to true iff nkfor precedes nkbak in the
c     adjacency list for nk -- the first neighbor is visi-
c     ble and is not nkbak.
c
          lpk = lpkl
    6     lpk = lptr(lpk)
            nc = abs(list(lpk))
            if (nc .ne. nkfor  .and.  nc .ne. nkbak) go to 6
          vis = nc .eq. nkfor
        endif
c
c loop on neighbors nc of nk, bypassing marked and nonvis-
c   ible neighbors.
c
        lpk = lpkl
    7   lpk = lptr(lpk)
          nc = abs(list(lpk))
          if (nc .eq. nkbak) vis = .true.
c
c   vis = .false. iff nk-nc is an interior constraint arc
c     (nk is a constraint node and nc lies strictly between
c     nkfor and nkbak).
c
          if (.not. vis) go to 15
          if (nc .eq. nkfor) vis = .false.
          if (lend(nc) .lt. 0) go to 15
c
c initialize distance dc between n1 and nc to euclidean
c   distance.
c
          xc = x(nc)
          yc = y(nc)
          dc = sqrt((xc-x1)*(xc-x1) + (yc-y1)*(yc-y1))
          if (isw  .and.  dc .ge. dl) go to 15
          if (k .eq. 1) go to 14
c
c k .ge. 2.  store the pointer lpcl to the last neighbor
c   of nc.
c
          lpcl = lend(nc)
c
c set dc to the length of the shortest path from n1 to nc
c   which has not previously been encountered and which is
c   a viable candidate for the shortest path from n1 to nl.
c   this is euclidean distance iff nc is visible from n1.
c   since the shortest path from n1 to nl contains only ele-
c   ments of npts which are constraint nodes (in addition to
c   n1 and nl), only these need be considered for the path
c   from n1 to nc.  thus, for distance function d(a,b) and
c   j = 1,...,k, dc = min(d(n1,nj) + d(nj,nc)) over con-
c   straint nodes nj = npts(j) which are visible from nc.
c
          do 13 j = 1,km1
            nj = npts(j)
            if (j .gt. 1  .and.  nj .lt. lcc1) go to 13
c
c if nc is a visible neighbor of nj, a path from n1 to nc
c   containing nj has already been considered.  thus, nj may
c   be bypassed if it is adjacent to nc.
c
            lp = lpcl
    8       lp = lptr(lp)
              if ( nj .eq. abs(list(lp)) ) go to 12
              if (lp .ne. lpcl) go to 8
c
c nj is a constraint node (unless j=1) not adjacent to nc,
c   and is visible from nc iff nj-nc is not intersected by
c   a constraint arc.  loop on constraints i in reverse
c   order --
c
            xj = x(nj)
            yj = y(nj)
            ifrst = nn+1
            do 11 i = ncc,1,-1
              ilast = ifrst - 1
              ifrst = lcc(i)
              nf1 = ilast
              ncf = nf1 .eq. nc
              njf = nf1 .eq. nj
              skip = ncf  .or.  njf
c
c loop on boundary constraint arcs nf1-nf2 which contain
c   neither nc nor nj.  ncf and njf are true iff nc (or nj)
c   has been encountered in the constraint, and skip =
c   .true. iff nf1 = nc or nf1 = nj.
c
              do 10 nf2 = ifrst,ilast
                if (nf2 .eq. nc) ncf = .true.
                if (nf2 .eq. nj) njf = .true.
                sksav = skip
                skip = nf2 .eq. nc  .or.  nf2 .eq. nj
c
c   the last constraint arc in the constraint need not be
c     tested if none of the arcs have been skipped.
c
                if ( sksav  .or.  skip  .or.
     .               (nf2 .eq. ilast  .and.
     .               .not. ncf  .and.  .not. njf) ) go to 9
                if ( intsec(x(nf1),y(nf1),x(nf2),y(nf2),
     .                      xc,yc,xj,yj) ) go to 12
    9           nf1 = nf2
   10           continue
              if (.not. ncf  .or.  .not. njf) go to 11
c
c nc and nj are constraint nodes in the same constraint.
c   nc-nj is intersected by an interior constraint arc iff
c   1)  nc left nf2->nf1 and (nj left nf1->nc and nj left
c         nc->nf2) or
c   2)  nc .not. left nf2->nf1 and (nj left nf1->nc or
c         nj left nc->nf2),
c   where nf1, nc, nf2 are consecutive constraint nodes.
c
              if (nc .ne. ifrst) then
                nf1 = nc - 1
              else
                nf1 = ilast
              endif
              if (nc .ne. ilast) then
                nf2 = nc + 1
              else
                nf2 = ifrst
              endif
              lft1 = (xc-x(nf1))*(yj-y(nf1)) .ge.
     .               (xj-x(nf1))*(yc-y(nf1))
              lft2 = (x(nf2)-xc)*(yj-yc) .ge.
     .               (xj-xc)*(y(nf2)-yc)
              lft12 = (x(nf1)-x(nf2))*(yc-y(nf2)) .ge.
     .                (xc-x(nf2))*(y(nf1)-y(nf2))
              if ( (lft1  .and.  lft2)  .or.  (.not. lft12
     .             .and.  (lft1  .or.  lft2)) ) go to 12
   11         continue
c
c nj is visible from nc.  exit the loop with dc = euclidean
c   distance if j = 1.
c
            if (j .eq. 1) go to 14
            dc = min(dc,ds(j) + sqrt((xc-xj)*(xc-xj) +
     .                  (yc-yj)*(yc-yj)))
            go to 13
c
c nj is not visible from nc or is adjacent to nc.  initial-
c   ize dc with d(n1,nk) + d(nk,nc) if j = 1.
c
   12       if (j .eq. 1) dc = ds(k) + sqrt((xc-xk)*(xc-xk)
     .                         + (yc-yk)*(yc-yk))
   13       continue
c
c compare dc with dl.
c
          if (isw  .and.  dc .ge. dl) go to 15
c
c the first (or a closer) candidate for nl has been
c   encountered.
c
   14     nl = nc
          dl = dc
          isw = .true.
   15     if (lpk .ne. lpkl) go to 7
   16   continue
c
c unmark the elements of npts and store nl and dl.
c
      do 17 k = 1,lm1
        nk = npts(k)
        lend(nk) = -lend(nk)
   17   continue
      npts(l) = nl
      ds(l) = dl
      ier = 0
      return
      end
      integer function indxcc (ncc,lcc,n,list,lend)
      integer ncc, lcc(*), n, list(*), lend(n)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   08/25/91
c
c   given a constrained delaunay triangulation, this func-
c tion returns the index, if any, of an exterior constraint
c curve (an unbounded constraint region).  an exterior con-
c straint curve is assumed to be present if and only if the
c clockwise-ordered sequence of boundary nodes is a subse-
c quence of a constraint node sequence.  the triangulation
c adjacencies corresponding to constraint edges may or may
c not have been forced by a call to addcst, and the con-
c straint region may or may not be valid (contain no nodes).
c
c
c on input:
c
c       ncc = number of constraints.  ncc .ge. 0.
c
c       lcc = list of constraint curve starting indexes (or
c             dummy array of length 1 if ncc = 0).  refer to
c             subroutine addcst.
c
c       n = number of nodes in the triangulation.  n .ge. 3.
c
c       list,lend = data structure defining the triangula-
c                   tion.  refer to subroutine trmesh.
c
c   input parameters are not altered by this function.  note
c that the parameters are not tested for validity.
c
c on output:
c
c       indxcc = index of the exterior constraint curve, if
c                present, or 0 otherwise.
c
c modules required by indxcc:  none
c
c***********************************************************
c
      integer i, ifrst, ilast, lp, n0, nst, nxt
      indxcc = 0
      if (ncc .lt. 1) return
c
c set n0 to the boundary node with smallest index.
c
      n0 = 0
    1 n0 = n0 + 1
        lp = lend(n0)
        if (list(lp) .gt. 0) go to 1
c
c search in reverse order for the constraint i, if any, that
c   contains n0.  ifrst and ilast index the first and last
c   nodes in constraint i.
c
      i = ncc
      ilast = n
    2 ifrst = lcc(i)
        if (n0 .ge. ifrst) go to 3
        if (i .eq. 1) return
        i = i - 1
        ilast = ifrst - 1
        go to 2
c
c n0 is in constraint i which indexes an exterior constraint
c   curve iff the clockwise-ordered sequence of boundary
c   node indexes beginning with n0 is increasing and bounded
c   above by ilast.
c
    3 nst = n0
c
    4 nxt = -list(lp)
        if (nxt .eq. nst) go to 5
        if (nxt .le. n0  .or.  nxt .gt. ilast) return
        n0 = nxt
        lp = lend(n0)
        go to 4
c
c constraint i contains the boundary node sequence as a
c   subset.
c
    5 indxcc = i
      return
      end
      subroutine insert (k,lp, list,lptr,lnew )
      integer k, lp, list(*), lptr(*), lnew
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   09/01/88
c
c   this subroutine inserts k as a neighbor of n1 following
c n2, where lp is the list pointer of n2 as a neighbor of
c n1.  note that, if n2 is the last neighbor of n1, k will
c become the first neighbor (even if n1 is a boundary node).
c
c
c on input:
c
c       k = index of the node to be inserted.
c
c       lp = list pointer of n2 as a neighbor of n1.
c
c the above parameters are not altered by this routine.
c
c       list,lptr,lnew = data structure defining the trian-
c                        gulation.  refer to subroutine
c                        trmesh.
c
c on output:
c
c       list,lptr,lnew = data structure updated with the
c                        addition of node k.
c
c modules required by insert:  none
c
c***********************************************************
c
      integer lsav
c
      lsav = lptr(lp)
      lptr(lp) = lnew
      list(lnew) = k
      lptr(lnew) = lsav
      lnew = lnew + 1
      return
      end
      subroutine intadd (kk,i1,i2,i3, list,lptr,lend,lnew )
      integer kk, i1, i2, i3, list(*), lptr(*), lend(*),
     .        lnew
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   02/22/91
c
c   this subroutine adds an interior node to a triangulation
c of a set of points in the plane.  the data structure is
c updated with the insertion of node kk into the triangle
c whose vertices are i1, i2, and i3.  no optimization of the
c triangulation is performed.
c
c
c on input:
c
c       kk = index of the node to be inserted.  kk .ge. 1
c            and kk must not be equal to i1, i2, or i3.
c
c       i1,i2,i3 = indexes of the counterclockwise-ordered
c                  sequence of vertices of a triangle which
c                  contains node kk.
c
c the above parameters are not altered by this routine.
c
c       list,lptr,lend,lnew = data structure defining the
c                             triangulation.  refer to sub-
c                             routine trmesh.  triangle
c                             (i1,i2,i3) must be included
c                             in the triangulation.
c
c on output:
c
c       list,lptr,lend,lnew = data structure updated with
c                             the addition of node kk.  kk
c                             will be connected to nodes i1,
c                             i2, and i3.
c
c modules required by intadd:  insert, lstptr
c
c***********************************************************
c
      integer lstptr
      integer k, lp, n1, n2, n3
      k = kk
c
c initialization.
c
      n1 = i1
      n2 = i2
      n3 = i3
c
c add k as a neighbor of i1, i2, and i3.
c
      lp = lstptr(lend(n1),n2,list,lptr)
      call insert (k,lp,list,lptr,lnew)
      lp = lstptr(lend(n2),n3,list,lptr)
      call insert (k,lp,list,lptr,lnew)
      lp = lstptr(lend(n3),n1,list,lptr)
      call insert (k,lp,list,lptr,lnew)
c
c add i1, i2, and i3 as neighbors of k.
c
      list(lnew) = n1
      list(lnew+1) = n2
      list(lnew+2) = n3
      lptr(lnew) = lnew + 1
      lptr(lnew+1) = lnew + 2
      lptr(lnew+2) = lnew
      lend(k) = lnew + 2
      lnew = lnew + 3
      return
      end
      logical function intsec (x1,y1,x2,y2,x3,y3,x4,y4)
      real(8) x1, y1, x2, y2, x3, y3, x4, y4
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   09/01/88
c
c   given a pair of line segments p1-p2 and p3-p4, this
c function returns the value .true. if and only if p1-p2
c shares one or more points with p3-p4.  the line segments
c include their endpoints, and the four points need not be
c distinct.  thus, either line segment may consist of a
c single point, and the segments may meet in a v (which is
c treated as an intersection).  note that an incorrect
c decision may result from floating point error if the four
c endpoints are nearly collinear.
c
c
c on input:
c
c       x1,y1 = coordinates of p1.
c
c       x2,y2 = coordinates of p2.
c
c       x3,y3 = coordinates of p3.
c
c       x4,y4 = coordinates of p4.
c
c input parameters are not altered by this function.
c
c on output:
c
c       intsec = logical value defined above.
c
c modules required by intsec:  none
c
c***********************************************************
c
      real(8) a, b, d, dx12, dx31, dx34, dy12, dy31, dy34
c
c test for overlap between the smallest rectangles that
c   contain the line segments and have sides parallel to
c   the axes.
c
      if ((x1 .lt. x3  .and.  x1 .lt. x4  .and.  x2 .lt. x3
     .     .and.  x2 .lt. x4)  .or.
     .    (x1 .gt. x3  .and.  x1 .gt. x4  .and.  x2 .gt. x3
     .     .and.  x2 .gt. x4)  .or.
     .    (y1 .lt. y3  .and.  y1 .lt. y4  .and.  y2 .lt. y3
     .     .and.  y2 .lt. y4)  .or.
     .    (y1 .gt. y3  .and.  y1 .gt. y4  .and.  y2 .gt. y3
     .     .and.  y2 .gt. y4)) then
        intsec = .false.
        return
      endif
c
c compute a = p4-p3 x p1-p3, b = p2-p1 x p1-p3, and
c   d = p2-p1 x p4-p3 (z components).
c
      dx12 = x2 - x1
      dy12 = y2 - y1
      dx34 = x4 - x3
      dy34 = y4 - y3
      dx31 = x1 - x3
      dy31 = y1 - y3
      a = dx34*dy31 - dx31*dy34
      b = dx12*dy31 - dx31*dy12
      d = dx12*dy34 - dx34*dy12
      if (d .eq. 0.) go to 1
c
c d .ne. 0 and the point of intersection of the lines de-
c   fined by the line segments is p = p1 + (a/d)*(p2-p1) =
c   p3 + (b/d)*(p4-p3).
c
      intsec = a/d .ge. 0.  .and.  a/d .le. 1.  .and.
     .         b/d .ge. 0.  .and.  b/d .le. 1.
      return
c
c d .eq. 0 and thus either the line segments are parallel,
c   or one (or both) of them is a single point.
c
    1 intsec = a .eq. 0.  .and.  b .eq. 0.
      return
      end
      integer function jrand (n, ix,iy,iz )
      integer n, ix, iy, iz
c
c***********************************************************
c
c                                              from stripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   07/28/98
c
c   this function returns a uniformly distributed pseudo-
c random integer in the range 1 to n.
c
c
c on input:
c
c       n = maximum value to be returned.
c
c n is not altered by this function.
c
c       ix,iy,iz = integer seeds initialized to values in
c                  the range 1 to 30,000 before the first
c                  call to jrand, and not altered between
c                  subsequent calls (unless a sequence of
c                  random numbers is to be repeated by
c                  reinitializing the seeds).
c
c on output:
c
c       ix,iy,iz = updated integer seeds.
c
c       jrand = random integer in the range 1 to n.
c
c reference:  b. a. wichmann and i. d. hill, "an efficient
c             and portable pseudo-random number generator",
c             applied statistics, vol. 31, no. 2, 1982,
c             pp. 188-190.
c
c modules required by jrand:  none
c
c intrinsic functions called by jrand:  int, mod, real
c
c***********************************************************
c
      real(8) u, x
c
c local parameters:
c
c u = pseudo-random number uniformly distributed in the
c     interval (0,1).
c x = pseudo-random number in the range 0 to 3 whose frac-
c       tional part is u.
c
      ix = mod(171*ix,30269)
      iy = mod(172*iy,30307)
      iz = mod(170*iz,30323)
      x = (real(ix)/30269.) + (real(iy)/30307.) +
     .    (real(iz)/30323.)
      u = x - int(x)
      jrand = real(n)*u + 1.
      return
      end
      logical function left (x1,y1,x2,y2,x0,y0)
      real(8)    x1, y1, x2, y2, x0, y0
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   09/01/88
c
c   this function determines whether node n0 is to the left
c or to the right of the line through n1-n2 as viewed by an
c observer at n1 facing n2.
c
c
c on input:
c
c       x1,y1 = coordinates of n1.
c
c       x2,y2 = coordinates of n2.
c
c       x0,y0 = coordinates of n0.
c
c input parameters are not altered by this function.
c
c on output:
c
c       left = .true. if and only if (x0,y0) is on or to the
c              left of the directed line n1->n2.
c
c modules required by left:  none
c
c***********************************************************
c
      real(8) dx1, dy1, dx2, dy2
c
c local parameters:
c
c dx1,dy1 = x,y components of the vector n1->n2
c dx2,dy2 = x,y components of the vector n1->n0
c
      dx1 = x2-x1
      dy1 = y2-y1
      dx2 = x0-x1
      dy2 = y0-y1
c
c if the sign of the vector cross product of n1->n2 and
c   n1->n0 is positive, then sin(a) > 0, where a is the
c   angle between the vectors, and thus a is in the range
c   (0,180) degrees.
c
      left = dx1*dy2 .ge. dx2*dy1
      return
      end
      integer function lstptr (lpl,nb,list,lptr)
      integer lpl, nb, list(*), lptr(*)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   09/01/88
c
c   this function returns the index (list pointer) of nb in
c the adjacency list for n0, where lpl = lend(n0).
c
c
c on input:
c
c       lpl = lend(n0)
c
c       nb = index of the node whose pointer is to be re-
c            turned.  nb must be connected to n0.
c
c       list,lptr = data structure defining the triangula-
c                   tion.  refer to subroutine trmesh.
c
c input parameters are not altered by this function.
c
c on output:
c
c       lstptr = pointer such that list(lstptr) = nb or
c                list(lstptr) = -nb, unless nb is not a
c                neighbor of n0, in which case lstptr = lpl.
c
c modules required by lstptr:  none
c
c***********************************************************
c
      integer lp, nd
c
      lp = lptr(lpl)
    1 nd = list(lp)
        if (nd .eq. nb) go to 2
        lp = lptr(lp)
        if (lp .ne. lpl) go to 1
c
    2 lstptr = lp
      return
      end
      integer function nbcnt (lpl,lptr)
      integer lpl, lptr(*)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   09/01/88
c
c   this function returns the number of neighbors of a node
c n0 in a triangulation created by subroutine trmesh (or
c trmshr).
c
c
c on input:
c
c       lpl = list pointer to the last neighbor of n0 --
c             lpl = lend(n0).
c
c       lptr = array of pointers associated with list.
c
c input parameters are not altered by this function.
c
c on output:
c
c       nbcnt = number of neighbors of n0.
c
c modules required by nbcnt:  none
c
c***********************************************************
c
      integer k, lp
c
      lp = lpl
      k = 1
c
    1 lp = lptr(lp)
        if (lp .eq. lpl) go to 2
        k = k + 1
        go to 1
c
    2 nbcnt = k
      return
      end
      integer function nearnd (xp,yp,ist,n,x,y,list,lptr,
     .                         lend, dsq)
      integer ist, n, list(*), lptr(*), lend(n)
      real(8)    xp, yp, x(n), y(n), dsq
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   06/27/98
c
c   given a point p in the plane and a delaunay triangula-
c tion created by subroutine trmesh or trmshr, this function
c returns the index of the nearest triangulation node to p.
c
c   the algorithm consists of implicitly adding p to the
c triangulation, finding the nearest neighbor to p, and
c implicitly deleting p from the triangulation.  thus, it
c is based on the fact that, if p is a node in a delaunay
c triangulation, the nearest node to p is a neighbor of p.
c
c
c on input:
c
c       xp,yp = cartesian coordinates of the point p to be
c               located relative to the triangulation.
c
c       ist = index of a node at which trfind begins the
c             search.  search time depends on the proximity
c             of this node to p.
c
c       n = number of nodes in the triangulation.  n .ge. 3.
c
c       x,y = arrays of length n containing the cartesian
c             coordinates of the nodes.
c
c       list,lptr,lend = data structure defining the trian-
c                        gulation.  refer to trmesh.
c
c input parameters are not altered by this function.
c
c on output:
c
c       nearnd = nodal index of the nearest node to p, or 0
c                if n < 3 or the triangulation data struc-
c                ture is invalid.
c
c       dsq = squared distance between p and nearnd unless
c             nearnd = 0.
c
c       note that the number of candidates for nearnd
c       (neighbors of p) is limited to lmax defined in
c       the parameter statement below.
c
c modules required by nearnd:  jrand, left, lstptr, trfind
c
c intrinsic function called by nearnd:  abs
c
c***********************************************************
c
      integer   lstptr
      integer   lmax
      parameter (lmax=25)
      integer   i1, i2, i3, l, listp(lmax), lp, lp1, lp2,
     .          lpl, lptrp(lmax), n1, n2, n3, nr, nst
      real(8)      cos1, cos2, ds1, dsr, dx11, dx12, dx21,
     .          dx22, dy11, dy12, dy21, dy22, sin1, sin2
c
c store local parameters and test for n invalid.
c
      if (n .lt. 3) go to 7
      nst = ist
      if (nst .lt. 1  .or.  nst .gt. n) nst = 1
c
c find a triangle (i1,i2,i3) containing p, or the rightmost
c   (i1) and leftmost (i2) visible boundary nodes as viewed
c   from p.
c
      call trfind (nst,xp,yp,n,x,y,list,lptr,lend, i1,i2,i3)
c
c test for collinear nodes.
c
      if (i1 .eq. 0) go to 7
c
c store the linked list of 'neighbors' of p in listp and
c   lptrp.  i1 is the first neighbor, and 0 is stored as
c   the last neighbor if p is not contained in a triangle.
c   l is the length of listp and lptrp, and is limited to
c   lmax.
c
      if (i3 .ne. 0) then
        listp(1) = i1
        lptrp(1) = 2
        listp(2) = i2
        lptrp(2) = 3
        listp(3) = i3
        lptrp(3) = 1
        l = 3
      else
        n1 = i1
        l = 1
        lp1 = 2
        listp(l) = n1
        lptrp(l) = lp1
c
c   loop on the ordered sequence of visible boundary nodes
c     n1 from i1 to i2.
c
    1   lpl = lend(n1)
          n1 = -list(lpl)
          l = lp1
          lp1 = l+1
          listp(l) = n1
          lptrp(l) = lp1
          if (n1 .ne. i2  .and.  lp1 .lt. lmax) go to 1
        l = lp1
        listp(l) = 0
        lptrp(l) = 1
      endif
c
c initialize variables for a loop on arcs n1-n2 opposite p
c   in which new 'neighbors' are 'swapped' in.  n1 follows
c   n2 as a neighbor of p, and lp1 and lp2 are the listp
c   indexes of n1 and n2.
c
      lp2 = 1
      n2 = i1
      lp1 = lptrp(1)
      n1 = listp(lp1)
c
c begin loop:  find the node n3 opposite n1->n2.
c
    2 lp = lstptr(lend(n1),n2,list,lptr)
        if (list(lp) .lt. 0) go to 4
        lp = lptr(lp)
        n3 = abs(list(lp))
c
c swap test:  exit the loop if l = lmax.
c
        if (l .eq. lmax) go to 5
        dx11 = x(n1) - x(n3)
        dx12 = x(n2) - x(n3)
        dx22 = x(n2) - xp
        dx21 = x(n1) - xp
c
        dy11 = y(n1) - y(n3)
        dy12 = y(n2) - y(n3)
        dy22 = y(n2) - yp
        dy21 = y(n1) - yp
c
        cos1 = dx11*dx12 + dy11*dy12
        cos2 = dx22*dx21 + dy22*dy21
        if (cos1 .ge. 0.  .and.  cos2 .ge. 0.) go to 4
        if (cos1 .lt. 0.  .and.  cos2 .lt. 0.) go to 3
c
        sin1 = dx11*dy12 - dx12*dy11
        sin2 = dx22*dy21 - dx21*dy22
        if (sin1*cos2 + cos1*sin2 .ge. 0.) go to 4
c
c swap:  insert n3 following n2 in the adjacency list for p.
c        the two new arcs opposite p must be tested.
c
    3   l = l+1
        lptrp(lp2) = l
        listp(l) = n3
        lptrp(l) = lp1
        lp1 = l
        n1 = n3
        go to 2
c
c no swap:  advance to the next arc and test for termination
c           on n1 = i1 (lp1 = 1) or n1 followed by 0.
c
    4   if (lp1 .eq. 1) go to 5
        lp2 = lp1
        n2 = n1
        lp1 = lptrp(lp1)
        n1 = listp(lp1)
        if (n1 .eq. 0) go to 5
        go to 2
c
c set nr and dsr to the index of the nearest node to p and
c   its squared distance from p, respectively.
c
    5 nr = i1
      dsr = (x(nr)-xp)**2 + (y(nr)-yp)**2
      do 6 lp = 2,l
        n1 = listp(lp)
        if (n1 .eq. 0) go to 6
        ds1 = (x(n1)-xp)**2 + (y(n1)-yp)**2
        if (ds1 .lt. dsr) then
          nr = n1
          dsr = ds1
        endif
    6   continue
      dsq = dsr
      nearnd = nr
      return
c
c invalid input.
c
    7 nearnd = 0
      return
      end
      subroutine optim (x,y,na, list,lptr,lend,nit,iwk, ier)
      integer na, list(*), lptr(*), lend(*), nit, iwk(2,na),
     .        ier
      real(8)    x(*), y(*)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   06/27/98
c
c   given a set of na triangulation arcs, this subroutine
c optimizes the portion of the triangulation consisting of
c the quadrilaterals (pairs of adjacent triangles) which
c have the arcs as diagonals by applying the circumcircle
c test and appropriate swaps to the arcs.
c
c   an iteration consists of applying the swap test and
c swaps to all na arcs in the order in which they are
c stored.  the iteration is repeated until no swap occurs
c or nit iterations have been performed.  the bound on the
c number of iterations may be necessary to prevent an
c infinite loop caused by cycling (reversing the effect of a
c previous swap) due to floating point inaccuracy when four
c or more nodes are nearly cocircular.
c
c
c on input:
c
c       x,y = arrays containing the nodal coordinates.
c
c       na = number of arcs in the set.  na .ge. 0.
c
c the above parameters are not altered by this routine.
c
c       list,lptr,lend = data structure defining the trian-
c                        gulation.  refer to subroutine
c                        trmesh.
c
c       nit = maximum number of iterations to be performed.
c             a reasonable value is 3*na.  nit .ge. 1.
c
c       iwk = integer array dimensioned 2 by na containing
c             the nodal indexes of the arc endpoints (pairs
c             of endpoints are stored in columns).
c
c on output:
c
c       list,lptr,lend = updated triangulation data struc-
c                        ture reflecting the swaps.
c
c       nit = number of iterations performed.
c
c       iwk = endpoint indexes of the new set of arcs
c             reflecting the swaps.
c
c       ier = error indicator:
c             ier = 0 if no errors were encountered.
c             ier = 1 if a swap occurred on the last of
c                     maxit iterations, where maxit is the
c                     value of nit on input.  the new set
c                     of arcs in not necessarily optimal
c                     in this case.
c             ier = 2 if na < 0 or nit < 1 on input.
c             ier = 3 if iwk(2,i) is not a neighbor of
c                     iwk(1,i) for some i in the range 1
c                     to na.  a swap may have occurred in
c                     this case.
c             ier = 4 if a zero pointer was returned by
c                     subroutine swap.
c
c modules required by optim:  lstptr, swap, swptst
c
c intrinsic function called by optim:  abs
c
c***********************************************************
c
      logical swptst
      integer i, io1, io2, iter, lp, lp21, lpl, lpp, maxit,
     .        n1, n2, nna
      logical swp
c
c local parameters:
c
c i =       column index for iwk
c io1,io2 = nodal indexes of the endpoints of an arc in iwk
c iter =    iteration count
c lp =      list pointer
c lp21 =    parameter returned by swap (not used)
c lpl =     pointer to the last neighbor of io1
c lpp =     pointer to the node preceding io2 as a neighbor
c             of io1
c maxit =   input value of nit
c n1,n2 =   nodes opposite io1->io2 and io2->io1,
c             respectively
c nna =     local copy of na
c swp =     flag set to true iff a swap occurs in the
c             optimization loop
c
      nna = na
      maxit = nit
      if (nna .lt. 0  .or.  maxit .lt. 1) go to 7
c
c initialize iteration count iter and test for na = 0.
c
      iter = 0
      if (nna .eq. 0) go to 5
c
c top of loop --
c   swp = true iff a swap occurred in the current iteration.
c
    1 if (iter .eq. maxit) go to 6
      iter = iter + 1
      swp = .false.
c
c   inner loop on arcs io1-io2 --
c
      do 4 i = 1,nna
        io1 = iwk(1,i)
        io2 = iwk(2,i)
c
c   set n1 and n2 to the nodes opposite io1->io2 and
c     io2->io1, respectively.  determine the following:
c
c     lpl = pointer to the last neighbor of io1,
c     lp = pointer to io2 as a neighbor of io1, and
c     lpp = pointer to the node n2 preceding io2.
c
        lpl = lend(io1)
        lpp = lpl
        lp = lptr(lpp)
    2   if (list(lp) .eq. io2) go to 3
          lpp = lp
          lp = lptr(lpp)
          if (lp .ne. lpl) go to 2
c
c   io2 should be the last neighbor of io1.  test for no
c     arc and bypass the swap test if io1 is a boundary
c     node.
c
        if (abs(list(lp)) .ne. io2) go to 8
        if (list(lp) .lt. 0) go to 4
c
c   store n1 and n2, or bypass the swap test if io1 is a
c     boundary node and io2 is its first neighbor.
c
    3   n2 = list(lpp)
        if (n2 .lt. 0) go to 4
        lp = lptr(lp)
        n1 = abs(list(lp))
c
c   test io1-io2 for a swap, and update iwk if necessary.
c
        if ( .not. swptst(n1,n2,io1,io2,x,y) ) go to 4
        call swap (n1,n2,io1,io2, list,lptr,lend, lp21)
        if (lp21 .eq. 0) go to 9
        swp = .true.
        iwk(1,i) = n1
        iwk(2,i) = n2
    4   continue
      if (swp) go to 1
c
c successful termination.
c
    5 nit = iter
      ier = 0
      return
c
c maxit iterations performed without convergence.
c
    6 nit = maxit
      ier = 1
      return
c
c invalid input parameter.
c
    7 nit = 0
      ier = 2
      return
c
c io2 is not a neighbor of io1.
c
    8 nit = iter
      ier = 3
      return
c
c zero pointer returned by swap.
c
    9 nit = iter
      ier = 4
      return
      end
      real(8) function store (x)
      real(8) x
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   03/18/90
c
c   this function forces its argument x to be stored in a
c memory location, thus providing a means of determining
c floating point number characteristics (such as the machine
c precision) when it is necessary to avoid computation in
c high precision registers.
c
c
c on input:
c
c       x = value to be stored.
c
c x is not altered by this function.
c
c on output:
c
c       store = value of x after it has been stored and
c               possibly truncated or rounded to the single
c               precision word length.
c
c modules required by store:  none
c
c***********************************************************
c
      real(8) y
      common/stcom/y
c
      y = x
      store = y
      return
      end
      subroutine swap (in1,in2,io1,io2, list,lptr,
     .                 lend, lp21)
      integer in1, in2, io1, io2, list(*), lptr(*), lend(*),
     .        lp21
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   06/22/98
c
c   given a triangulation of a set of points on the unit
c sphere, this subroutine replaces a diagonal arc in a
c strictly convex quadrilateral (defined by a pair of adja-
c cent triangles) with the other diagonal.  equivalently, a
c pair of adjacent triangles is replaced by another pair
c having the same union.
c
c
c on input:
c
c       in1,in2,io1,io2 = nodal indexes of the vertices of
c                         the quadrilateral.  io1-io2 is re-
c                         placed by in1-in2.  (io1,io2,in1)
c                         and (io2,io1,in2) must be trian-
c                         gles on input.
c
c the above parameters are not altered by this routine.
c
c       list,lptr,lend = data structure defining the trian-
c                        gulation.  refer to subroutine
c                        trmesh.
c
c on output:
c
c       list,lptr,lend = data structure updated with the
c                        swap -- triangles (io1,io2,in1) and
c                        (io2,io1,in2) are replaced by
c                        (in1,in2,io2) and (in2,in1,io1)
c                        unless lp21 = 0.
c
c       lp21 = index of in1 as a neighbor of in2 after the
c              swap is performed unless in1 and in2 are
c              adjacent on input, in which case lp21 = 0.
c
c module required by swap:  lstptr
c
c intrinsic function called by swap:  abs
c
c***********************************************************
c
      integer lstptr
      integer lp, lph, lpsav
c
c local parameters:
c
c lp,lph,lpsav = list pointers
c
c
c test for in1 and in2 adjacent.
c
      lp = lstptr(lend(in1),in2,list,lptr)
      if (abs(list(lp)) .eq. in2) then
        lp21 = 0
        return
      endif
c
c delete io2 as a neighbor of io1.
c
      lp = lstptr(lend(io1),in2,list,lptr)
      lph = lptr(lp)
      lptr(lp) = lptr(lph)
c
c if io2 is the last neighbor of io1, make in2 the
c   last neighbor.
c
      if (lend(io1) .eq. lph) lend(io1) = lp
c
c insert in2 as a neighbor of in1 following io1
c   using the hole created above.
c
      lp = lstptr(lend(in1),io1,list,lptr)
      lpsav = lptr(lp)
      lptr(lp) = lph
      list(lph) = in2
      lptr(lph) = lpsav
c
c delete io1 as a neighbor of io2.
c
      lp = lstptr(lend(io2),in1,list,lptr)
      lph = lptr(lp)
      lptr(lp) = lptr(lph)
c
c if io1 is the last neighbor of io2, make in1 the
c   last neighbor.
c
      if (lend(io2) .eq. lph) lend(io2) = lp
c
c insert in1 as a neighbor of in2 following io2.
c
      lp = lstptr(lend(in2),io2,list,lptr)
      lpsav = lptr(lp)
      lptr(lp) = lph
      list(lph) = in1
      lptr(lph) = lpsav
      lp21 = lph
      return
      end
      logical function swptst (in1,in2,io1,io2,x,y)
      integer in1, in2, io1, io2
      real(8)    x(*), y(*)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   09/01/88
c
c   this function applies the circumcircle test to a quadri-
c lateral defined by a pair of adjacent triangles.  the
c diagonal arc (shared triangle side) should be swapped for
c the other diagonl if and only if the fourth vertex is
c strictly interior to the circumcircle of one of the
c triangles (the decision is independent of the choice of
c triangle).  equivalently, the diagonal is chosen to maxi-
c mize the smallest of the six interior angles over the two
c pairs of possible triangles (the decision is for no swap
c if the quadrilateral is not strictly convex).
c
c   when the four vertices are nearly cocircular (the
c neutral case), the preferred decision is no swap -- in
c order to avoid unnecessary swaps and, more important, to
c avoid cycling in subroutine optim which is called by
c delnod and edge.  thus, a tolerance swtol (stored in
c swpcom by trmesh or trmshr) is used to define 'nearness'
c to the neutral case.
c
c
c on input:
c
c       in1,in2,io1,io2 = nodal indexes of the vertices of
c                         the quadrilateral.  io1-io2 is the
c                         triangulation arc (shared triangle
c                         side) to be replaced by in1-in2 if
c                         the decision is to swap.  the
c                         triples (io1,io2,in1) and (io2,
c                         io1,in2) must define triangles (be
c                         in counterclockwise order) on in-
c                         put.
c
c       x,y = arrays containing the nodal coordinates.
c
c input parameters are not altered by this routine.
c
c on output:
c
c       swptst = .true. if and only if the arc connecting
c                io1 and io2 is to be replaced.
c
c modules required by swptst:  none
c
c***********************************************************
c
      real(8) dx11, dx12, dx22, dx21, dy11, dy12, dy22, dy21,
     .     sin1, sin2, cos1, cos2, sin12, swtol
c
c tolerance stored by trmesh or trmshr.
c
      common/swpcom/swtol
c
c local parameters:
c
c dx11,dy11 = x,y components of the vector in1->io1
c dx12,dy12 = x,y components of the vector in1->io2
c dx22,dy22 = x,y components of the vector in2->io2
c dx21,dy21 = x,y components of the vector in2->io1
c sin1 =      cross product of the vectors in1->io1 and
c               in1->io2 -- proportional to sin(t1), where
c               t1 is the angle at in1 formed by the vectors
c cos1 =      inner product of the vectors in1->io1 and
c               in1->io2 -- proportional to cos(t1)
c sin2 =      cross product of the vectors in2->io2 and
c               in2->io1 -- proportional to sin(t2), where
c               t2 is the angle at in2 formed by the vectors
c cos2 =      inner product of the vectors in2->io2 and
c               in2->io1 -- proportional to cos(t2)
c sin12 =     sin1*cos2 + cos1*sin2 -- proportional to
c               sin(t1+t2)
c
c
c compute the vectors containing the angles t1 and t2.
c
      dx11 = x(io1) - x(in1)
      dx12 = x(io2) - x(in1)
      dx22 = x(io2) - x(in2)
      dx21 = x(io1) - x(in2)
c
      dy11 = y(io1) - y(in1)
      dy12 = y(io2) - y(in1)
      dy22 = y(io2) - y(in2)
      dy21 = y(io1) - y(in2)
c
c compute inner products.
c
      cos1 = dx11*dx12 + dy11*dy12
      cos2 = dx22*dx21 + dy22*dy21
c
c the diagonals should be swapped iff (t1+t2) > 180
c   degrees.  the following two tests ensure numerical
c   stability:  the decision must be false when both
c   angles are close to 0, and true when both angles
c   are close to 180 degrees.
c
      if (cos1 .ge. 0.  .and.  cos2 .ge. 0.) go to 2
      if (cos1 .lt. 0.  .and.  cos2 .lt. 0.) go to 1
c
c compute vector cross products (z-components).
c
      sin1 = dx11*dy12 - dx12*dy11
      sin2 = dx22*dy21 - dx21*dy22
      sin12 = sin1*cos2 + cos1*sin2
      if (sin12 .ge. -swtol) go to 2
c
c swap.
c
    1 swptst = .true.
      return
c
c no swap.
c
    2 swptst = .false.
      return
      end
      subroutine trfind (nst,px,py,n,x,y,list,lptr,lend, i1,
     .                   i2,i3)
      integer nst, n, list(*), lptr(*), lend(n), i1, i2, i3
      real(8)    px, py, x(n), y(n)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   07/28/98
c
c   this subroutine locates a point p relative to a triangu-
c lation created by subroutine trmesh or trmshr.  if p is
c contained in a triangle, the three vertex indexes are
c returned.  otherwise, the indexes of the rightmost and
c leftmost visible boundary nodes are returned.
c
c
c on input:
c
c       nst = index of a node at which trfind begins the
c             search.  search time depends on the proximity
c             of this node to p.
c
c       px,py = x and y coordinates of the point p to be
c               located.
c
c       n = number of nodes in the triangulation.  n .ge. 3.
c
c       x,y = arrays of length n containing the coordinates
c             of the nodes in the triangulation.
c
c       list,lptr,lend = data structure defining the trian-
c                        gulation.  refer to subroutine
c                        trmesh.
c
c input parameters are not altered by this routine.
c
c on output:
c
c       i1,i2,i3 = nodal indexes, in counterclockwise order,
c                  of the vertices of a triangle containing
c                  p if p is contained in a triangle.  if p
c                  is not in the convex hull of the nodes,
c                  i1 indexes the rightmost visible boundary
c                  node, i2 indexes the leftmost visible
c                  boundary node, and i3 = 0.  rightmost and
c                  leftmost are defined from the perspective
c                  of p, and a pair of points are visible
c                  from each other if and only if the line
c                  segment joining them intersects no trian-
c                  gulation arc.  if p and all of the nodes
c                  lie on a common line, then i1 = i2 = i3 =
c                  0 on output.
c
c modules required by trfind:  jrand, left, lstptr, store
c
c intrinsic function called by trfind:  abs
c
c***********************************************************
c
      integer jrand, lstptr
      logical left
      real(8)    store
      integer ix, iy, iz, lp, n0, n1, n1s, n2, n2s, n3, n4,
     .        nb, nf, nl, np, npp
      logical frwrd
      real(8)    b1, b2, xa, xb, xc, xp, ya, yb, yc, yp
c
      save    ix, iy, iz
      data    ix/1/, iy/2/, iz/3/
c
c local parameters:
c
c b1,b2 =    unnormalized barycentric coordinates of p with
c              respect to (n1,n2,n3)
c ix,iy,iz = integer seeds for jrand
c lp =       list pointer
c n0,n1,n2 = nodes in counterclockwise order defining a
c              cone (with vertex n0) containing p
c n1s,n2s =  saved values of n1 and n2
c n3,n4 =    nodes opposite n1->n2 and n2->n1, respectively
c nb =       index of a boundary node -- first neighbor of
c              nf or last neighbor of nl in the boundary
c              traversal loops
c nf,nl =    first and last neighbors of n0, or first
c              (rightmost) and last (leftmost) nodes
c              visible from p when p is exterior to the
c              triangulation
c np,npp =   indexes of boundary nodes used in the boundary
c              traversal loops
c xa,xb,xc = dummy arguments for frwrd
c ya,yb,yc = dummy arguments for frwrd
c xp,yp =    local variables containing the components of p
c
c statement function:
c
c frwrd = true iff c is forward of a->b
c              iff <a->b,a->c> .ge. 0.
c
      frwrd(xa,ya,xb,yb,xc,yc) = (xb-xa)*(xc-xa) +
     .                           (yb-ya)*(yc-ya) .ge. 0.
c
c initialize variables.
c
      xp = px
      yp = py
      n0 = nst
      if (n0 .lt. 1  .or.  n0 .gt. n)
     .  n0 = jrand(n, ix,iy,iz )
c
c set nf and nl to the first and last neighbors of n0, and
c   initialize n1 = nf.
c
    1 lp = lend(n0)
      nl = list(lp)
      lp = lptr(lp)
      nf = list(lp)
      n1 = nf
c
c find a pair of adjacent neighbors n1,n2 of n0 that define
c   a wedge containing p:  p left n0->n1 and p right n0->n2.
c
      if (nl .gt. 0) go to 2
c
c   n0 is a boundary node.  test for p exterior.
c
      nl = -nl
      if ( .not. left(x(n0),y(n0),x(nf),y(nf),xp,yp) ) then
        nl = n0
        go to 9
      endif
      if ( .not. left(x(nl),y(nl),x(n0),y(n0),xp,yp) ) then
        nb = nf
        nf = n0
        np = nl
        npp = n0
        go to 11
      endif
      go to 3
c
c   n0 is an interior node.  find n1.
c
    2 if ( left(x(n0),y(n0),x(n1),y(n1),xp,yp) ) go to 3
        lp = lptr(lp)
        n1 = list(lp)
        if (n1 .eq. nl) go to 6
        go to 2
c
c   p is to the left of edge n0->n1.  initialize n2 to the
c     next neighbor of n0.
c
    3 lp = lptr(lp)
        n2 = abs(list(lp))
        if ( .not. left(x(n0),y(n0),x(n2),y(n2),xp,yp) )
     .    go to 7
        n1 = n2
        if (n1 .ne. nl) go to 3
      if ( .not. left(x(n0),y(n0),x(nf),y(nf),xp,yp) )
     .  go to 6
      if (xp .eq. x(n0) .and. yp .eq. y(n0)) go to 5
c
c   p is left of or on edges n0->nb for all neighbors nb
c     of n0.
c   all points are collinear iff p is left of nb->n0 for
c     all neighbors nb of n0.  search the neighbors of n0.
c     note -- n1 = nl and lp points to nl.
c
    4 if ( .not. left(x(n1),y(n1),x(n0),y(n0),xp,yp) )
     .  go to 5
        lp = lptr(lp)
        n1 = abs(list(lp))
        if (n1 .eq. nl) go to 17
        go to 4
c
c   p is to the right of n1->n0, or p=n0.  set n0 to n1 and
c     start over.
c
    5 n0 = n1
      go to 1
c
c   p is between edges n0->n1 and n0->nf.
c
    6 n2 = nf
c
c p is contained in the wedge defined by line segments
c   n0->n1 and n0->n2, where n1 is adjacent to n2.  set
c   n3 to the node opposite n1->n2, and save n1 and n2 to
c   test for cycling.
c
    7 n3 = n0
      n1s = n1
      n2s = n2
c
c top of edge hopping loop.  test for termination.
c
    8 if ( left(x(n1),y(n1),x(n2),y(n2),xp,yp) ) then
c
c   p left n1->n2 and hence p is in (n1,n2,n3) unless an
c     error resulted from floating point inaccuracy and
c     collinearity.  compute the unnormalized barycentric
c     coordinates of p with respect to (n1,n2,n3).
c
        b1 = (x(n3)-x(n2))*(yp-y(n2)) -
     .       (xp-x(n2))*(y(n3)-y(n2))
        b2 = (x(n1)-x(n3))*(yp-y(n3)) -
     .       (xp-x(n3))*(y(n1)-y(n3))
        if (store(b1+1.) .ge. 1.  .and.
     .      store(b2+1.) .ge. 1.) go to 16
c
c   restart with n0 randomly selected.
c
        n0 = jrand(n, ix,iy,iz )
        go to 1
      endif
c
c   set n4 to the neighbor of n2 which follows n1 (node
c     opposite n2->n1) unless n1->n2 is a boundary edge.
c
      lp = lstptr(lend(n2),n1,list,lptr)
      if (list(lp) .lt. 0) then
        nf = n2
        nl = n1
        go to 9
      endif
      lp = lptr(lp)
      n4 = abs(list(lp))
c
c   select the new edge n1->n2 which intersects the line
c     segment n0-p, and set n3 to the node opposite n1->n2.
c
      if ( left(x(n0),y(n0),x(n4),y(n4),xp,yp) ) then
        n3 = n1
        n1 = n4
        n2s = n2
        if (n1 .ne. n1s  .and.  n1 .ne. n0) go to 8
      else
        n3 = n2
        n2 = n4
        n1s = n1
        if (n2 .ne. n2s  .and.  n2 .ne. n0) go to 8
      endif
c
c   the starting node n0 or edge n1-n2 was encountered
c     again, implying a cycle (infinite loop).  restart
c     with n0 randomly selected.
c
      n0 = jrand(n, ix,iy,iz )
      go to 1
c
c boundary traversal loops.  nl->nf is a boundary edge and
c   p right nl->nf.  save nl and nf.

    9 np = nl
      npp = nf
c
c find the first (rightmost) visible boundary node nf.  nb
c   is set to the first neighbor of nf, and np is the last
c   neighbor.
c
   10 lp = lend(nf)
      lp = lptr(lp)
      nb = list(lp)
      if ( .not. left(x(nf),y(nf),x(nb),y(nb),xp,yp) )
     .  go to 12
c
c   p left nf->nb and thus nb is not visible unless an error
c     resulted from floating point inaccuracy and collinear-
c     ity of the 4 points np, nf, nb, and p.
c
   11 if ( frwrd(x(nf),y(nf),x(np),y(np),xp,yp)  .or.
     .     frwrd(x(nf),y(nf),x(np),y(np),x(nb),y(nb)) ) then
        i1 = nf
        go to 13
      endif
c
c   bottom of loop.
c
   12 np = nf
      nf = nb
      go to 10
c
c find the last (leftmost) visible boundary node nl.  nb
c   is set to the last neighbor of nl, and npp is the first
c   neighbor.
c
   13 lp = lend(nl)
      nb = -list(lp)
      if ( .not. left(x(nb),y(nb),x(nl),y(nl),xp,yp) )
     .  go to 14
c
c   p left nb->nl and thus nb is not visible unless an error
c     resulted from floating point inaccuracy and collinear-
c     ity of the 4 points p, nb, nl, and npp.
c
      if ( frwrd(x(nl),y(nl),x(npp),y(npp),xp,yp)  .or.
     .     frwrd(x(nl),y(nl),x(npp),y(npp),x(nb),y(nb)) )
     .  go to 15
c
c   bottom of loop.
c
   14 npp = nl
      nl = nb
      go to 13
c
c nl is the leftmost visible boundary node.
c
   15 i2 = nl
      i3 = 0
      return
c
c p is in the triangle (n1,n2,n3).
c
   16 i1 = n1
      i2 = n2
      i3 = n3
      return
c
c all points are collinear.
c
   17 i1 = 0
      i2 = 0
      i3 = 0
      return
      end
      subroutine trlist (ncc,lcc,n,list,lptr,lend,nrow, nt,
     .                   ltri,lct,ier)
      integer ncc, lcc(*), n, list(*), lptr(*), lend(n),
     .        nrow, nt, ltri(nrow,*), lct(*), ier
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   03/22/97
c
c   this subroutine converts a triangulation data structure
c from the linked list created by subroutine trmesh or
c trmshr to a triangle list.
c
c on input:
c
c       ncc = number of constraints.  ncc .ge. 0.
c
c       lcc = list of constraint curve starting indexes (or
c             dummy array of length 1 if ncc = 0).  refer to
c             subroutine addcst.
c
c       n = number of nodes in the triangulation.  n .ge. 3.
c
c       list,lptr,lend = linked list data structure defin-
c                        ing the triangulation.  refer to
c                        subroutine trmesh.
c
c       nrow = number of rows (entries per triangle) re-
c              served for the triangle list ltri.  the value
c              must be 6 if only the vertex indexes and
c              neighboring triangle indexes are to be
c              stored, or 9 if arc indexes are also to be
c              assigned and stored.  refer to ltri.
c
c the above parameters are not altered by this routine.
c
c       ltri = integer array of length at least nrow*nt,
c              where nt is at most 2n-5.  (a sufficient
c              length is 12n if nrow=6 or 18n if nrow=9.)
c
c       lct = integer array of length ncc or dummy array of
c             length 1 if ncc = 0.
c
c on output:
c
c       nt = number of triangles in the triangulation unless
c            ier .ne. 0, in which case nt = 0.  nt = 2n - nb
c            - 2, where nb is the number of boundary nodes.
c
c       ltri = nrow by nt array whose j-th column contains
c              the vertex nodal indexes (first three rows),
c              neighboring triangle indexes (second three
c              rows), and, if nrow = 9, arc indexes (last
c              three rows) associated with triangle j for
c              j = 1,...,nt.  the vertices are ordered
c              counterclockwise with the first vertex taken
c              to be the one with smallest index.  thus,
c              ltri(2,j) and ltri(3,j) are larger than
c              ltri(1,j) and index adjacent neighbors of
c              node ltri(1,j).  for i = 1,2,3, ltri(i+3,j)
c              and ltri(i+6,j) index the triangle and arc,
c              respectively, which are opposite (not shared
c              by) node ltri(i,j), with ltri(i+3,j) = 0 if
c              ltri(i+6,j) indexes a boundary arc.  vertex
c              indexes range from 1 to n, triangle indexes
c              from 0 to nt, and, if included, arc indexes
c              from 1 to na = nt+n-1.  the triangles are or-
c              dered on first (smallest) vertex indexes,
c              except that the sets of constraint triangles
c              (triangles contained in the closure of a con-
c              straint region) follow the non-constraint
c              triangles.
c
c       lct = array of length ncc containing the triangle
c             index of the first triangle of constraint j in
c             lct(j).  thus, the number of non-constraint
c             triangles is lct(1)-1, and constraint j con-
c             tains lct(j+1)-lct(j) triangles, where
c             lct(ncc+1) = nt+1.
c
c       ier = error indicator.
c             ier = 0 if no errors were encountered.
c             ier = 1 if ncc, n, nrow, or an lcc entry is
c                     outside its valid range on input.
c             ier = 2 if the triangulation data structure
c                     (list,lptr,lend) is invalid.  note,
c                     however, that these arrays are not
c                     completely tested for validity.
c
c modules required by trlist:  none
c
c intrinsic function called by trlist:  abs
c
c***********************************************************
c
      integer i, i1, i2, i3, isv, j, jlast, ka, kn, kt, l,
     .        lcc1, lp, lp2, lpl, lpln1, n1, n1st, n2, n3,
     .        nm2, nn
      logical arcs, cstri, pass2
c
c test for invalid input parameters and store the index
c   lcc1 of the first constraint node (if any).
c
      nn = n
      if (ncc .lt. 0  .or.  (nrow .ne. 6  .and.
     .    nrow .ne. 9)) go to 12
      lcc1 = nn+1
      if (ncc .eq. 0) then
        if (nn .lt. 3) go to 12
      else
        do 1 i = ncc,1,-1
          if (lcc1-lcc(i) .lt. 3) go to 12
          lcc1 = lcc(i)
    1     continue
        if (lcc1 .lt. 1) go to 12
      endif
c
c initialize parameters for loop on triangles kt = (n1,n2,
c   n3), where n1 < n2 and n1 < n3.  this requires two
c   passes through the nodes with all non-constraint
c   triangles stored on the first pass, and the constraint
c   triangles stored on the second.
c
c   arcs = true iff arc indexes are to be stored.
c   ka,kt = numbers of currently stored arcs and triangles.
c   n1st = starting index for the loop on nodes (n1st = 1 on
c            pass 1, and n1st = lcc1 on pass 2).
c   nm2 = upper bound on candidates for n1.
c   pass2 = true iff constraint triangles are to be stored.
c
      arcs = nrow .eq. 9
      ka = 0
      kt = 0
      n1st = 1
      nm2 = nn-2
      pass2 = .false.
c
c loop on nodes n1:  j = constraint containing n1,
c                    jlast = last node in constraint j.
c
    2 j = 0
      jlast = lcc1 - 1
      do 11 n1 = n1st,nm2
        if (n1 .gt. jlast) then
c
c n1 is the first node in constraint j+1.  update j and
c   jlast, and store the first constraint triangle index
c   if in pass 2.
c
          j = j + 1
          if (j .lt. ncc) then
            jlast = lcc(j+1) - 1
          else
            jlast = nn
          endif
          if (pass2) lct(j) = kt + 1
        endif
c
c loop on pairs of adjacent neighbors (n2,n3).  lpln1 points
c   to the last neighbor of n1, and lp2 points to n2.
c
        lpln1 = lend(n1)
        lp2 = lpln1
    3     lp2 = lptr(lp2)
          n2 = list(lp2)
          lp = lptr(lp2)
          n3 = abs(list(lp))
          if (n2 .lt. n1  .or.  n3 .lt. n1) go to 10
c
c (n1,n2,n3) is a constraint triangle iff the three nodes
c   are in the same constraint and n2 < n3.  bypass con-
c   straint triangles on pass 1 and non-constraint triangles
c   on pass 2.
c
          cstri = n1 .ge. lcc1  .and.  n2 .lt. n3  .and.
     .            n3 .le. jlast
          if ((cstri  .and.  .not. pass2)  .or.
     .        (.not. cstri  .and.  pass2)) go to 10
c
c add a new triangle kt = (n1,n2,n3).
c
          kt = kt + 1
          ltri(1,kt) = n1
          ltri(2,kt) = n2
          ltri(3,kt) = n3
c
c loop on triangle sides (i1,i2) with neighboring triangles
c   kn = (i1,i2,i3).
c
          do 9 i = 1,3
            if (i .eq. 1) then
              i1 = n3
              i2 = n2
            elseif (i .eq. 2) then
              i1 = n1
              i2 = n3
            else
              i1 = n2
              i2 = n1
            endif
c
c set i3 to the neighbor of i1 which follows i2 unless
c   i2->i1 is a boundary arc.
c
            lpl = lend(i1)
            lp = lptr(lpl)
    4       if (list(lp) .eq. i2) go to 5
              lp = lptr(lp)
              if (lp .ne. lpl) go to 4
c
c   i2 is the last neighbor of i1 unless the data structure
c     is invalid.  bypass the search for a neighboring
c     triangle if i2->i1 is a boundary arc.
c
            if (abs(list(lp)) .ne. i2) go to 13
            kn = 0
            if (list(lp) .lt. 0) go to 8
c
c   i2->i1 is not a boundary arc, and lp points to i2 as
c     a neighbor of i1.
c
    5       lp = lptr(lp)
            i3 = abs(list(lp))
c
c find l such that ltri(l,kn) = i3 (not used if kn > kt),
c   and permute the vertex indexes of kn so that i1 is
c   smallest.
c
            if (i1 .lt. i2  .and.  i1 .lt. i3) then
              l = 3
            elseif (i2 .lt. i3) then
              l = 2
              isv = i1
              i1 = i2
              i2 = i3
              i3 = isv
            else
              l = 1
              isv = i1
              i1 = i3
              i3 = i2
              i2 = isv
            endif
c
c test for kn > kt (triangle index not yet assigned).
c
            if (i1 .gt. n1  .and.  .not. pass2) go to 9
c
c find kn, if it exists, by searching the triangle list in
c   reverse order.
c
            do 6 kn = kt-1,1,-1
              if (ltri(1,kn) .eq. i1  .and.  ltri(2,kn) .eq.
     .            i2  .and.  ltri(3,kn) .eq. i3) go to 7
    6         continue
            go to 9
c
c store kt as a neighbor of kn.
c
    7       ltri(l+3,kn) = kt
c
c store kn as a neighbor of kt, and add a new arc ka.
c
    8       ltri(i+3,kt) = kn
            if (arcs) then
              ka = ka + 1
              ltri(i+6,kt) = ka
              if (kn .ne. 0) ltri(l+6,kn) = ka
            endif
    9       continue
c
c bottom of loop on triangles.
c
   10     if (lp2 .ne. lpln1) go to 3
   11     continue
c
c bottom of loop on nodes.
c
      if (.not. pass2  .and.  ncc .gt. 0) then
        pass2 = .true.
        n1st = lcc1
        go to 2
      endif
c
c no errors encountered.
c
      nt = kt
      ier = 0
      return
c
c invalid input parameter.
c
   12 nt = 0
      ier = 1
      return
c
c invalid triangulation data structure:  i1 is a neighbor of
c   i2, but i2 is not a neighbor of i1.
c
   13 nt = 0
      ier = 2
      return
      end
      subroutine trlprt (ncc,lct,n,x,y,nrow,nt,ltri,lout,
     .                   prntx)
      integer ncc, lct(*), n, nrow, nt, ltri(nrow,nt),
     .        lout
      logical prntx
      real(8)    x(n), y(n)
c
c***********************************************************
c
c                                               from trlpack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   07/02/98
c
c   given a triangulation of a set of points in the plane,
c this subroutine prints the triangle list created by
c subroutine trlist and, optionally, the nodal coordinates
c on logical unit lout.  the numbers of boundary nodes,
c triangles, and arcs, and the constraint region triangle
c indexes, if any, are also printed.
c
c   all parameters other than lout and prntx should be
c unaltered from their values on output from trlist.
c
c
c on input:
c
c       ncc = number of constraints.
c
c       lct = list of constraint triangle starting indexes
c             (or dummy array of length 1 if ncc = 0).
c
c       n = number of nodes in the triangulation.
c           3 .le. n .le. 9999.
c
c       x,y = arrays of length n containing the coordinates
c             of the nodes in the triangulation -- not used
c             unless prntx = true.
c
c       nrow = number of rows (entries per triangle) re-
c              served for the triangle list ltri.  the value
c              must be 6 if only the vertex indexes and
c              neighboring triangle indexes are stored, or 9
c              if arc indexes are also stored.
c
c       nt = number of triangles in the triangulation.
c            1 .le. nt .le. 9999.
c
c       ltri = nrow by nt array whose j-th column contains
c              the vertex nodal indexes (first three rows),
c              neighboring triangle indexes (second three
c              rows), and, if nrow = 9, arc indexes (last
c              three rows) associated with triangle j for
c              j = 1,...,nt.
c
c       lout = logical unit number for output.  0 .le. lout
c              .le. 99.  output is printed on unit 6 if lout
c              is outside its valid range on input.
c
c       prntx = logical variable with value true if and only
c               if x and y are to be printed (to 6 decimal
c               places).
c
c none of the parameters are altered by this routine.
c
c modules required by trlprt:  none
c
c***********************************************************
c
      integer i, k, lun, na, nb, nl, nlmax, nmax
      data    nmax/9999/,  nlmax/60/
c
c local parameters:
c
c   i = do-loop, nodal index, and row index for ltri
c   k = do-loop and triangle index
c   lun = logical unit number for output
c   na = number of triangulation arcs
c   nb = number of boundary nodes
c   nl = number of lines printed on the current page
c   nlmax = maximum number of print lines per page
c   nmax = maximum value of n and nt (4-digit format)
c
      lun = lout
      if (lun .lt. 0  .or.  lun .gt. 99) lun = 6
c
c print a heading and test for invalid input.
c
      write (lun,100)
      nl = 1
      if (n .lt. 3  .or.  n .gt. nmax  .or.
     .    (nrow .ne. 6  .and.  nrow .ne. 9)  .or.
     .    nt .lt. 1  .or.  nt .gt. nmax) then
c
c print an error message and bypass the loops.
c
        write (lun,110) n, nrow, nt
        go to 3
      endif
      if (prntx) then
c
c print x and y.
c
        write (lun,101)
        nl = 6
        do 1 i = 1,n
          if (nl .ge. nlmax) then
            write (lun,106)
            nl = 0
          endif
          write (lun,102) i, x(i), y(i)
          nl = nl + 1
    1     continue
      endif
c
c print the triangulation ltri.
c
      if (nl .gt. nlmax/2) then
        write (lun,106)
        nl = 0
      endif
      if (nrow .eq. 6) then
        write (lun,103)
      else
        write (lun,104)
      endif
      nl = nl + 5
      do 2 k = 1,nt
        if (nl .ge. nlmax) then
          write (lun,106)
          nl = 0
        endif
        write (lun,105) k, (ltri(i,k), i = 1,nrow)
        nl = nl + 1
    2   continue
c
c print nb, na, and nt (boundary nodes, arcs, and
c   triangles).
c
      nb = 2*n - nt - 2
      na = nt + n - 1
      if (nl .gt. nlmax-6) write (lun,106)
      write (lun,107) nb, na, nt
c
c print ncc and lct.
c
    3 write (lun,108) ncc
      if (ncc .gt. 0) write (lun,109) (lct(i), i = 1,ncc)
      return
c
c print formats:
c
  100 format (///,24x,'tripack (trlist) output')
  101 format (//16x,'node',7x,'x(node)',10x,'y(node)'//)
  102 format (16x,i4,2e17.6)
  103 format (//1x,'triangle',8x,'vertices',12x,'neighbors'/
     .        4x,'kt',7x,'n1',5x,'n2',5x,'n3',4x,'kt1',4x,
     .        'kt2',4x,'kt3'/)
  104 format (//1x,'triangle',8x,'vertices',12x,'neighbors',
     .        14x,'arcs'/
     .        4x,'kt',7x,'n1',5x,'n2',5x,'n3',4x,'kt1',4x,
     .        'kt2',4x,'kt3',4x,'ka1',4x,'ka2',4x,'ka3'/)
  105 format (2x,i4,2x,6(3x,i4),3(2x,i5))
  106 format (///)
  107 format (/1x,'nb = ',i4,' boundary nodes',5x,
     .        'na = ',i5,' arcs',5x,'nt = ',i5,
     .        ' triangles')
  108 format (/1x,'ncc =',i3,' constraint curves')
  109 format (1x,9x,14i5)
  110 format (//1x,10x,'*** invalid parameter:  n =',i5,
     .        ', nrow =',i5,', nt =',i5,' ***')
      end
      subroutine trmesh (n,x,y, list,lptr,lend,lnew,near,
     .                   next,dist,ier)
      integer n, list(*), lptr(*), lend(n), lnew, near(n),
     .        next(n), ier
      real(8)    x(n), y(n), dist(n)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   06/28/98
c
c   this subroutine creates a delaunay triangulation of a
c set of n arbitrarily distributed points in the plane re-
c ferred to as nodes.  the delaunay triangulation is defined
c as a set of triangles with the following five properties:
c
c  1)  the triangle vertices are nodes.
c  2)  no triangle contains a node other than its vertices.
c  3)  the interiors of the triangles are pairwise disjoint.
c  4)  the union of triangles is the convex hull of the set
c        of nodes (the smallest convex set which contains
c        the nodes).
c  5)  the interior of the circumcircle of each triangle
c        contains no node.
c
c the first four properties define a triangulation, and the
c last property results in a triangulation which is as close
c as possible to equiangular in a certain sense and which is
c uniquely defined unless four or more nodes lie on a common
c circle.  this property makes the triangulation well-suited
c for solving closest point problems and for triangle-based
c interpolation.
c
c   the triangulation can be generalized to a constrained
c delaunay triangulation by a call to subroutine addcst.
c this allows for user-specified boundaries defining a non-
c convex and/or multiply connected region.
c
c   the algorithm for constructing the triangulation has
c expected time complexity o(n*log(n)) for most nodal dis-
c tributions.  also, since the algorithm proceeds by adding
c nodes incrementally, the triangulation may be updated with
c the addition (or deletion) of a node very efficiently.
c the adjacency information representing the triangulation
c is stored as a linked list requiring approximately 13n
c storage locations.
c
c
c   the following is a list of the software package modules
c which a user may wish to call directly:
c
c  addcst - generalizes the delaunay triangulation to allow
c             for user-specified constraints.
c
c  addnod - updates the triangulation by appending or
c             inserting a new node.
c
c  areap  - computes the area bounded by a closed polygonal
c             curve such as the boundary of the triangula-
c             tion or of a constraint region.
c
c  bnodes - returns an array containing the indexes of the
c             boundary nodes in counterclockwise order.
c             counts of boundary nodes, triangles, and arcs
c             are also returned.
c
c  circum - computes the area, circumcenter, circumradius,
c             and, optionally, the aspect ratio of a trian-
c             gle defined by user-specified vertices.
c
c  delarc - deletes a boundary arc from the triangulation.
c
c  delnod - updates the triangulation with the deletion of a
c             node.
c
c  edge   - forces a pair of nodes to be connected by an arc
c             in the triangulation.
c
c  getnp  - determines the ordered sequence of l closest
c             nodes to a given node, along with the associ-
c             ated distances.  the distance between nodes is
c             taken to be the length of the shortest connec-
c             ting path which intersects no constraint
c             region.
c
c  intsec - determines whether or not an arbitrary pair of
c             line segments share a common point.
c
c  jrand  - generates a uniformly distributed pseudo-random
c             integer.
c
c  left   - locates a point relative to a line.
c
c  nearnd - returns the index of the nearest node to an
c             arbitrary point, along with its squared
c             distance.
c
c  store  - forces a value to be stored in main memory so
c             that the precision of floating point numbers
c             in memory locations rather than registers is
c             computed.
c
c  trlist - converts the triangulation data structure to a
c             triangle list more suitable for use in a fin-
c             ite element code.
c
c  trlprt - prints the triangle list created by subroutine
c             trlist.
c
c  trmesh - creates a delaunay triangulation of a set of
c             nodes.
c
c  trmshr - creates a delaunay triangulation (more effici-
c             ently than trmesh) of a set of nodes lying at
c             the vertices of a (possibly skewed) rectangu-
c             lar grid.
c
c  trplot - creates a level-2 encapsulated postscript (eps)
c             file containing a triangulation plot.
c
c  trprnt - prints the triangulation data structure and,
c             optionally, the nodal coordinates.
c
c
c on input:
c
c       n = number of nodes in the triangulation.  n .ge. 3.
c
c       x,y = arrays of length n containing the cartesian
c             coordinates of the nodes.  (x(k),y(k)) is re-
c             ferred to as node k, and k is referred to as
c             a nodal index.  the first three nodes must not
c             be collinear.
c
c the above parameters are not altered by this routine.
c
c       list,lptr = arrays of length at least 6n-12.
c
c       lend = array of length at least n.
c
c       near,next,dist = work space arrays of length at
c                        least n.  the space is used to
c                        efficiently determine the nearest
c                        triangulation node to each un-
c                        processed node for use by addnod.
c
c on output:
c
c       list = set of nodal indexes which, along with lptr,
c              lend, and lnew, define the triangulation as a
c              set of n adjacency lists -- counterclockwise-
c              ordered sequences of neighboring nodes such
c              that the first and last neighbors of a bound-
c              ary node are boundary nodes (the first neigh-
c              bor of an interior node is arbitrary).  in
c              order to distinguish between interior and
c              boundary nodes, the last neighbor of each
c              boundary node is represented by the negative
c              of its index.
c
c       lptr = set of pointers (list indexes) in one-to-one
c              correspondence with the elements of list.
c              list(lptr(i)) indexes the node which follows
c              list(i) in cyclical counterclockwise order
c              (the first neighbor follows the last neigh-
c              bor).
c
c       lend = set of pointers to adjacency lists.  lend(k)
c              points to the last neighbor of node k for
c              k = 1,...,n.  thus, list(lend(k)) < 0 if and
c              only if k is a boundary node.
c
c       lnew = pointer to the first empty location in list
c              and lptr (list length plus one).  list, lptr,
c              lend, and lnew are not altered if ier < 0,
c              and are incomplete if ier > 0.
c
c       near,next,dist = garbage.
c
c       ier = error indicator:
c             ier =  0 if no errors were encountered.
c             ier = -1 if n < 3 on input.
c             ier = -2 if the first three nodes are
c                      collinear.
c             ier = -4 if an error flag was returned by a
c                      call to swap in addnod.  this is an
c                      internal error and should be reported
c                      to the programmer.
c             ier =  l if nodes l and m coincide for some
c                      m > l.  the linked list represents
c                      a triangulation of nodes 1 to m-1
c                      in this case.
c
c modules required by trmesh:  addnod, bdyadd, insert,
c                                intadd, jrand, left,
c                                lstptr, store, swap,
c                                swptst, trfind
c
c intrinsic function called by trmesh:  abs
c
c***********************************************************
c
      logical left
      real(8)    store
      integer i, i0, j, k, km1, lcc(1), lp, lpl, ncc, nexti,
     .        nn
      real(8)    d, d1, d2, d3, eps, swtol
      common/swpcom/swtol
c
c local parameters:
c
c d =        squared distance from node k to node i
c d1,d2,d3 = squared distances from node k to nodes 1, 2,
c              and 3, respectively
c eps =      half the machine precision
c i,j =      nodal indexes
c i0 =       index of the node preceding i in a sequence of
c              unprocessed nodes:  i = next(i0)
c k =        index of node to be added and do-loop index:
c              k > 3
c km1 =      k-1
c lcc(1) =   dummy array
c lp =       list index (pointer) of a neighbor of k
c lpl =      pointer to the last neighbor of k
c ncc =      number of constraint curves
c nexti =    next(i)
c nn =       local copy of n
c swtol =    tolerance for function swptst
c
      nn = n
      if (nn .lt. 3) then
        ier = -1
        return
      endif
c
c compute a tolerance for function swptst:  swtol = 10*
c   (machine precision)
c
      eps = 1.
    1 eps = eps/2.
        swtol = store(eps + 1.)
        if (swtol .gt. 1.) go to 1
      swtol = eps*20.
c
c store the first triangle in the linked list.
c
      if ( .not. left(x(1),y(1),x(2),y(2),x(3),y(3)) ) then
c
c   the initial triangle is (3,2,1) = (2,1,3) = (1,3,2).
c
        list(1) = 3
        lptr(1) = 2
        list(2) = -2
        lptr(2) = 1
        lend(1) = 2
c
        list(3) = 1
        lptr(3) = 4
        list(4) = -3
        lptr(4) = 3
        lend(2) = 4
c
        list(5) = 2
        lptr(5) = 6
        list(6) = -1
        lptr(6) = 5
        lend(3) = 6
c
      elseif ( .not. left(x(2),y(2),x(1),y(1),x(3),y(3)) )
     .       then
c
c   the initial triangle is (1,2,3).
c
        list(1) = 2
        lptr(1) = 2
        list(2) = -3
        lptr(2) = 1
        lend(1) = 2
c
        list(3) = 3
        lptr(3) = 4
        list(4) = -1
        lptr(4) = 3
        lend(2) = 4
c
        list(5) = 1
        lptr(5) = 6
        list(6) = -2
        lptr(6) = 5
        lend(3) = 6
c
      else
c
c   the first three nodes are collinear.
c
        ier = -2
        return
      endif
c
c initialize lnew and test for n = 3.
c
      lnew = 7
      if (nn .eq. 3) then
        ier = 0
        return
      endif
c
c a nearest-node data structure (near, next, and dist) is
c   used to obtain an expected-time (n*log(n)) incremental
c   algorithm by enabling constant search time for locating
c   each new node in the triangulation.
c
c for each unprocessed node k, near(k) is the index of the
c   triangulation node closest to k (used as the starting
c   point for the search in subroutine trfind) and dist(k)
c   is an increasing function of the distance between nodes
c   k and near(k).
c
c since it is necessary to efficiently find the subset of
c   unprocessed nodes associated with each triangulation
c   node j (those that have j as their near entries), the
c   subsets are stored in near and next as follows:  for
c   each node j in the triangulation, i = near(j) is the
c   first unprocessed node in j's set (with i = 0 if the
c   set is empty), l = next(i) (if i > 0) is the second,
c   next(l) (if l > 0) is the third, etc.  the nodes in each
c   set are initially ordered by increasing indexes (which
c   maximizes efficiency) but that ordering is not main-
c   tained as the data structure is updated.
c
c initialize the data structure for the single triangle.
c
      near(1) = 0
      near(2) = 0
      near(3) = 0
      do 2 k = nn,4,-1
        d1 = (x(k)-x(1))**2 + (y(k)-y(1))**2
        d2 = (x(k)-x(2))**2 + (y(k)-y(2))**2
        d3 = (x(k)-x(3))**2 + (y(k)-y(3))**2
        if (d1 .le. d2  .and.  d1 .le. d3) then
          near(k) = 1
          dist(k) = d1
          next(k) = near(1)
          near(1) = k
        elseif (d2 .le. d1  .and.  d2 .le. d3) then
          near(k) = 2
          dist(k) = d2
          next(k) = near(2)
          near(2) = k
        else
          near(k) = 3
          dist(k) = d3
          next(k) = near(3)
          near(3) = k
        endif
    2   continue
c
c add the remaining nodes.  parameters for addnod are as
c   follows:
c
c   k = index of the node to be added.
c   near(k) = index of the starting node for the search in
c             trfind.
c   ncc = number of constraint curves.
c   lcc = dummy array (since ncc = 0).
c   km1 = number of nodes in the triangulation.
c
      ncc = 0
      do 7 k = 4,nn
        km1 = k-1
        call addnod (k,x(k),y(k),near(k),ncc, lcc,km1,x,y,
     .               list,lptr,lend,lnew, ier)
        if (ier .ne. 0) return
c
c remove k from the set of unprocessed nodes associated
c   with near(k).
c
        i = near(k)
        if (near(i) .eq. k) then
          near(i) = next(k)
        else
          i = near(i)
    3     i0 = i
            i = next(i0)
            if (i .ne. k) go to 3
          next(i0) = next(k)
        endif
        near(k) = 0
c
c loop on neighbors j of node k.
c
        lpl = lend(k)
        lp = lpl
    4   lp = lptr(lp)
          j = abs(list(lp))
c
c loop on elements i in the sequence of unprocessed nodes
c   associated with j:  k is a candidate for replacing j
c   as the nearest triangulation node to i.  the next value
c   of i in the sequence, next(i), must be saved before i
c   is moved because it is altered by adding i to k's set.
c
          i = near(j)
    5     if (i .eq. 0) go to 6
          nexti = next(i)
c
c test for the distance from i to k less than the distance
c   from i to j.
c
          d = (x(k)-x(i))**2 + (y(k)-y(i))**2
          if (d .lt. dist(i)) then
c
c replace j by k as the nearest triangulation node to i:
c   update near(i) and dist(i), and remove i from j's set
c   of unprocessed nodes and add it to k's set.
c
            near(i) = k
            dist(i) = d
            if (i .eq. near(j)) then
              near(j) = nexti
            else
              next(i0) = nexti
            endif
            next(i) = near(k)
            near(k) = i
          else
            i0 = i
          endif
c
c bottom of loop on i.
c
          i = nexti
          go to 5
c
c bottom of loop on neighbors j.
c
    6     if (lp .ne. lpl) go to 4
    7   continue
      return
      end
      subroutine trmshr (n,nx,x,y, nit, list,lptr,lend,lnew,
     .                   ier)
      integer  n, nx, nit, list(*), lptr(*), lend(n), lnew,
     .         ier
      real(8)     x(n), y(n)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   06/27/98
c
c   this subroutine creates a delaunay triangulation of a
c set of n nodes in the plane, where the nodes are the vert-
c ices of an nx by ny skewed rectangular grid with the
c natural ordering.  thus, n = nx*ny, and the nodes are
c ordered from left to right beginning at the top row so
c that adjacent nodes have indexes which differ by 1 in the
c x-direction and by nx in the y-direction.  a skewed rec-
c tangular grid is defined as one in which each grid cell is
c a strictly convex quadrilateral (and is thus the convex
c hull of its four vertices).  equivalently, any transfor-
c mation from a rectangle to a grid cell which is bilinear
c in both components has an invertible jacobian.
c
c   if the nodes are not distributed and ordered as defined
c above, subroutine trmesh must be called in place of this
c routine.  refer to subroutine addcst for the treatment of
c constraints.
c
c   the first phase of the algorithm consists of construc-
c ting a triangulation by choosing a diagonal arc in each
c grid cell.  if nit = 0, all diagonals connect lower left
c to upper right corners and no error checking or additional
c computation is performed.  otherwise, each diagonal arc is
c chosen to be locally optimal, and boundary arcs are added
c where necessary in order to cover the convex hull of the
c nodes.  (this is the first iteration.)  if nit > 1 and no
c error was detected, the triangulation is then optimized by
c a sequence of up to nit-1 iterations in which interior
c arcs of the triangulation are tested and swapped if appro-
c priate.  the algorithm terminates when an iteration
c results in no swaps and/or when the allowable number of
c iterations has been performed.  nit = 0 is sufficient to
c produce a delaunay triangulation if the original grid is
c actually rectangular, and nit = 1 is sufficient if it is
c close to rectangular.  note, however, that the ordering
c and distribution of nodes is not checked for validity in
c the case nit = 0, and the triangulation will not be valid
c unless the rectangular grid covers the convex hull of the
c nodes.
c
c
c on input:
c
c       n = number of nodes in the grid.  n = nx*ny for some
c           ny .ge. 2.
c
c       nx = number of grid points in the x-direction.  nx
c            .ge. 2.
c
c       x,y = arrays of length n containing coordinates of
c             the nodes with the ordering and distribution
c             defined in the header comments above.
c             (x(k),y(k)) is referred to as node k.
c
c the above parameters are not altered by this routine.
c
c       nit = nonnegative integer specifying the maximum
c             number of iterations to be employed.  refer
c             to the header comments above.
c
c       list,lptr = arrays of length at least 6n-12.
c
c       lend = array of length at least n.
c
c on output:
c
c       nit = number of iterations employed.
c
c       list,lptr,lend,lnew = data structure defining the
c                             triangulation.  refer to sub-
c                             routine trmesh.
c
c       ier = error indicator:
c             ier = 0 if no errors were encountered.
c             ier = k if the grid element with upper left
c                     corner at node k is not a strictly
c                     convex quadrilateral.  the algorithm
c                     is terminated when the first such
c                     occurrence is detected.  note that
c                     this test is not performed if nit = 0
c                     on input.
c             ier = -1 if n, nx, or nit is outside its valid
c                      range on input.
c             ier = -2 if nit > 1 on input, and the optimi-
c                      zation loop failed to converge within
c                      the allowable number of iterations.
c                      the triangulation is valid but not
c                      optimal in this case.
c
c modules required by trmshr:  insert, left, lstptr, nbcnt,
c                                store, swap, swptst
c
c intrinsic function called by trmshr:  abs
c
c***********************************************************
c
      integer lstptr, nbcnt
      logical left, swptst
      real(8)    store
      integer i, iter, j, k, kp1, lp, lpf, lpk, lpl, lpp,
     .        m1, m2, m3, m4, maxit, n0, n1, n2, n3, n4, ni,
     .        nj, nm1, nn, nnb
      logical tst
      real(8)    eps, swtol
      common/swpcom/swtol
c
c store local variables and test for errors in input
c   parameters.
c
      ni = nx
      nj = n/ni
      nn = ni*nj
      maxit = nit
      nit = 0
      if (n .ne. nn  .or.  nj .lt. 2  .or.  ni .lt. 2  .or.
     .    maxit .lt. 0) then
        ier = -1
        return
      endif
      ier = 0
c
c compute a tolerance for function swptst:  swtol = 10*
c   (machine precision)
c
      eps = 1.
    1 eps = eps/2.
        swtol = store(eps + 1.)
        if (swtol .gt. 1.) go to 1
      swtol = eps*20.
c
c loop on grid points (i,j) corresponding to nodes k =
c   (j-1)*ni + i.  tst = true iff diagonals are to be
c   chosen by the swap test.  m1, m2, m3, and m4 are the
c   slopes (-1, 0, or 1) of the diagonals in quadrants 1
c   to 4 (counterclockwise beginning with the upper right)
c   for a coordinate system with origin at node k.
c
      tst = maxit .gt. 0
      m1 = 0
      m4 = 0
      lp = 0
      kp1 = 1
      do 6 j = 1,nj
        do 5 i = 1,ni
          m2 = m1
          m3 = m4
          k = kp1
          kp1 = k + 1
          lpf = lp + 1
          if (j .eq. nj  .and.  i .ne. ni) go to 2
          if (i .ne. 1) then
            if (j .ne. 1) then
c
c   k is not in the top row, leftmost column, or bottom row
c     (unless k is the lower right corner).  take the first
c     neighbor to be the node above k.
c
              lp = lp + 1
              list(lp) = k - ni
              lptr(lp) = lp + 1
              if (m2 .le. 0) then
                lp = lp + 1
                list(lp) = k - 1 - ni
                lptr(lp) = lp + 1
              endif
            endif
c
c   k is not in the leftmost column.  the next (or first)
c     neighbor is to the left of k.
c
            lp = lp + 1
            list(lp) = k - 1
            lptr(lp) = lp + 1
            if (j .eq. nj) go to 3
            if (m3 .ge. 0) then
              lp = lp + 1
              list(lp) = k - 1 + ni
              lptr(lp) = lp + 1
            endif
          endif
c
c   k is not in the bottom row.  the next (or first)
c     neighbor is below k.
c
          lp = lp + 1
          list(lp) = k + ni
          lptr(lp) = lp + 1
c
c   test for a negative diagonal in quadrant 4 unless k is
c     in the rightmost column.  the quadrilateral associated
c     with the quadrant is tested for strict convexity un-
c     less nit = 0 on input.
c
          if (i .eq. ni) go to 3
          m4 = 1
          if (.not. tst) go to 2
          if ( left(x(kp1),y(kp1),x(k+ni),y(k+ni),x(k),y(k))
     .         .or.  left(x(k),y(k),x(kp1+ni),y(kp1+ni),
     .                    x(k+ni),y(k+ni))
     .         .or.  left(x(k+ni),y(k+ni),x(kp1),y(kp1),
     .                    x(kp1+ni),y(kp1+ni))
     .         .or.  left(x(kp1+ni),y(kp1+ni),x(k),y(k),
     .                    x(kp1),y(kp1)) )          go to 12
          if ( swptst(kp1,k+ni,k,kp1+ni,x,y) ) go to 2
          m4 = -1
          lp = lp + 1
          list(lp) = kp1 + ni
          lptr(lp) = lp + 1
c
c   the next (or first) neighbor is to the right of k.
c
    2     lp = lp + 1
          list(lp) = kp1
          lptr(lp) = lp + 1
c
c   test for a positive diagonal in quadrant 1 (the neighbor
c     of k-ni which follows k is not k+1) unless k is in the
c     top row.
c
          if (j .eq. 1) go to 3
          if (tst) then
            m1 = -1
            lpk = lstptr(lend(k-ni),k,list,lptr)
            lpk = lptr(lpk)
            if (list(lpk) .ne. kp1) then
              m1 = 1
              lp = lp + 1
              list(lp) = kp1 - ni
              lptr(lp) = lp + 1
            endif
          endif
c
c   if k is in the leftmost column (and not the top row) or
c     in the bottom row (and not the rightmost column), then
c     the next neighbor is the node above k.
c
          if (i .ne. 1  .and.  j .ne. nj) go to 4
          lp = lp + 1
          list(lp) = k - ni
          lptr(lp) = lp + 1
          if (i .eq. 1) go to 3
c
c   k is on the bottom row (and not the leftmost or right-
c     most column).
c
          if (m2 .le. 0) then
            lp = lp + 1
            list(lp) = k - 1 - ni
            lptr(lp) = lp + 1
          endif
          lp = lp + 1
          list(lp) = k - 1
          lptr(lp) = lp + 1
c
c   k is a boundary node.
c
    3     list(lp) = -list(lp)
c
c   bottom of loop.  store lend and correct lptr(lp).
c     lpf and lp point to the first and last neighbors
c     of k.
c
    4     lend(k) = lp
          lptr(lp) = lpf
    5     continue
    6   continue
c
c store lnew, and terminate the algorithm if nit = 0 on
c   input.
c
      lnew = lp + 1
      if (maxit .eq. 0) return
c
c add boundary arcs where necessary in order to cover the
c   convex hull of the nodes.  n1, n2, and n3 are consecu-
c   tive boundary nodes in counterclockwise order, and n0
c   is the starting point for each loop around the boundary.
c
      n0 = 1
      n1 = n0
      n2 = ni + 1
c
c   tst is set to true if an arc is added.  the boundary
c     loop is repeated until a traversal results in no
c     added arcs.
c
    7 tst = .false.
c
c   top of boundary loop.  set n3 to the first neighbor of
c     n2, and test for n3 left n1 -> n2.
c
    8   lpl = lend(n2)
          lp = lptr(lpl)
          n3 = list(lp)
          if ( left(x(n1),y(n1),x(n2),y(n2),x(n3),y(n3)) )
     .       n1 = n2
          if (n1 .ne. n2) then
c
c   add the boundary arc n1-n3.  if n0 = n2, the starting
c     point is changed to n3, since n2 will be removed from
c     the boundary.  n3 is inserted as the first neighbor of
c     n1, n2 is changed to an interior node, and n1 is
c     inserted as the last neighbor of n3.
c
            tst = .true.
            if (n2 .eq. n0) n0 = n3
            lp = lend(n1)
            call insert (n3,lp, list,lptr,lnew )
            list(lpl) = -list(lpl)
            lp = lend(n3)
            list(lp) = n2
            call insert (-n1,lp, list,lptr,lnew )
            lend(n3) = lnew - 1
          endif
c
c   bottom of loops.  test for termination.
c
          n2 = n3
          if (n1 .ne. n0) go to 8
        if (tst) go to 7
c
c terminate the algorithm if nit = 1 on input.
c
      nit = 1
      if (maxit .eq. 1) return
c
c optimize the triangulation by applying the swap test and
c   appropriate swaps to the interior arcs.  the loop is
c   repeated until no swaps are performed or maxit itera-
c   tions have been applied.  iter is the current iteration,
c   and tst is set to true if a swap occurs.
c
      iter = 1
      nm1 = nn - 1
    9 iter = iter + 1
        tst = .false.
c
c   loop on interior arcs n1-n2, where n2 > n1 and
c     (n1,n2,n3) and (n2,n1,n4) are adjacent triangles.
c
c   top of loop on nodes n1.
c
        do 11 n1 = 1,nm1
          lpl = lend(n1)
          n4 = list(lpl)
          lpf = lptr(lpl)
          n2 = list(lpf)
          lp = lptr(lpf)
          n3 = list(lp)
          nnb = nbcnt(lpl,lptr)
c
c   top of loop on neighbors n2 of n1.  nnb is the number of
c                                       neighbors of n1.
c
          do 10 i = 1,nnb
c
c   bypass the swap test if n1 is a boundary node and n2 is
c     the first neighbor (n4 < 0), n2 < n1, or n1-n2 is a
c     diagonal arc (already locally optimal) when iter = 2.
c
            if ( n4 .gt. 0  .and.  n2 .gt. n1  .and.
     .          (iter .ne. 2  .or.  abs(n1+ni-n2) .ne. 1) )
     .          then
              if (swptst(n3,n4,n1,n2,x,y) ) then
c
c   swap diagonal n1-n2 for n3-n4, set tst to true, and set
c     n2 to n4 (the neighbor preceding n3).
c
                call swap (n3,n4,n1,n2, list,lptr,lend, lpp)
                if (lpp .ne. 0) then
                  tst = .true.
                  n2 = n4
                endif
              endif
            endif
c
c   bottom of neighbor loop.
c
            if (list(lpl) .eq. -n3) go to 11
            n4 = n2
            n2 = n3
            lp = lstptr(lpl,n2,list,lptr)
            lp = lptr(lp)
            n3 = abs(list(lp))
   10       continue
   11     continue
c
c   test for termination.
c
        if (tst  .and.  iter .lt. maxit) go to 9
      nit = iter
      if (tst) ier = -2
      return
c
c invalid grid cell encountered.
c
   12 ier = k
      return
      end
      subroutine trplot (lun,pltsiz,wx1,wx2,wy1,wy2,ncc,lcc,
     .                   n,x,y,list,lptr,lend,title,
     .                   numbr, ier)
      character*(*) title
      integer lun, ncc, lcc(*), n, list(*), lptr(*),
     .        lend(n), ier
      logical numbr
      real(8)    pltsiz, wx1, wx2, wy1, wy2, x(n), y(n)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   07/15/98
c
c   this subroutine creates a level-2 encapsulated post-
c script (eps) file containing a triangulation plot.
c
c
c on input:
c
c       lun = logical unit number in the range 0 to 99.
c             the unit should be opened with an appropriate
c             file name before the call to this routine.
c
c       pltsiz = plot size in inches.  the window is mapped,
c                with aspect ratio preserved, to a rectangu-
c                lar viewport with maximum side-length equal
c                to .88*pltsiz (leaving room for labels out-
c                side the viewport).  the viewport is
c                centered on the 8.5 by 11 inch page, and
c                its boundary is drawn.  1.0 .le. pltsiz
c                .le. 8.5.
c
c       wx1,wx2,wy1,wy2 = parameters defining a rectangular
c                         window against which the triangu-
c                         lation is clipped.  (only the
c                         portion of the triangulation that
c                         lies in the window is drawn.)
c                         (wx1,wy1) and (wx2,wy2) are the
c                         lower left and upper right cor-
c                         ners, respectively.  wx1 < wx2 and
c                         wy1 < wy2.
c
c       ncc = number of constraint curves.  refer to subrou-
c             tine addcst.  ncc .ge. 0.
c
c       lcc = array of length ncc (or dummy parameter if
c             ncc = 0) containing the index of the first
c             node of constraint i in lcc(i).  for i = 1 to
c             ncc, lcc(i+1)-lcc(i) .ge. 3, where lcc(ncc+1)
c             = n+1.
c
c       n = number of nodes in the triangulation.  n .ge. 3.
c
c       x,y = arrays of length n containing the coordinates
c             of the nodes with non-constraint nodes in the
c             first lcc(1)-1 locations.
c
c       list,lptr,lend = data structure defining the trian-
c                        gulation.  refer to subroutine
c                        trmesh.
c
c       title = type character variable or constant contain-
c               ing a string to be centered above the plot.
c               the string must be enclosed in parentheses;
c               i.e., the first and last characters must be
c               '(' and ')', respectively, but these are not
c               displayed.  title may have at most 80 char-
c               acters including the parentheses.
c
c       numbr = option indicator:  if numbr = true, the
c               nodal indexes are plotted next to the nodes.
c
c input parameters are not altered by this routine.
c
c on output:
c
c       ier = error indicator:
c             ier = 0 if no errors were encountered.
c             ier = 1 if lun, pltsiz, ncc, or n is outside
c                     its valid range.  lcc is not tested
c                     for validity.
c             ier = 2 if wx1 >= wx2 or wy1 >= wy2.
c             ier = 3 if an error was encountered in writing
c                     to unit lun.
c
c   various plotting options can be controlled by altering
c the data statement below.
c
c modules required by trplot:  none
c
c intrinsic functions called by trplot:  abs, char, nint,
c                                          real
c
c***********************************************************
c
      integer i, ifrst, ih, ilast, ipx1, ipx2, ipy1, ipy2,
     .        iw, lp, lpl, n0, n0bak, n0for, n1, nls
      logical annot, cnstr, pass1
      real(8)    dashl, dx, dy, fsizn, fsizt, r, sfx, sfy, t,
     .        tx, ty, x0, y0
c
      data    annot/.true./,  dashl/4.0/,  fsizn/10.0/,
     .        fsizt/16.0/
c
c local parameters:
c
c annot =     logical variable with value true iff the plot
c               is to be annotated with the values of wx1,
c               wx2, wy1, and wy2
c cnstr       logical variable used to flag constraint arcs:
c               true iff n0-n1 lies in a constraint region
c dashl =     length (in points, at 72 points per inch) of
c               dashes and spaces in a dashed line pattern
c               used for drawing constraint arcs
c dx =        window width wx2-wx1
c dy =        window height wy2-wy1
c fsizn =     font size in points for labeling nodes with
c               their indexes if numbr = true
c fsizt =     font size in points for the title (and
c               annotation if annot = true)
c i =         constraint index (1 to ncc)
c ifrst =     index of the first node in constraint i
c ih =        height of the viewport in points
c ilast =     index of the last node in constraint i
c ipx1,ipy1 = x and y coordinates (in points) of the lower
c               left corner of the bounding box or viewport
c ipx2,ipy2 = x and y coordinates (in points) of the upper
c               right corner of the bounding box or viewport
c iw =        width of the viewport in points
c lp =        list index (pointer)
c lpl =       pointer to the last neighbor of n0
c n0 =        nodal index and do-loop index
c n0bak =     predecessor of n0 in a constraint curve
c               (sequence of adjacent constraint nodes)
c n0for =     successor to n0 in a constraint curve
c n1 =        index of a neighbor of n0
c nls =       index of the last non-constraint node
c pass1 =     logical variable used to flag the first pass
c               through the constraint nodes
c r =         aspect ratio dx/dy
c sfx,sfy =   scale factors for mapping world coordinates
c               (window coordinates in [wx1,wx2] x [wy1,wy2])
c               to viewport coordinates in [ipx1,ipx2] x
c               [ipy1,ipy2]
c t =         temporary variable
c tx,ty =     translation vector for mapping world coordi-
c               nates to viewport coordinates
c x0,y0 =     x(n0),y(n0) or label location
c
c
c test for error 1, and set nls to the last non-constraint
c   node.
c
      if (lun .lt. 0  .or.  lun .gt. 99  .or.
     .    pltsiz .lt. 1.0  .or.  pltsiz .gt. 8.5  .or.
     .    ncc .lt. 0  .or.  n .lt. 3) go to 11
      nls = n
      if (ncc .gt. 0) nls = lcc(1)-1
c
c compute the aspect ratio of the window.
c
      dx = wx2 - wx1
      dy = wy2 - wy1
      if (dx .le. 0.0  .or.  dy .le. 0.0) go to 12
      r = dx/dy
c
c compute the lower left (ipx1,ipy1) and upper right
c   (ipx2,ipy2) corner coordinates of the bounding box.
c   the coordinates, specified in default user space units
c   (points, at 72 points/inch with origin at the lower
c   left corner of the page), are chosen to preserve the
c   aspect ratio r, and to center the plot on the 8.5 by 11
c   inch page.  the center of the page is (306,396), and
c   t = pltsiz/2 in points.
c
      t = 36.0*pltsiz
      if (r .ge. 1.0) then
        ipx1 = 306 - nint(t)
        ipx2 = 306 + nint(t)
        ipy1 = 396 - nint(t/r)
        ipy2 = 396 + nint(t/r)
      else
        ipx1 = 306 - nint(t*r)
        ipx2 = 306 + nint(t*r)
        ipy1 = 396 - nint(t)
        ipy2 = 396 + nint(t)
      endif
c
c output header comments.
c
      write (lun,100,err=13) ipx1, ipy1, ipx2, ipy2
  100 format ('%!ps-adobe-3.0 epsf-3.0'/
     .        '%%boundingbox:',4i4/
     .        '%%title:  triangulation'/
     .        '%%creator:  tripack'/
     .        '%%endcomments')
c
c set (ipx1,ipy1) and (ipx2,ipy2) to the corner coordinates
c   of a viewport obtained by shrinking the bounding box by
c   12% in each dimension.
c
      iw = nint(0.88*real(ipx2-ipx1))
      ih = nint(0.88*real(ipy2-ipy1))
      ipx1 = 306 - iw/2
      ipx2 = 306 + iw/2
      ipy1 = 396 - ih/2
      ipy2 = 396 + ih/2
c
c set the line thickness to 2 points, and draw the
c   viewport boundary.
c
      t = 2.0
      write (lun,110,err=13) t
      write (lun,120,err=13) ipx1, ipy1
      write (lun,130,err=13) ipx1, ipy2
      write (lun,130,err=13) ipx2, ipy2
      write (lun,130,err=13) ipx2, ipy1
      write (lun,140,err=13)
      write (lun,150,err=13)
  110 format (f12.6,' setlinewidth')
  120 format (2i4,' moveto')
  130 format (2i4,' lineto')
  140 format ('closepath')
  150 format ('stroke')
c
c set up a mapping from the window to the viewport.
c
      sfx = real(iw)/dx
      sfy = real(ih)/dy
      tx = ipx1 - sfx*wx1
      ty = ipy1 - sfy*wy1
      write (lun,160,err=13) tx, ty, sfx, sfy
  160 format (2f12.6,' translate'/
     .        2f12.6,' scale')
c
c the line thickness (believe it or fucking not) must be
c   changed to reflect the new scaling which is applied to
c   all subsequent output.  set it to 1.0 point.
c
      t = 2.0/(sfx+sfy)
      write (lun,110,err=13) t
c
c save the current graphics state, and set the clip path to
c   the boundary of the window.
c
      write (lun,170,err=13)
      write (lun,180,err=13) wx1, wy1
      write (lun,190,err=13) wx2, wy1
      write (lun,190,err=13) wx2, wy2
      write (lun,190,err=13) wx1, wy2
      write (lun,200,err=13)
  170 format ('gsave')
  180 format (2f12.6,' moveto')
  190 format (2f12.6,' lineto')
  200 format ('closepath clip newpath')
c
c draw the edges n0->n1, where n1 > n0, beginning with a
c   loop on non-constraint nodes n0.  lpl points to the
c   last neighbor of n0.
c
      do 3 n0 = 1,nls
        x0 = x(n0)
        y0 = y(n0)
        lpl = lend(n0)
        lp = lpl
c
c   loop on neighbors n1 of n0.
c
    2   lp = lptr(lp)
          n1 = abs(list(lp))
          if (n1 .gt. n0) then
c
c   add the edge to the path.
c
            write (lun,210,err=13) x0, y0, x(n1), y(n1)
  210       format (2f12.6,' moveto',2f12.6,' lineto')
          endif
          if (lp .ne. lpl) go to 2
    3   continue
c
c loop through the constraint nodes twice.  the non-
c   constraint arcs incident on constraint nodes are
c   drawn (with solid lines) on the first pass, and the
c   constraint arcs (both boundary and interior, if any)
c   are drawn (with dashed lines) on the second pass.
c
      pass1 = .true.
c
c loop on constraint nodes n0 with (n0bak,n0,n0for) a sub-
c   sequence of constraint i.  the outer loop is on
c   constraints i with first and last nodes ifrst and ilast.
c
    4 ifrst = n+1
      do 8 i = ncc,1,-1
        ilast = ifrst - 1
        ifrst = lcc(i)
        n0bak = ilast
        do 7 n0 = ifrst,ilast
          n0for = n0 + 1
          if (n0 .eq. ilast) n0for = ifrst
          lpl = lend(n0)
          x0 = x(n0)
          y0 = y(n0)
          lp = lpl
c
c   loop on neighbors n1 of n0.  cnstr = true iff n0-n1 is a
c     constraint arc.
c
c   initialize cnstr to true iff the first neighbor of n0
c     strictly follows n0for and precedes or coincides with
c     n0bak (in counterclockwise order).
c
    5     lp = lptr(lp)
            n1 = abs(list(lp))
            if (n1 .ne. n0for  .and.  n1 .ne. n0bak) go to 5
          cnstr = n1 .eq. n0bak
          lp = lpl
c
c   loop on neighbors n1 of n0.  update cnstr and test for
c     n1 > n0.
c
    6     lp = lptr(lp)
            n1 = abs(list(lp))
            if (n1 .eq. n0for) cnstr = .true.
            if (n1 .gt. n0) then
c
c   draw the edge iff (pass1=true and cnstr=false) or
c     (pass1=false and cnstr=true); i.e., cnstr and pass1
c     have opposite values.
c
              if (cnstr .neqv. pass1)
     .          write (lun,210,err=13) x0, y0, x(n1), y(n1)
            endif
            if (n1 .eq. n0bak) cnstr = .false.
c
c   bottom of loops.
c
            if (lp .ne. lpl) go to 6
          n0bak = n0
    7     continue
    8   continue
      if (pass1) then
c
c end of first pass:  paint the path and change to dashed
c   lines for subsequent drawing.  since the scale factors
c   are applied to everything, the dash length must be
c   specified in world coordinates.
c
        pass1 = .false.
        write (lun,150,err=13)
        t = dashl*2.0/(sfx+sfy)
        write (lun,220,err=13) t
  220   format ('[',f12.6,'] 0 setdash')
        go to 4
      endif
c
c paint the path and restore the saved graphics state (with
c   no clip path).
c
      write (lun,150,err=13)
      write (lun,230,err=13)
  230 format ('grestore')
      if (numbr) then
c
c nodes in the window are to be labeled with their indexes.
c   convert fsizn from points to world coordinates, and
c   output the commands to select a font and scale it.
c
        t = fsizn*2.0/(sfx+sfy)
        write (lun,240,err=13) t
  240   format ('/helvetica findfont'/
     .          f12.6,' scalefont setfont')
c
c   loop on nodes n0 with coordinates (x0,y0).
c
        do 9 n0 = 1,n
          x0 = x(n0)
          y0 = y(n0)
          if (x0 .lt. wx1  .or.  x0 .gt. wx2  .or.
     .        y0 .lt. wy1  .or.  y0 .gt. wy2) go to 9
c
c   move to (x0,y0), and draw the label n0.  the first char-
c     acter will have its lower left corner about one
c     character width to the right of the nodal position.
c
          write (lun,180,err=13) x0, y0
          write (lun,250,err=13) n0
  250     format ('(',i3,') show')
    9     continue
      endif
c
c convert fsizt from points to world coordinates, and output
c   the commands to select a font and scale it.
c
      t = fsizt*2.0/(sfx+sfy)
      write (lun,240,err=13) t
c
c display title centered above the plot:
c
      y0 = wy2 + 3.0*t
      write (lun,260,err=13) title, (wx1+wx2)/2.0, y0
  260 format (a80/'  stringwidth pop 2 div neg ',f12.6,
     .        ' add ',f12.6,' moveto')
      write (lun,270,err=13) title
  270 format (a80/'  show')
      if (annot) then
c
c display the window extrema below the plot.
c
        x0 = wx1
        y0 = wy1 - 100.0/(sfx+sfy)
        write (lun,180,err=13) x0, y0
        write (lun,280,err=13) wx1, wx2
        y0 = y0 - 2.0*t
        write (lun,290,err=13) x0, y0, wy1, wy2
  280   format ('(window:   wx1 = ',e9.3,',   wx2 = ',e9.3,
     .          ') show')
  290   format ('(window:  ) stringwidth pop ',f12.6,' add',
     .          f12.6,' moveto'/
     .          '( wy1 = ',e9.3,',   wy2 = ',e9.3,') show')
      endif
c
c paint the path and output the showpage command and
c   end-of-file indicator.
c
      write (lun,300,err=13)
  300 format ('stroke'/
     .        'showpage'/
     .        '%%eof')
c
c hp's interpreters require a one-byte end-of-postscript-job
c   indicator (to eliminate a timeout error message):
c   ascii 4.
c
      write (lun,310,err=13) char(4)
  310 format (a1)
c
c no error encountered.
c
      ier = 0
      return
c
c invalid input parameter.
c
   11 ier = 1
      return
c
c dx or dy is not positive.
c
   12 ier = 2
      return
c
c error writing to unit lun.
c
   13 ier = 3
      return
      end
      subroutine trprnt (ncc,lcc,n,x,y,list,lptr,lend,lout,
     .                   prntx)
      integer ncc, lcc(*), n, list(*), lptr(*), lend(n),
     .        lout
      logical prntx
      real(8)    x(n), y(n)
c
c***********************************************************
c
c                                               from tripack
c                                            robert j. renka
c                                  dept. of computer science
c                                       univ. of north texas
c                                           renka@cs.unt.edu
c                                                   07/30/98
c
c   given a triangulation of a set of points in the plane,
c this subroutine prints the adjacency lists and, option-
c ally, the nodal coordinates on logical unit lout.  the
c list of neighbors of a boundary node is followed by index
c 0.  the numbers of boundary nodes, triangles, and arcs,
c and the constraint curve starting indexes, if any, are
c also printed.
c
c
c on input:
c
c       ncc = number of constraints.
c
c       lcc = list of constraint curve starting indexes (or
c             dummy array of length 1 if ncc = 0).
c
c       n = number of nodes in the triangulation.
c           3 .le. n .le. 9999.
c
c       x,y = arrays of length n containing the coordinates
c             of the nodes in the triangulation -- not used
c             unless prntx = true.
c
c       list,lptr,lend = data structure defining the trian-
c                        gulation.  refer to subroutine
c                        trmesh.
c
c       lout = logical unit number for output.  0 .le. lout
c              .le. 99.  output is printed on unit 6 if lout
c              is outside its valid range on input.
c
c       prntx = logical variable with value true if and only
c               if x and y are to be printed (to 6 decimal
c               places).
c
c none of the parameters are altered by this routine.
c
c modules required by trprnt:  none
c
c***********************************************************
c
      integer i, inc, k, lp, lpl, lun, na, nabor(100), nb,
     .        nd, nl, nlmax, nmax, node, nn, nt
      data  nmax/9999/,  nlmax/60/
c
      nn = n
      lun = lout
      if (lun .lt. 0  .or.  lun .gt. 99) lun = 6
c
c print a heading and test the range of n.
c
      write (lun,100) nn
      if (nn .lt. 3  .or.  nn .gt. nmax) then
c
c n is outside its valid range.
c
        write (lun,110)
        go to 5
      endif
c
c initialize nl (the number of lines printed on the current
c   page) and nb (the number of boundary nodes encountered).
c
      nl = 6
      nb = 0
      if (.not. prntx) then
c
c print list only.  k is the number of neighbors of node
c   which are stored in nabor.
c
        write (lun,101)
        do 2 node = 1,nn
          lpl = lend(node)
          lp = lpl
          k = 0
c
    1     k = k + 1
            lp = lptr(lp)
            nd = list(lp)
            nabor(k) = nd
            if (lp .ne. lpl) go to 1
          if (nd .le. 0) then
c
c   node is a boundary node.  correct the sign of the last
c     neighbor, add 0 to the end of the list, and increment
c     nb.
c
            nabor(k) = -nd
            k = k + 1
            nabor(k) = 0
            nb = nb + 1
          endif
c
c   increment nl and print the list of neighbors.
c
          inc = (k-1)/14 + 2
          nl = nl + inc
          if (nl .gt. nlmax) then
            write (lun,106)
            nl = inc
          endif
          write (lun,103) node, (nabor(i), i = 1,k)
          if (k .ne. 14) write (lun,105)
    2     continue
      else
c
c print x, y, and list.
c
        write (lun,102)
        do 4 node = 1,nn
          lpl = lend(node)
          lp = lpl
          k = 0
    3     k = k + 1
            lp = lptr(lp)
            nd = list(lp)
            nabor(k) = nd
            if (lp .ne. lpl) go to 3
          if (nd .le. 0) then
c
c   node is a boundary node.
c
            nabor(k) = -nd
            k = k + 1
            nabor(k) = 0
            nb = nb + 1
          endif
c
c   increment nl and print x, y, and nabor.
c
          inc = (k-1)/8 + 2
          nl = nl + inc
          if (nl .gt. nlmax) then
            write (lun,106)
            nl = inc
          endif
          write (lun,104) node, x(node), y(node),
     .                    (nabor(i), i = 1,k)
          if (k .ne. 8) write (lun,105)
    4     continue
      endif
c
c print nb, na, and nt (boundary nodes, arcs, and
c   triangles).
c
      nt = 2*nn - nb - 2
      na = nt + nn - 1
      if (nl .gt. nlmax-6) write (lun,106)
      write (lun,107) nb, na, nt
c
c print ncc and lcc.
c
    5 write (lun,108) ncc
      if (ncc .gt. 0) write (lun,109) (lcc(i), i = 1,ncc)
      return
c
c print formats:
c
  100 format (///,26x,'adjacency sets,    n = ',i5//)
  101 format (1x,'node',32x,'neighbors of node'//)
  102 format (1x,'node',5x,'x(node)',8x,'y(node)',
     .        20x,'neighbors of node'//)
  103 format (1x,i4,5x,14i5/(1x,9x,14i5))
  104 format (1x,i4,2e15.6,5x,8i5/(1x,39x,8i5))
  105 format (1x)
  106 format (///)
  107 format (/1x,'nb = ',i4,' boundary nodes',5x,
     .        'na = ',i5,' arcs',5x,'nt = ',i5,
     .        ' triangles')
  108 format (/1x,'ncc =',i3,' constraint curves')
  109 format (1x,9x,14i5)
  110 format (1x,10x,'*** n is outside its valid',
     .        ' range ***')
      end