      subroutine     gxbrkcmt(e, de, nu, dnu, c, dc)
      
      double precision e, nu, c(6,6)
      double precision de, dnu, dc(6,6)
      double precision c1, c2
      integer i, j
c
      do 20  i=1,6
        do 10  j=1,6
           c(i,j) = 0.0d0
          dc(i,j) = 0.0d0
 10     continue
 20   continue
C 
      c1 = 1.0d0/(1.0d0-2.0d0*nu)
      c2 = 1.0d0/(1.0d0+nu)
C
      c(1,1) =  e * (1.0d0-nu)*c1*c2
      c(2,2) =  c(1,1)
      c(3,3) =  c(1,1)
      c(1,2) =  c(1,1)*nu/(1.0d0-nu)
      c(2,1) =  c(1,2)
      c(1,3) =  c(1,2)
      c(2,3) =  c(1,2)
      c(3,1) =  c(1,2)
      c(3,2) =  c(1,2)
      c(4,4) =  e/2.0d0*c2
      c(5,5) =  c(4,4)
      c(6,6) =  c(4,4)
c
      dc(1,1) =  c1*c2*(de*(1.0d0-nu)-e*dnu
     $           +e*(1.0d0-nu)*(1.0d0+4.0d0*nu)*dnu*c1*c2)
      dc(2,2) =  dc(1,1)
      dc(3,3) =  dc(1,1)
      dc(1,2) =  c1*c2*(de*nu+e*dnu+e*nu*(1.0d0+4.0d0*nu)*dnu*c1*c2)
      dc(2,1) =  dc(1,2)
      dc(1,3) =  dc(1,2)
      dc(2,3) =  dc(1,2)
      dc(3,1) =  dc(1,2)
      dc(3,2) =  dc(1,2)
      dc(4,4) =  c2*(de-e*dnu*c2)/2.0d0
      dc(5,5) =  dc(4,4)
      dc(6,6) =  dc(4,4)
c
      return
      end
C=END FORTRAN
