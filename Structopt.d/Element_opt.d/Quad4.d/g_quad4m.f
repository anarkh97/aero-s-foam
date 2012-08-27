      subroutine gxquad4m(x, g_x, y, g_y, h, g_h, c, g_c, p, 
     &                    sm, g_sm, m, rip)
C
C                   A R G U M E N T S
C
        integer p, m
        double precision x(4), y(4), h(4), c(3,3)
        double precision sm(8,8)
C
	double precision g_x(4), g_y(4), g_h(4), g_c(3,3)
        double precision g_sm(8,8)
C
        double precision rip
C
C                   L O C A L   V A R I A B L E S
C
        double precision q(4), qx(4), qy(4)
        double precision xi, eta, det, w, weight
        double precision c1x, c1y, c2x, c2y, c3x, c3y
        integer i, ix, iy, j, jx, jy, k, l
        integer ls(8), ids(3)
C
C                   D A T A
C
        double precision g_q(4), g_qx(4), g_qy(4)
	double precision g_c1x, g_c2x, g_c3x, g_c1y, g_c2y, g_c3y
	double precision g_w, g_det

        double precision d14_b, d13_b, d12_b, d7_v, d13_v, d9_b, d8_b 
	double precision d2_v, d7_b, d6_b

        data ls  /1, 3, 5, 7, 2, 4, 6, 8/
        data ids /1,2,4/
C
        do j = 1, 8
          do i = 1, 8
              sm(i, j) = 0.0d0
            g_sm(i, j) = 0.0d0
          enddo
	enddo
c
c       Gauss quadrature loop
C
        do 99994 k = 1, p
          do 99995 l = 1, p
c	  
            call qgauss(p, k, p, l, xi, eta, weight)
            call gxq4shpe(xi, eta, x, g_x, y, g_y, q ,qx, g_qx,  
     &                    qy, g_qy, det, g_det)
     
            if (det .le. 0.0) then
              write (6,*) 'Negative Jacobian determinant in 4 node quad'

              if (det .eq. 0.0) then
                write (6, *) 'Zero Jacobian determinant in 4 node quad'
              endif
              stop
            endif

            d2_v  = weight * det
            d13_v = h(1)*q(1) + h(2)*q(2) + h(3)*q(3) + h(4)*q(4)
            d6_b  = d2_v * q(4)
            d9_b  = d2_v * q(3)
            d12_b = d2_v * q(2)
            d13_b = d2_v * q(1)
            d14_b = d13_v * weight

            g_w   = d6_b * g_h(4) + d9_b * g_h(3) + d12_b * g_h(2) 
     &            + d13_b * g_h(1) + d14_b * g_det

            w = d2_v * d13_v

C
            do 99996 j = 1, 4

              jx = ls(j)
              jy = ls(j + 4)
	      
              d7_v = c(1, 1) * qx(j) + c(1, 3) * qy(j)
              d6_b = w * qy(j)
              d7_b = w * c(1, 3)
              d8_b = w * qx(j)
              d9_b = w * c(1, 1)

              g_c1x = d7_v * g_w + d7_b * g_qy(j) + d6_b * g_c(1, 3)
     &              + d9_b * g_qx(j) + d8_b * g_c(1,1)

              c1x = d7_v * w

              d7_v = c(1, 3) * qx(j) + c(1, 2) * qy(j)
              d6_b = w * qy(j)
              d7_b = w * c(1, 2)
              d8_b = w * qx(j)
              d9_b = w * c(1, 3)

              g_c1y = d7_v * g_w + d7_b * g_qy(j) + d6_b * g_c(1, 2) 
     &              + d9_b * g_qx(j) + d8_b * g_c(1, 3)

              c1y = d7_v * w

              d7_v = c(1, 2) * qx(j) + c(2, 3) * qy(j)
              d6_b = w * qy(j)
              d7_b = w * c(2, 3)
              d8_b = w * qx(j)
              d9_b = w * c(1, 2)

              g_c2x = d7_v * g_w + d7_b * g_qy(j) + d6_b * g_c(2, 3) 
     &              + d9_b * g_qx(j) + d8_b * g_c(1, 2)

              c2x = d7_v * w

              d7_v = c(2, 3) * qx(j) + c(2, 2) * qy(j)
              d6_b = w * qy(j)
              d7_b = w * c(2, 2)
              d8_b = w * qx(j)
              d9_b = w * c(2, 3)

              g_c2y = d7_v * g_w + d7_b * g_qy(j) + d6_b * g_c(2, 2) 
     &              + d9_b * g_qx(j) + d8_b * g_c(2, 3)

              c2y = d7_v * w

              d7_v = c(1, 3) * qx(j) + c(3, 3) * qy(j)
              d6_b = w * qy(j)
              d7_b = w * c(3, 3)
              d8_b = w * qx(j)
              d9_b = w * c(1, 3)

              g_c3x = d7_v * g_w + d7_b * g_qy(j) + d6_b * g_c(3, 3) 
     &              + d9_b * g_qx(j) + d8_b * g_c(1, 3)

              c3x = d7_v * w

              d7_v = c(3, 3) * qx(j) + c(2, 3) * qy(j)
              d6_b = w * qy(j)
              d7_b = w * c(2, 3)
              d8_b = w * qx(j)
              d9_b = w * c(3, 3)

              g_c3y = d7_v * g_w + d7_b * g_qy(j) + d6_b * g_c(2, 3) 
     &              + d9_b * g_qx(j) + d8_b * g_c(3, 3)

              c3y = d7_v * w

              do 99997 i = j, 4

                ix = ls(i)
                iy = ls(i + 4)

                g_sm(ix, jx) = qy(i) * g_c3x + c3x * g_qy(i) 
     &                       + qx(i) * g_c1x + c1x * g_qx(i) 
     &                       + g_sm(ix, jx)

                sm(ix, jx) = sm(ix, jx) + qx(i) * c1x + qy(i) * c3x

                g_sm(jx, ix) = g_sm(ix, jx)

                sm(jx, ix) = sm(ix, jx)

                g_sm(iy, jy) = qy(i) * g_c2y + c2y * g_qy(i) 
     &                       + qx(i) * g_c3y + c3y * g_qx(i) 
     &                       + g_sm(iy, jy)

                sm(iy, jy) = sm(iy, jy) + qx(i) * c3y + qy(i) * c2y

                g_sm(jy, iy) = g_sm(iy, jy)

                sm(jy, iy) = sm(iy, jy)

                g_sm(ix, jy) = qy(i) * g_c3y + c3y * g_qy(i) 
     &                       + qx(i) * g_c1y + c1y * g_qx(i) 
     &                       + g_sm(ix, jy)

                sm(ix, jy) = sm(ix, jy) + qx(i) * c1y + qy(i) * c3y

                g_sm(iy, jx) = qy(i) * g_c2x + c2x * g_qy(i) 
     &                       + qx(i) * g_c3x + c3x * g_qx(i) 
     &                       + g_sm(iy, jx)
     
                sm(iy, jx) = sm(iy, jx) + qx(i) * c3x + qy(i) * c2x

                g_sm(jy, ix) = g_sm(ix, jy)

                sm(jy, ix) = sm(ix, jy)

                g_sm(jx, iy) = g_sm(iy, jx)

                sm(jx, iy) = sm(iy, jx)

99997         continue
99996       continue
99995     continue
99994   continue
C
        return
      end
