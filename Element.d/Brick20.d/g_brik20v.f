C Gradient of stiffness matrix of 20-node hexahedron in 3-D  stress
C
      subroutine  gxbrik20v(x, dx, y, dy, z, dz, c, dc, p, sm,status)
C
C                   A R G U M E N T S
C
      integer           p,status
      double precision  x(20), y(20), z(20), c(6,6)
      double precision  dx(20), dy(20), dz(20), dc(6,6)
      double precision  sm(60,60)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  q(20), qx(20), qy(20), qz(20)
      double precision  dqx(20), dqy(20), dqz(20)
      double precision  xi, eta, mu, det, w, weight, ddet, dw
      double precision  c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z
      double precision  c4x, c4y, c4z, c5x, c5y, c5z, c6x, c6y, c6z
      double precision  dc1x, dc1y, dc1z, dc2x, dc2y, dc2z, dc3x
      double precision  dc3y, dc3z, dc4x, dc4y, dc4z, dc5x, dc5y
      double precision  dc5z, dc6x, dc6y, dc6z, cvar
      integer           i, ix, iy, iz, j, jx, jy, jz, k, l, m
      integer           ls(60)
C
C                   D A T A
C
      data              ls 
     $    /1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,
     $     2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,
     $     3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60 /
C
C                   L O G I C
C
      status =  0
      do 12  j = 1,60
        do 11  i = 1,60
          sm(i,j) = 0.0d0
 11     continue
 12   continue

C     check if any variation
      cvar = 0.0
      
      do 15 i = 1, 8
        cvar = cvar + dx(i)**2 + dy(i)**2 +dz(i)**2 
 15   continue
 	
      if (cvar.eq.0) then
C
      do 3000  k = 1,p
        do 2500  l = 1,p
          do 2400  m = 1,p
            call hxgaus20 (p, k, p, l, p, m, xi, eta, mu, weight)
            call h20shpe (xi,eta,mu,x,y,z,q,qx,qy,qz,det)
            if (det .le. 0.0d0) then
               write(6,*) 'Negative Jacobian determinant in brik20v.f'
               status = -1
               return
            end if
C
            w  =    weight * det
C
            do 2000  j = 1,20
              jx =    ls(j)
              jy =    ls(j+20)
              jz =    ls(j+40)
              c1x = (c(1,1)*qx(j) + c(1,4)*qy(j) + c(1,6)*qz(j))*w
              c1y = (c(1,4)*qx(j) + c(1,2)*qy(j) + c(1,5)*qz(j))*w
              c1z = (c(1,6)*qx(j) + c(1,5)*qy(j) + c(1,3)*qz(j))*w
              c2x = (c(2,1)*qx(j) + c(2,4)*qy(j) + c(2,6)*qz(j))*w
              c2y = (c(2,4)*qx(j) + c(2,2)*qy(j) + c(2,5)*qz(j))*w
              c2z = (c(2,6)*qx(j) + c(2,5)*qy(j) + c(2,3)*qz(j))*w
              c3x = (c(3,1)*qx(j) + c(3,4)*qy(j) + c(3,6)*qz(j))*w
              c3y = (c(3,4)*qx(j) + c(3,2)*qy(j) + c(3,5)*qz(j))*w
              c3z = (c(3,6)*qx(j) + c(3,5)*qy(j) + c(3,3)*qz(j))*w
              c4x = (c(4,1)*qx(j) + c(4,4)*qy(j) + c(4,6)*qz(j))*w
              c4y = (c(4,4)*qx(j) + c(4,2)*qy(j) + c(4,5)*qz(j))*w
              c4z = (c(4,6)*qx(j) + c(4,5)*qy(j) + c(4,3)*qz(j))*w
              c5x = (c(5,1)*qx(j) + c(5,4)*qy(j) + c(5,6)*qz(j))*w
              c5y = (c(5,4)*qx(j) + c(5,2)*qy(j) + c(5,5)*qz(j))*w
              c5z = (c(5,6)*qx(j) + c(5,5)*qy(j) + c(5,3)*qz(j))*w
              c6x = (c(6,1)*qx(j) + c(6,4)*qy(j) + c(6,6)*qz(j))*w
              c6y = (c(6,4)*qx(j) + c(6,2)*qy(j) + c(6,5)*qz(j))*w
              c6z = (c(6,6)*qx(j) + c(6,5)*qy(j) + c(6,3)*qz(j))*w
	      
              dc1x = (dc(1,1)*qx(j) + dc(1,4)*qy(j) + dc(1,6)*qz(j))*w
              dc1y = (dc(1,4)*qx(j) + dc(1,2)*qy(j) + dc(1,5)*qz(j))*w
              dc1z = (dc(1,6)*qx(j) + dc(1,5)*qy(j) + dc(1,3)*qz(j))*w
              dc2x = (dc(2,1)*qx(j) + dc(2,4)*qy(j) + dc(2,6)*qz(j))*w
              dc2y = (dc(2,4)*qx(j) + dc(2,2)*qy(j) + dc(2,5)*qz(j))*w
              dc2z = (dc(2,6)*qx(j) + dc(2,5)*qy(j) + dc(2,3)*qz(j))*w
              dc3x = (dc(3,1)*qx(j) + dc(3,4)*qy(j) + dc(3,6)*qz(j))*w
              dc3y = (dc(3,4)*qx(j) + dc(3,2)*qy(j) + dc(3,5)*qz(j))*w
              dc3z = (dc(3,6)*qx(j) + dc(3,5)*qy(j) + dc(3,3)*qz(j))*w
              dc4x = (dc(4,1)*qx(j) + dc(4,4)*qy(j) + dc(4,6)*qz(j))*w
              dc4y = (dc(4,4)*qx(j) + dc(4,2)*qy(j) + dc(4,5)*qz(j))*w
              dc4z = (dc(4,6)*qx(j) + dc(4,5)*qy(j) + dc(4,3)*qz(j))*w
              dc5x = (dc(5,1)*qx(j) + dc(5,4)*qy(j) + dc(5,6)*qz(j))*w
              dc5y = (dc(5,4)*qx(j) + dc(5,2)*qy(j) + dc(5,5)*qz(j))*w
              dc5z = (dc(5,6)*qx(j) + dc(5,5)*qy(j) + dc(5,3)*qz(j))*w
              dc6x = (dc(6,1)*qx(j) + dc(6,4)*qy(j) + dc(6,6)*qz(j))*w
              dc6y = (dc(6,4)*qx(j) + dc(6,2)*qy(j) + dc(6,5)*qz(j))*w
              dc6z = (dc(6,6)*qx(j) + dc(6,5)*qy(j) + dc(6,3)*qz(j))*w
	      
              do 1500  i = j,20
                ix =     ls(i)
                iy =     ls(i+20)
                iz =     ls(i+40)
                sm(ix,jx) =sm(ix,jx)+qx(i)*dc1x+qy(i)*dc4x +qz(i)*dc6x
                sm(jx,ix) =sm(ix,jx)
                sm(iy,jy) =sm(iy,jy)+qx(i)*dc4y+qy(i)*dc2y +qz(i)*dc5y
                sm(jy,iy) =sm(iy,jy)
                sm(iz,jz) =sm(iz,jz)+qx(i)*dc6z+qy(i)*dc5z +qz(i)*dc3z
                sm(jz,iz) =sm(iz,jz)
                sm(ix,jy) =sm(ix,jy)+qx(i)*dc1y+qy(i)*dc4y +qz(i)*dc6y
                sm(iy,jx) =sm(iy,jx)+qx(i)*dc4x+qy(i)*dc2x +qz(i)*dc5x
                sm(jy,ix) =sm(ix,jy)
                sm(jx,iy) =sm(iy,jx)
                sm(ix,jz) =sm(ix,jz)+qx(i)*dc1z+qy(i)*dc4z +qz(i)*dc6z
                sm(iz,jx) =sm(iz,jx)+qx(i)*dc6x+qy(i)*dc5x +qz(i)*dc3x
                sm(jz,ix) =sm(ix,jz)
                sm(jx,iz) =sm(iz,jx)
                sm(iy,jz) =sm(iy,jz)+qx(i)*dc4z+qy(i)*dc2z +qz(i)*dc5z
                sm(iz,jy) =sm(iz,jy)+qx(i)*dc6y+qy(i)*dc5y +qz(i)*dc3y
                sm(jz,iy) =sm(iy,jz)
                sm(jy,iz) =sm(iz,jy)
 1500           continue
 2000         continue
 2400       continue
 2500     continue
 3000   continue
 
      else
C
      do 6000  k = 1,p
        do 5500  l = 1,p
          do 5400  m = 1,p
            call hxgaus20 (p, k, p, l, p, m, xi, eta, mu, weight)
            call gxh20shpe (xi,eta,mu,x,y,z,dx,dy,dz,q,qx,qy,qz,dqx,
     *                      dqy,dqz,det,ddet)
            if (det .le. 0.0d0) then
               write(6,*) 'Negative Jacobian determinant in brik20v.f'
               status = -1
               return
            end if
C
            w  =    weight * det 
            dw =    weight * ddet
C
            do 5000  j = 1,20
              jx =    ls(j)
              jy =    ls(j+20)
              jz =    ls(j+40)
              c1x = (c(1,1)*qx(j) + c(1,4)*qy(j) + c(1,6)*qz(j))*w
              c1y = (c(1,4)*qx(j) + c(1,2)*qy(j) + c(1,5)*qz(j))*w
              c1z = (c(1,6)*qx(j) + c(1,5)*qy(j) + c(1,3)*qz(j))*w
              c2x = (c(2,1)*qx(j) + c(2,4)*qy(j) + c(2,6)*qz(j))*w
              c2y = (c(2,4)*qx(j) + c(2,2)*qy(j) + c(2,5)*qz(j))*w
              c2z = (c(2,6)*qx(j) + c(2,5)*qy(j) + c(2,3)*qz(j))*w
              c3x = (c(3,1)*qx(j) + c(3,4)*qy(j) + c(3,6)*qz(j))*w
              c3y = (c(3,4)*qx(j) + c(3,2)*qy(j) + c(3,5)*qz(j))*w
              c3z = (c(3,6)*qx(j) + c(3,5)*qy(j) + c(3,3)*qz(j))*w
              c4x = (c(4,1)*qx(j) + c(4,4)*qy(j) + c(4,6)*qz(j))*w
              c4y = (c(4,4)*qx(j) + c(4,2)*qy(j) + c(4,5)*qz(j))*w
              c4z = (c(4,6)*qx(j) + c(4,5)*qy(j) + c(4,3)*qz(j))*w
              c5x = (c(5,1)*qx(j) + c(5,4)*qy(j) + c(5,6)*qz(j))*w
              c5y = (c(5,4)*qx(j) + c(5,2)*qy(j) + c(5,5)*qz(j))*w
              c5z = (c(5,6)*qx(j) + c(5,5)*qy(j) + c(5,3)*qz(j))*w
              c6x = (c(6,1)*qx(j) + c(6,4)*qy(j) + c(6,6)*qz(j))*w
              c6y = (c(6,4)*qx(j) + c(6,2)*qy(j) + c(6,5)*qz(j))*w
              c6z = (c(6,6)*qx(j) + c(6,5)*qy(j) + c(6,3)*qz(j))*w
	      
              dc1x = (dc(1,1)*qx(j) + dc(1,4)*qy(j) + dc(1,6)*qz(j))*w
              dc1y = (dc(1,4)*qx(j) + dc(1,2)*qy(j) + dc(1,5)*qz(j))*w
              dc1z = (dc(1,6)*qx(j) + dc(1,5)*qy(j) + dc(1,3)*qz(j))*w
              dc2x = (dc(2,1)*qx(j) + dc(2,4)*qy(j) + dc(2,6)*qz(j))*w
              dc2y = (dc(2,4)*qx(j) + dc(2,2)*qy(j) + dc(2,5)*qz(j))*w
              dc2z = (dc(2,6)*qx(j) + dc(2,5)*qy(j) + dc(2,3)*qz(j))*w
              dc3x = (dc(3,1)*qx(j) + dc(3,4)*qy(j) + dc(3,6)*qz(j))*w
              dc3y = (dc(3,4)*qx(j) + dc(3,2)*qy(j) + dc(3,5)*qz(j))*w
              dc3z = (dc(3,6)*qx(j) + dc(3,5)*qy(j) + dc(3,3)*qz(j))*w
              dc4x = (dc(4,1)*qx(j) + dc(4,4)*qy(j) + dc(4,6)*qz(j))*w
              dc4y = (dc(4,4)*qx(j) + dc(4,2)*qy(j) + dc(4,5)*qz(j))*w
              dc4z = (dc(4,6)*qx(j) + dc(4,5)*qy(j) + dc(4,3)*qz(j))*w
              dc5x = (dc(5,1)*qx(j) + dc(5,4)*qy(j) + dc(5,6)*qz(j))*w
              dc5y = (dc(5,4)*qx(j) + dc(5,2)*qy(j) + dc(5,5)*qz(j))*w
              dc5z = (dc(5,6)*qx(j) + dc(5,5)*qy(j) + dc(5,3)*qz(j))*w
              dc6x = (dc(6,1)*qx(j) + dc(6,4)*qy(j) + dc(6,6)*qz(j))*w
              dc6y = (dc(6,4)*qx(j) + dc(6,2)*qy(j) + dc(6,5)*qz(j))*w
              dc6z = (dc(6,6)*qx(j) + dc(6,5)*qy(j) + dc(6,3)*qz(j))*w
	      
              dc1x =dc1x+(c(1,1)*dqx(j)+c(1,4)*dqy(j)+c(1,6)*dqz(j))*w
              dc1y =dc1y+(c(1,4)*dqx(j)+c(1,2)*dqy(j)+c(1,5)*dqz(j))*w
              dc1z =dc1z+(c(1,6)*dqx(j)+c(1,5)*dqy(j)+c(1,3)*dqz(j))*w
              dc2x =dc2x+(c(2,1)*dqx(j)+c(2,4)*dqy(j)+c(2,6)*dqz(j))*w
              dc2y =dc2y+(c(2,4)*dqx(j)+c(2,2)*dqy(j)+c(2,5)*dqz(j))*w
              dc2z =dc2z+(c(2,6)*dqx(j)+c(2,5)*dqy(j)+c(2,3)*dqz(j))*w
              dc3x =dc3x+(c(3,1)*dqx(j)+c(3,4)*dqy(j)+c(3,6)*dqz(j))*w
              dc3y =dc3y+(c(3,4)*dqx(j)+c(3,2)*dqy(j)+c(3,5)*dqz(j))*w
              dc3z =dc3z+(c(3,6)*dqx(j)+c(3,5)*dqy(j)+c(3,3)*dqz(j))*w
              dc4x =dc4x+(c(4,1)*dqx(j)+c(4,4)*dqy(j)+c(4,6)*dqz(j))*w
              dc4y =dc4y+(c(4,4)*dqx(j)+c(4,2)*dqy(j)+c(4,5)*dqz(j))*w
              dc4z =dc4z+(c(4,6)*dqx(j)+c(4,5)*dqy(j)+c(4,3)*dqz(j))*w
              dc5x =dc5x+(c(5,1)*dqx(j)+c(5,4)*dqy(j)+c(5,6)*dqz(j))*w
              dc5y =dc5y+(c(5,4)*dqx(j)+c(5,2)*dqy(j)+c(5,5)*dqz(j))*w
              dc5z =dc5z+(c(5,6)*dqx(j)+c(5,5)*dqy(j)+c(5,3)*dqz(j))*w
              dc6x =dc6x+(c(6,1)*dqx(j)+c(6,4)*dqy(j)+c(6,6)*dqz(j))*w
              dc6y =dc6y+(c(6,4)*dqx(j)+c(6,2)*dqy(j)+c(6,5)*dqz(j))*w
              dc6z =dc6z+(c(6,6)*dqx(j)+c(6,5)*dqy(j)+c(6,3)*dqz(j))*w
	      
              dc1x =dc1x+(c(1,1)*qx(j) +c(1,4)*qy(j) +c(1,6)*qz(j))*dw
              dc1y =dc1y+(c(1,4)*qx(j) +c(1,2)*qy(j) +c(1,5)*qz(j))*dw
              dc1z =dc1z+(c(1,6)*qx(j) +c(1,5)*qy(j) +c(1,3)*qz(j))*dw
              dc2x =dc2x+(c(2,1)*qx(j) +c(2,4)*qy(j) +c(2,6)*qz(j))*dw
              dc2y =dc2y+(c(2,4)*qx(j) +c(2,2)*qy(j) +c(2,5)*qz(j))*dw
              dc2z =dc2z+(c(2,6)*qx(j) +c(2,5)*qy(j) +c(2,3)*qz(j))*dw
              dc3x =dc3x+(c(3,1)*qx(j) +c(3,4)*qy(j) +c(3,6)*qz(j))*dw
              dc3y =dc3y+(c(3,4)*qx(j) +c(3,2)*qy(j) +c(3,5)*qz(j))*dw
              dc3z =dc3z+(c(3,6)*qx(j) +c(3,5)*qy(j) +c(3,3)*qz(j))*dw
              dc4x =dc4x+(c(4,1)*qx(j) +c(4,4)*qy(j) +c(4,6)*qz(j))*dw
              dc4y =dc4y+(c(4,4)*qx(j) +c(4,2)*qy(j) +c(4,5)*qz(j))*dw
              dc4z =dc4z+(c(4,6)*qx(j) +c(4,5)*qy(j) +c(4,3)*qz(j))*dw
              dc5x =dc5x+(c(5,1)*qx(j) +c(5,4)*qy(j) +c(5,6)*qz(j))*dw
              dc5y =dc5y+(c(5,4)*qx(j) +c(5,2)*qy(j) +c(5,5)*qz(j))*dw
              dc5z =dc5z+(c(5,6)*qx(j) +c(5,5)*qy(j) +c(5,3)*qz(j))*dw
              dc6x =dc6x+(c(6,1)*qx(j) +c(6,4)*qy(j) +c(6,6)*qz(j))*dw
              dc6y =dc6y+(c(6,4)*qx(j) +c(6,2)*qy(j) +c(6,5)*qz(j))*dw
              dc6z =dc6z+(c(6,6)*qx(j) +c(6,5)*qy(j) +c(6,3)*qz(j))*dw
	      
              do 4500  i = j,20
                ix =     ls(i)
                iy =     ls(i+20)
                iz =     ls(i+40)
                sm(ix,jx) =sm(ix,jx)+qx(i)*dc1x+qy(i)*dc4x +qz(i)*dc6x
                sm(ix,jx) =sm(ix,jx)+dqx(i)*c1x+dqy(i)*c4x +dqz(i)*c6x
                sm(jx,ix) =sm(ix,jx)
                sm(iy,jy) =sm(iy,jy)+qx(i)*dc4y+qy(i)*dc2y +qz(i)*dc5y
                sm(iy,jy) =sm(iy,jy)+dqx(i)*c4y+dqy(i)*c2y +dqz(i)*c5y
                sm(jy,iy) =sm(iy,jy)
                sm(iz,jz) =sm(iz,jz)+qx(i)*dc6z+qy(i)*dc5z +qz(i)*dc3z
                sm(iz,jz) =sm(iz,jz)+dqx(i)*c6z+dqy(i)*c5z +dqz(i)*c3z
                sm(jz,iz) =sm(iz,jz)
                sm(ix,jy) =sm(ix,jy)+qx(i)*dc1y+qy(i)*dc4y +qz(i)*dc6y
                sm(ix,jy) =sm(ix,jy)+dqx(i)*c1y+dqy(i)*c4y +dqz(i)*c6y
                sm(iy,jx) =sm(iy,jx)+qx(i)*dc4x+qy(i)*dc2x +qz(i)*dc5x
                sm(iy,jx) =sm(iy,jx)+dqx(i)*c4x+dqy(i)*c2x +dqz(i)*c5x
                sm(jy,ix) =sm(ix,jy)
                sm(jx,iy) =sm(iy,jx)
                sm(ix,jz) =sm(ix,jz)+qx(i)*dc1z+qy(i)*dc4z +qz(i)*dc6z
                sm(ix,jz) =sm(ix,jz)+dqx(i)*c1z+dqy(i)*c4z +dqz(i)*c6z
                sm(iz,jx) =sm(iz,jx)+qx(i)*dc6x+qy(i)*dc5x +qz(i)*dc3x
                sm(iz,jx) =sm(iz,jx)+dqx(i)*c6x+dqy(i)*c5x +dqz(i)*c3x
                sm(jz,ix) =sm(ix,jz)
                sm(jx,iz) =sm(iz,jx)
                sm(iy,jz) =sm(iy,jz)+qx(i)*dc4z+qy(i)*dc2z +qz(i)*dc5z
                sm(iy,jz) =sm(iy,jz)+dqx(i)*c4z+dqy(i)*c2z +dqz(i)*c5z
                sm(iz,jy) =sm(iz,jy)+qx(i)*dc6y+qy(i)*dc5y +qz(i)*dc3y
                sm(iz,jy) =sm(iz,jy)+dqx(i)*c6y+dqy(i)*c5y +dqz(i)*c3y
                sm(jz,iy) =sm(iy,jz)
                sm(jy,iz) =sm(iz,jy)
 4500           continue
 5000         continue
 5400       continue
 5500     continue
 6000   continue
      endif

C
      return
      end
