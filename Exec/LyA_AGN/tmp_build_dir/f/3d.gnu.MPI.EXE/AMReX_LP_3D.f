

c-----------------------------------------------------------------------
c      
c     Gauss-Seidel Red-Black (GSRB):
c     Apply the GSRB relaxation to the state phi for the equation
c     L(phi) = Div(Grad(phi(x))) = rhs(x) central differenced, according
c     to the arrays of boundary masks (m#) and auxiliary data (f#).
c     
c     In general, if the linear operator L=gamma*y-rho, the GS relaxation
c     is y = (R - rho)/gamma.  Near a boundary, the ghost data is filled
c     using a polynomial interpolant based on the "old" phi values, so
c     L=(gamma-delta)*y - rho + delta*yOld.   The resulting iteration is
c     
c     y = (R - delta*yOld + rho)/(gamma - delta)
c     
c     This expression is valid additionally in the interior provided
c     delta->0 there.  delta is constructed by summing all the
c     contributions to the central stencil element coming from boundary 
c     interpolants.  The f#s contain the corresponding coefficient of 
c     the interpolating polynomial.  The masks are set > 0 if the boundary 
c     value was filled with an interpolant involving the central stencil 
c     element.
c     
c-----------------------------------------------------------------------
      subroutine gsrb3dsim (phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2,
     & phi_h3,rhs, rhs_l1, rhs_l2, rhs_l3, rhs_h1, rhs_h2, rhs_h3,f0, 
     &f0_l1, f0_l2, f0_l3, f0_h1, f0_h2, f0_h3, m0, m0_l1, m0_l2, m0_l
     &3, m0_h1, m0_h2, m0_h3,f1, f1_l1, f1_l2, f1_l3, f1_h1, f1_h2, f1
     &_h3, m1, m1_l1, m1_l2, m1_l3, m1_h1, m1_h2, m1_h3,f2, f2_l1, f2_
     &l2, f2_l3, f2_h1, f2_h2, f2_h3, m2, m2_l1, m2_l2, m2_l3, m2_h1, 
     &m2_h2, m2_h3,f3, f3_l1, f3_l2, f3_l3, f3_h1, f3_h2, f3_h3, m3, m
     &3_l1, m3_l2, m3_l3, m3_h1, m3_h2, m3_h3,f4, f4_l1, f4_l2, f4_l3,
     & f4_h1, f4_h2, f4_h3, m4, m4_l1, m4_l2, m4_l3, m4_h1, m4_h2, m4_
     &h3,f5, f5_l1, f5_l2, f5_l3, f5_h1, f5_h2, f5_h3, m5, m5_l1, m5_l
     &2, m5_l3, m5_h1, m5_h2, m5_h3,lo, hi, blo, bhi,nc, h, redblack)
      implicit none
      integer nc
      integer phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3
      DOUBLE PRECISION phi(phi_l1:phi_h1, phi_l2:phi_h2, phi_l3:phi_h3,n
     &c)
      integer rhs_l1, rhs_l2, rhs_l3, rhs_h1, rhs_h2, rhs_h3
      DOUBLE PRECISION rhs(rhs_l1:rhs_h1, rhs_l2:rhs_h2, rhs_l3:rhs_h3,n
     &c)
      integer lo(3), hi(3)
      integer blo(3), bhi(3)
      integer redblack
      integer f0_l1, f0_l2, f0_l3, f0_h1, f0_h2, f0_h3
      integer f1_l1, f1_l2, f1_l3, f1_h1, f1_h2, f1_h3
      integer f2_l1, f2_l2, f2_l3, f2_h1, f2_h2, f2_h3
      integer f3_l1, f3_l2, f3_l3, f3_h1, f3_h2, f3_h3
      integer f4_l1, f4_l2, f4_l3, f4_h1, f4_h2, f4_h3
      integer f5_l1, f5_l2, f5_l3, f5_h1, f5_h2, f5_h3
      DOUBLE PRECISION  f0(f0_l1:f0_h1, f0_l2:f0_h2, f0_l3:f0_h3)
      DOUBLE PRECISION  f1(f1_l1:f1_h1, f1_l2:f1_h2, f1_l3:f1_h3)
      DOUBLE PRECISION  f2(f2_l1:f2_h1, f2_l2:f2_h2, f2_l3:f2_h3)
      DOUBLE PRECISION  f3(f3_l1:f3_h1, f3_l2:f3_h2, f3_l3:f3_h3)
      DOUBLE PRECISION  f4(f4_l1:f4_h1, f4_l2:f4_h2, f4_l3:f4_h3)
      DOUBLE PRECISION  f5(f5_l1:f5_h1, f5_l2:f5_h2, f5_l3:f5_h3)
      integer m0_l1, m0_l2, m0_l3, m0_h1, m0_h2, m0_h3
      integer m1_l1, m1_l2, m1_l3, m1_h1, m1_h2, m1_h3
      integer m2_l1, m2_l2, m2_l3, m2_h1, m2_h2, m2_h3
      integer m3_l1, m3_l2, m3_l3, m3_h1, m3_h2, m3_h3
      integer m4_l1, m4_l2, m4_l3, m4_h1, m4_h2, m4_h3
      integer m5_l1, m5_l2, m5_l3, m5_h1, m5_h2, m5_h3
      integer m0(m0_l1:m0_h1, m0_l2:m0_h2, m0_l3:m0_h3)
      integer m1(m1_l1:m1_h1, m1_l2:m1_h2, m1_l3:m1_h3)
      integer m2(m2_l1:m2_h1, m2_l2:m2_h2, m2_l3:m2_h3)
      integer m3(m3_l1:m3_h1, m3_l2:m3_h2, m3_l3:m3_h3)
      integer m4(m4_l1:m4_h1, m4_l2:m4_h2, m4_l3:m4_h3)
      integer m5(m5_l1:m5_h1, m5_l2:m5_h2, m5_l3:m5_h3)
      DOUBLE PRECISION  h

      integer  i, j, k, ioff, n

      DOUBLE PRECISION cf0, cf1, cf2, cf3, cf4, cf5
      DOUBLE PRECISION delta, gamma, rho

      parameter(gamma = 6.0D0)

      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               ioff = MOD(lo(1) + j + k + redblack,2)
               do i = lo(1) + ioff,hi(1),2

                  cf0 = merge(f0(blo(1),j,k), 0.0D0,(i .eq. blo(1)) .and
     &. (m0(blo(1)-1,j,k).gt.0))
                  cf1 = merge(f1(i,blo(2),k), 0.0D0,(j .eq. blo(2)) .and
     &. (m1(i,blo(2)-1,k).gt.0))
                  cf2 = merge(f2(i,j,blo(3)), 0.0D0,(k .eq. blo(3)) .and
     &. (m2(i,j,blo(3)-1).gt.0))
                  cf3 = merge(f3(bhi(1),j,k), 0.0D0,(i .eq. bhi(1)) .and
     &. (m3(bhi(1)+1,j,k).gt.0))
                  cf4 = merge(f4(i,bhi(2),k), 0.0D0,(j .eq. bhi(2)) .and
     &. (m4(i,bhi(2)+1,k).gt.0))
                  cf5 = merge(f5(i,j,bhi(3)), 0.0D0,(k .eq. bhi(3)) .and
     &. (m5(i,j,bhi(3)+1).gt.0))

                  delta = cf0 + cf1 + cf2 + cf3 + cf4 + cf5

                  rho = phi(i-1,j,k,n) + phi(i+1,j,k,n)+ phi(i,j-1,k,n) 
     &+ phi(i,j+1,k,n)+ phi(i,j,k-1,n) + phi(i,j,k+1,n)

                  phi(i,j,k,n) = (rhs(i,j,k,n)*h*h - rho + phi(i,j,k,n)*
     &delta)/ (delta - gamma)

               end do
            end do
         end do
      end do

      end
c-----------------------------------------------------------------------
c
c     Fill in a matrix x vector operator here
c
      subroutine adotx3dsim(y, y_l1, y_l2, y_l3, y_h1, y_h2, y_h3,x, x_l
     &1, x_l2, x_l3, x_h1, x_h2, x_h3,lo, hi, nc,h)
      implicit none
      integer nc
      integer lo(3)
      integer hi(3)
      integer y_l1, y_l2, y_l3, y_h1, y_h2, y_h3
      DOUBLE PRECISION  y(y_l1:y_h1, y_l2:y_h2, y_l3:y_h3,nc)
      integer x_l1, x_l2, x_l3, x_h1, x_h2, x_h3
      DOUBLE PRECISION  x(x_l1:x_h1, x_l2:x_h2, x_l3:x_h3,nc)
      DOUBLE PRECISION  h

      integer i, j, k, n
      DOUBLE PRECISION scal

      scal = 1.0D0/h**2

      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  y(i,j,k,n) = scal*( x(i-1,j,k,n) + x(i+1,j,k,n)+ x(i,j
     &-1,k,n) + x(i,j+1,k,n)+ x(i,j,k-1,n) + x(i,j,k+1,n)-
     & 6*x(i,j,k,n) )
               end do
            end do
         end do
      end do

      end

c-----------------------------------------------------------------------
c
c     Fill in fluxes
c
      subroutine flux3dsim(x,x_l1, x_l2, x_l3, x_h1, x_h2, x_h3,xlo,xhi,
     &ylo,yhi,zlo,zhi,nc,h,xflux,xflux_l1, xflux_l2, xflux_l3, xflux_h
     &1, xflux_h2, xflux_h3,yflux,yflux_l1, yflux_l2, yflux_l3, yflux_
     &h1, yflux_h2, yflux_h3,zflux,zflux_l1, zflux_l2, zflux_l3, zflux
     &_h1, zflux_h2, zflux_h3)
      implicit none
      integer xlo(3), xhi(3)
      integer ylo(3), yhi(3)
      integer zlo(3), zhi(3)
      integer nc
      integer x_l1, x_l2, x_l3, x_h1, x_h2, x_h3
      integer xflux_l1, xflux_l2, xflux_l3, xflux_h1, xflux_h2, xflux_h3
      integer yflux_l1, yflux_l2, yflux_l3, yflux_h1, yflux_h2, yflux_h3
      integer zflux_l1, zflux_l2, zflux_l3, zflux_h1, zflux_h2, zflux_h3
      DOUBLE PRECISION  x(x_l1:x_h1, x_l2:x_h2, x_l3:x_h3,nc)
      DOUBLE PRECISION xflux(xflux_l1:xflux_h1, xflux_l2:xflux_h2, xflux
     &_l3:xflux_h3,nc)
      DOUBLE PRECISION yflux(yflux_l1:yflux_h1, yflux_l2:yflux_h2, yflux
     &_l3:yflux_h3,nc)
      DOUBLE PRECISION zflux(zflux_l1:zflux_h1, zflux_l2:zflux_h2, zflux
     &_l3:zflux_h3,nc)
      DOUBLE PRECISION h(3)

      DOUBLE PRECISION dhx, dhy, dhz
      integer i,j,k,n

      dhx = 1.0D0/h(1)
      dhy = 1.0D0/h(2)
      dhz = 1.0D0/h(3)

      do n = 1, nc
         do k = xlo(3), xhi(3)
            do j = xlo(2), xhi(2)
               do i = xlo(1), xhi(1)
                  xflux(i,j,k,n) = - dhx*( x(i,j,k,n) - x(i-1,j,k,n) )
               end do
            end do
         end do

         do k = ylo(3), yhi(3)
            do j = ylo(2), yhi(2)
               do i = ylo(1), yhi(1)
                  yflux(i,j,k,n) = - dhy*( x(i,j,k,n) - x(i,j-1,k,n) )
               end do
            end do
         end do

         do k = zlo(3), zhi(3)
            do j = zlo(2), zhi(2)
               do i = zlo(1), zhi(1)
                  zflux(i,j,k,n) = - dhz*( x(i,j,k,n) - x(i,j,k-1,n) )
               end do
            end do
         end do
      end do

      end
